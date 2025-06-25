import asyncio, logging, json
from pathlib import Path
from .config import logging_config, path_config
from .client import UniProtClient, NCBIClient, EBIClient
from .util import file_util

logger = logging.getLogger(__name__)


async def collect_orthologous_sequences(uni_prot_client, ncbi_client):
    genes_to_proteins = file_util.open_file_return_as_json(
        f"{path_config.DATA_INPUT_DIR}/genes_to_proteins.json"
    )
    for gene_name, ortholog_list in genes_to_proteins.items():
        all_orthologs_as_fasta = ""
        for ortholog in ortholog_list:
            if "uniprot_id" in ortholog and ortholog["uniprot_id"] != "":
                all_orthologs_as_fasta += (
                    await uni_prot_client.fetch_protein_fasta_sequence_by_accession_id(
                        ortholog["uniprot_id"]
                    )
                )
            elif (
                "entrez_protein_id" in ortholog and ortholog["entrez_protein_id"] != ""
            ):
                all_orthologs_as_fasta += (
                    await ncbi_client.fetch_protein_fasta_by_entrez_id(
                        ortholog["entrez_protein_id"]
                    )
                )
        file_util.save_fasta_to_output_dir(
            gene_name, "orthologs", all_orthologs_as_fasta
        )


async def align_sequences_msa(ebi_client):
    orthologs_dir = Path(path_config.DATA_OUTPUT_DIR) / "orthologs"
    for file_path in orthologs_dir.glob("*.fasta"):
        ortholog_fasta = file_util.open_file_return_as_str(file_path)
        job_id = await ebi_client.submit_job(ortholog_fasta)
        gene_name = file_path.stem
        if not job_id:
            return
        status = ""
        while status not in ("FINISHED", "ERROR", "FAILURE"):
            status = await ebi_client.check_status(job_id)
            if status == "FINISHED":
                result = await ebi_client.get_result(job_id, "fa")
                file_util.save_fasta_to_output_dir(gene_name, "msa", result)
            elif status in ("ERROR", "FAILURE"):
                print(f"Job {job_id} failed with status: {status}")
                break
            else:
                await asyncio.sleep(5)  # wait before polling again


async def generate_phylogenetic_trees(ebi_client):
    orthologs_dir = Path(path_config.DATA_OUTPUT_DIR) / "orthologs"
    for file_path in orthologs_dir.glob("*.fasta"):
        with open(file_path, "r") as f:
            fasta_str = f.read()
        print(f"Submitting {file_path.name} for tree generation...")
        job_id = await ebi_client.submit_job(fasta_str)
        if not job_id:
            print(f"Submission failed for {file_path.name}")
            continue
        await asyncio.sleep(15)  # Wait for job to (hopefully) finish
        status = await ebi_client.check_status(job_id)
        print(f"Job {job_id} status: {status}")
        if status == "FINISHED":
            tree = await ebi_client.get_phylogenetic_tree(job_id)
            print(
                f"Tree for {file_path.name} (first 200 chars):\n{tree[:200] if tree else 'No tree returned'}\n"
            )
        else:
            print(f"Job {job_id} did not finish successfully.")


async def async_main():
    uni_prot_client = UniProtClient()
    ncbi_client = NCBIClient()
    ebi_client = EBIClient()

    await collect_orthologous_sequences(uni_prot_client, ncbi_client)
    await align_sequences_msa(ebi_client)
    await generate_phylogenetic_trees(ebi_client)
    await ebi_client.close()

    # Add more modular steps here as needed
    pass


def main():
    logging_config.setup_logging()
    logger.info("Starting pipeline...")
    asyncio.run(async_main())


if __name__ == "__main__":
    main()
