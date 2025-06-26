import asyncio
import logging
from .client import UniProtClient, NCBIClient, EBIClient
from .config import logging_config, path_config
from .util import file_util
from .service import biopython_service
from pathlib import Path

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
        # Use direct path reference for orthologs
        orthologs_dir = path_config.ORTHOLOGS_OUTPUT_DIR
        orthologs_dir.mkdir(parents=True, exist_ok=True)
        file_path = orthologs_dir / f"{gene_name}.fasta"
        with open(file_path, "w") as fasta_file:
            fasta_file.write(all_orthologs_as_fasta)


async def align_sequences_msa(ebi_client):
    # Use direct path reference for orthologs and msa
    orthologs_dir = path_config.ORTHOLOGS_OUTPUT_DIR
    msa_dir = path_config.MSA_OUTPUT_DIR
    msa_dir.mkdir(parents=True, exist_ok=True)
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
                msa_path = msa_dir / f"{gene_name}.fasta"
                with open(msa_path, "w") as msa_file:
                    msa_file.write(result)
            elif status in ("ERROR", "FAILURE"):
                print(f"Job {job_id} failed with status: {status}")
                break
            else:
                await asyncio.sleep(5)  # wait before polling again


async def generate_phylogenetic_trees(ebi_client):
    # Use direct path reference for orthologs and trees
    orthologs_dir = path_config.ORTHOLOGS_OUTPUT_DIR
    trees_dir = path_config.TREES_OUTPUT_DIR
    trees_dir.mkdir(parents=True, exist_ok=True)
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
            if tree:
                gene_name = file_path.stem
                tree_path = trees_dir / f"{gene_name}.nwk"
                with open(tree_path, "w") as tree_file:
                    tree_file.write(tree)
                print(f"Tree for {file_path.name} saved to {tree_path}")
            else:
                print(f"No tree returned for {file_path.name}")
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

    biopython_service.compute_conservation_for_all_msas()
    biopython_service.plot_all_conservation_scores()
    biopython_service.visualize_and_save_trees()

    # Fetch and save variants for P35498
    accession = "P35498"
    output_dir = str(path_config.VARIANTS_OUTPUT_DIR)
    logger.info(f"Fetching variants for {accession}...")
    await uni_prot_client.fetch_protein_variants_by_accession_id(accession, output_dir)
    # Overlay variants on conservation plot
    conservation_csv = str(Path(path_config.CONSERVATION_OUTPUT_DIR) / "P35498_conservation.csv")
    variants_csv = str(Path(output_dir) / f"{accession}_variants.csv")
    biopython_service.plot_variants_on_conservation(conservation_csv, variants_csv, output_dir)
    await uni_prot_client.close()
    pass


def main():
    logging_config.setup_logging()
    logger.info("Starting pipeline...")
    asyncio.run(async_main())


if __name__ == "__main__":
    main()
