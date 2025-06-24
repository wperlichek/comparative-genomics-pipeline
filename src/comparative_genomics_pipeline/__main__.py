import asyncio, logging, json
from .config import logging_config, path_config
from .client import UniProtClient, NCBIClient
from .util import file_util

logger = logging.getLogger(__name__)


async def async_main():

    uni_prot_client = UniProtClient()
    ncbi_client = NCBIClient()

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
        file_util.save_fasta_to_output_dir(gene_name, all_orthologs_as_fasta)
    pass


def main():
    logging_config.setup_logging()
    logger.info("Starting pipeline...")
    asyncio.run(async_main())


if __name__ == "__main__":
    main()
