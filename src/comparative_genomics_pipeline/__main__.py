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

    protein_ids_and_needed_sources = ""  # {id: source}

    # Iterate and make one or other call, concatenating a new fasta file and save it?
    fasta_sequence_accession = (
        await uni_prot_client.fetch_protein_fasta_sequence_by_accession_id("")
    )
    fasta_sequence_entrez = await ncbi_client.fetch_protein_fasta_by_entrez_id("")

    pass


def main():
    logging_config.setup_logging()
    logger.info("Starting pipeline...")
    asyncio.run(async_main())


if __name__ == "__main__":
    main()
