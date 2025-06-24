import asyncio, logging
from .config import logging_config
from .client import UniProtClient

logger = logging.getLogger(__name__)

async def async_main():
    accession_id = "G1SSP8"
    uni_prot_client = UniProtClient()
    fasta_sequence = await uni_prot_client.fetch_protein_fasta_sequence_by_accession_id(accession_id)
    print(f"FASTA sequence for {accession_id}:\n{fasta_sequence}")
    # TODO ::    One more organism in data/input/genes_to_proteins.json needed for MSA
    #            Need client for entrez_protein_id (Biopython could work)
    #            Need client for MSA (European bioinformatics Clustal Omega REST API)
    pass


def main():
    logging_config.setup_logging()
    logger.info("Starting pipeline...")
    asyncio.run(async_main())


if __name__ == "__main__":
    main()
