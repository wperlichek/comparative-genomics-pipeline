import asyncio, logging, json
from .config import logging_config, path_config
from .client import UniProtClient, NCBIClient

logger = logging.getLogger(__name__)

async def async_main():
    accession_id = "G1SSP8"
    entrez_protein_id = "XP_051685089.1"
    uni_prot_client = UniProtClient()
    ncbi_client = NCBIClient()

    print("Data input path: " + f"{path_config.DATA_INPUT_DIR}/genes_to_proteins.json")
    
    with open(f"{path_config.DATA_INPUT_DIR}/genes_to_proteins.json", "r") as f:
        data = json.load(f)
        print("Got data: " + str(data))

    fasta_sequence_accession = await uni_prot_client.fetch_protein_fasta_sequence_by_accession_id(accession_id)
    fasta_sequence_entrez = await ncbi_client.fetch_protein_fasta_by_entrez_id(entrez_protein_id)
    print(f"FASTA sequence for {accession_id}:\n{fasta_sequence_accession}")
    print(f"FASTA sequence for {entrez_protein_id}:\n{fasta_sequence_entrez}")
    pass


def main():
    logging_config.setup_logging()
    logger.info("Starting pipeline...")
    asyncio.run(async_main())


if __name__ == "__main__":
    main()
