# __init__.py

from .uniprot_client import UniProtClient
from .ncbi_client import NCBIClient

__all__ = ["UniProtClient", "NCBIClient"]