"""
Comparative Genomics Pipeline

A bioinformatics pipeline for comparative genomics analysis focusing on 
epilepsy-related genes across multiple species.
"""

__version__ = "0.1.0"
__author__ = "William P""
__email__ = "author@gmail.com"

from .client import UniProtClient, NCBIClient, EBIClient
from .service import biopython_service
from .config import logging_config, path_config
from .util import file_util

__all__ = [
    "UniProtClient",
    "NCBIClient", 
    "EBIClient",
    "biopython_service",
    "logging_config",
    "path_config",
    "file_util",
]