import json
from ..config import path_config


def open_file_return_as_json(file):
    ## TODO :: try/except protection
    ## TODO :: Use pathlib.Path to build paths (safer and more cross-platform)?
    with open(file, "r") as f:
        contents = json.load(f)
        return contents


def save_fasta_to_output_dir(file_name: str = "", fasta: str = ""):
    ## TODO :: try/except protection
    with open(f"{path_config.DATA_OUTPUT_DIR}/{file_name}.fasta", "w") as fasta_file:
        fasta_file.write(fasta)
