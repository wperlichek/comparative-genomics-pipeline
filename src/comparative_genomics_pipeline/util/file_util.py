import json
from ..config import path_config


def open_file_return_as_json(file):
    ## TODO :: try/except protection
    ## TODO :: Use pathlib.Path to build paths (safer and more cross-platform)?
    with open(file, "r") as f:
        contents = json.load(f)
        return contents


def open_file_return_as_str(file):
    with open(file, "r") as f:
        return f.read()


