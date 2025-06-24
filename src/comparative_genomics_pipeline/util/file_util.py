import json
from ..config import path_config


def open_file_return_as_json(file):
    ## TODO :: try/except protection
    with open(file, "r") as f:
        contents = json.load(f)
        return contents
