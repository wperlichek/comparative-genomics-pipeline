from pathlib import Path
from Bio import Phylo
import matplotlib.pyplot as plt
from ..config import path_config

def visualize_and_save_trees(tree_files=None, output_dir=None):
    """
    Visualize each Newick tree in the list and save as PNG to the trees output folder.
    Args:
        tree_files (list[Path] or None): List of .nwk tree file paths. If None, will glob all in TREES_OUTPUT_DIR.
        output_dir (Path or None): Where to save PNGs. If None, uses TREES_OUTPUT_DIR.
    """
    if output_dir is None:
        output_dir = path_config.TREES_OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    if tree_files is None:
        tree_files = list(output_dir.glob("*.nwk"))

    for tree_path in tree_files:
        tree = Phylo.read(tree_path, "newick")
        fig = plt.figure(figsize=(10, 5))
        Phylo.draw(tree, do_show=False)
        png_path = output_dir / f"{tree_path.stem}.png"
        plt.savefig(png_path, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved visualization for {tree_path.name} to {png_path}")
