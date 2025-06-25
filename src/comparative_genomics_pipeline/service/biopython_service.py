from pathlib import Path
from Bio import Phylo, AlignIO
from Bio.Align import AlignInfo
import matplotlib.pyplot as plt
import numpy as np
import csv
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


def compute_conservation_scores(msa_file, output_file=None):
    """
    Compute conservation (Shannon entropy) for each column in the MSA and save as CSV.
    Uses Biopython for alignment parsing and consensus extraction.
    Args:
        msa_file (Path): Path to the MSA FASTA file.
        output_file (Path or None): Where to save the CSV. If None, saves to CONSERVATION_OUTPUT_DIR.
    Returns:
        Path to the output CSV file.
    """
    alignment = AlignIO.read(str(msa_file), "fasta")
    aln_len = alignment.get_alignment_length()
    scores = []
    for i in range(aln_len):
        column = str(alignment[:, i]).replace("-", "")  # remove gaps
        total = len(column)
        if total == 0:
            entropy = 0.0
            consensus = "-"
        else:
            freqs = [column.count(aa) / total for aa in set(column)]
            entropy = -sum(p * np.log2(p) for p in freqs if p > 0)
            entropy = abs(entropy)  # Ensure non-negative
            consensus = max(set(column), key=column.count)
        scores.append((i + 1, entropy, consensus))
    if output_file is None:
        output_dir = path_config.CONSERVATION_OUTPUT_DIR
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / f"{msa_file.stem}_conservation.csv"
    with open(output_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Position", "ShannonEntropy", "ConsensusResidue"])
        writer.writerows(scores)
    print(f"Saved conservation scores to {output_file}")
    return output_file


def compute_conservation_for_all_msas(msa_dir=None):
    """
    Compute conservation scores for all MSA files in the given directory.
    """
    if msa_dir is None:
        msa_dir = path_config.MSA_OUTPUT_DIR
    for msa_file in Path(msa_dir).glob("*.fasta"):
        compute_conservation_scores(msa_file)
