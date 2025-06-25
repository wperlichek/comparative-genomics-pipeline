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
    Outputs both entropy with and without gaps.
    """
    alignment = AlignIO.read(str(msa_file), "fasta")
    aln_len = alignment.get_alignment_length()
    scores = []
    for i in range(aln_len):
        column_full = str(alignment[:, i])  # includes gaps
        column_nogap = column_full.replace("-", "")  # excludes gaps

        # With gaps
        total_full = len(column_full)
        if total_full == 0:
            entropy_full = 0.0
        else:
            freqs_full = [column_full.count(aa) / total_full for aa in set(column_full)]
            entropy_full = -sum(p * np.log2(p) for p in freqs_full if p > 0)
            entropy_full = abs(entropy_full)

        # Without gaps
        total_nogap = len(column_nogap)
        if total_nogap == 0:
            entropy_nogap = 0.0
            consensus = "-"
        else:
            freqs_nogap = [
                column_nogap.count(aa) / total_nogap for aa in set(column_nogap)
            ]
            entropy_nogap = -sum(p * np.log2(p) for p in freqs_nogap if p > 0)
            entropy_nogap = abs(entropy_nogap)
            consensus = max(set(column_nogap), key=column_nogap.count)

        scores.append((i + 1, entropy_full, entropy_nogap, consensus))

    if output_file is None:
        output_dir = path_config.CONSERVATION_OUTPUT_DIR
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / f"{msa_file.stem}_conservation.csv"
    with open(output_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(
            [
                "Position",
                "ShannonEntropy_WithGaps",
                "ShannonEntropy_NoGaps",
                "ConsensusResidue",
            ]
        )
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
