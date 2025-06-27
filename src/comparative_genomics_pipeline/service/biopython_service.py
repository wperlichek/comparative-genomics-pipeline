from pathlib import Path
from Bio import Phylo, AlignIO
from Bio.Align import AlignInfo
import matplotlib.pyplot as plt
import numpy as np
import csv
import pandas as pd
from ..config import path_config
from ..visualization import ConservationPlotter, PhylogeneticPlotter, VariantPlotter


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
    print(
        f"DEBUG: Alignment has {len(alignment)} sequences, length {alignment.get_alignment_length()}"
    )
    # Print first 5 columns for inspection
    for i in range(5):
        column_full = str(alignment[:, i])
        column_nogap = column_full.replace("-", "")
        print(f"Col {i+1}: full='{column_full}', nogap='{column_nogap}'")
        print(f"  unique (full): {set(column_full)}")
        print(f"  unique (nogap): {set(column_nogap)}")
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


def plot_conservation_scores(csv_file, output_dir=None):
    """
    Plot Shannon entropy (with and without gaps) from a conservation CSV file.
    Saves the plot as a PNG in the output_dir (defaults to CONSERVATION_OUTPUT_DIR).
    """
    if output_dir is None:
        output_dir = path_config.CONSERVATION_OUTPUT_DIR
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(csv_file)
    plt.figure(figsize=(12, 5))
    plt.plot(
        df["Position"], df["ShannonEntropy_WithGaps"], label="With Gaps", alpha=0.7
    )
    plt.plot(df["Position"], df["ShannonEntropy_NoGaps"], label="No Gaps", alpha=0.7)
    plt.xlabel("Alignment Position")
    plt.ylabel("Shannon Entropy")
    plt.title(f"Conservation (Shannon Entropy): {csv_file.stem}")
    plt.legend()
    plt.tight_layout()
    png_path = output_dir / f"{csv_file.stem}_entropy.png"
    plt.savefig(png_path)
    plt.close()
    print(f"Saved conservation plot to {png_path}")


def plot_all_conservation_scores(conservation_dir=None):
    """
    Plot conservation scores for all CSVs in the given directory.
    """
    if conservation_dir is None:
        conservation_dir = path_config.CONSERVATION_OUTPUT_DIR
    for csv_file in Path(conservation_dir).glob("*_conservation.csv"):
        plot_conservation_scores(csv_file)


def plot_variants_on_conservation(conservation_csv, variants_csv, output_dir=None):
    """
    Overlay variant positions on conservation plot.
    """
    consv = pd.read_csv(conservation_csv)
    vars = pd.read_csv(variants_csv)

    # Try to extract integer positions from variant CSV
    def parse_pos(x):
        try:
            if "{" in str(x):
                return int(eval(x)["value"])
            return int(x)
        except Exception:
            return None

    vars["Position"] = vars["position"].apply(parse_pos)
    plt.figure(figsize=(12, 4))
    plt.plot(consv["Position"], consv[consv.columns[1]], label="Conservation")
    plt.vlines(
        vars["Position"],
        ymin=consv[consv.columns[1]].min(),
        ymax=consv[consv.columns[1]].max(),
        color="red",
        alpha=0.5,
        label="Variants",
    )
    plt.xlabel("Protein Position")
    plt.ylabel("Conservation (entropy)")
    plt.title(f"Conservation and Variant Positions: {Path(conservation_csv).stem}")
    plt.legend()
    plt.tight_layout()
    if output_dir is None:
        output_dir = Path(variants_csv).parent
    out_path = Path(output_dir) / f"{Path(conservation_csv).stem}_with_variants.png"
    plt.savefig(out_path)
    plt.close()
    print(f"Saved conservation+variant plot to {out_path}")


# Scientific plotting functions with enhanced rigor and clarity

def plot_conservation_scientific(csv_file, output_dir=None):
    """
    Create publication-quality conservation plots with statistical rigor.
    
    Args:
        csv_file (Path): Path to conservation CSV file
        output_dir (Path, optional): Output directory for plots
        
    Returns:
        Path: Path to saved scientific plot
    """
    plotter = ConservationPlotter()
    return plotter.plot_conservation_with_confidence(Path(csv_file), 
                                                   Path(output_dir) if output_dir else None)


def plot_all_conservation_scientific(conservation_dir=None):
    """
    Create scientific conservation plots for all CSVs in the given directory.
    
    Args:
        conservation_dir (Path, optional): Directory containing conservation CSV files
    """
    if conservation_dir is None:
        conservation_dir = path_config.CONSERVATION_OUTPUT_DIR
    
    plotter = ConservationPlotter()
    csv_files = list(Path(conservation_dir).glob("*_conservation.csv"))
    
    for csv_file in csv_files:
        try:
            plotter.plot_conservation_with_confidence(csv_file, Path(conservation_dir))
            print(f"Created scientific conservation plot for {csv_file.name}")
        except Exception as e:
            print(f"Error creating scientific plot for {csv_file.name}: {e}")


def visualize_trees_scientific(tree_files=None, output_dir=None):
    """
    Create publication-quality phylogenetic tree visualizations.
    
    Args:
        tree_files (list[Path], optional): List of .nwk tree file paths
        output_dir (Path, optional): Where to save PNGs
    """
    if output_dir is None:
        output_dir = path_config.TREES_OUTPUT_DIR
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if tree_files is None:
        tree_files = list(output_dir.glob("*.nwk"))

    plotter = PhylogeneticPlotter()
    
    for tree_path in tree_files:
        try:
            plotter.plot_tree_scientific(Path(tree_path), output_dir)
            print(f"Created scientific tree visualization for {tree_path.name}")
        except Exception as e:
            print(f"Error creating scientific tree plot for {tree_path.name}: {e}")


def plot_variants_scientific(conservation_csv, variants_csv, output_dir=None):
    """
    Create publication-quality variant overlay plots with statistical analysis.
    
    Args:
        conservation_csv (Path): Path to conservation scores CSV
        variants_csv (Path): Path to variants CSV  
        output_dir (Path, optional): Output directory for plots
        
    Returns:
        Path: Path to saved scientific plot
    """
    plotter = VariantPlotter()
    return plotter.plot_variants_with_statistics(Path(conservation_csv), 
                                               Path(variants_csv),
                                               Path(output_dir) if output_dir else None)
