from pathlib import Path
from Bio import Phylo, AlignIO
from Bio.Align import AlignInfo
import matplotlib.pyplot as plt
import numpy as np
import csv
import pandas as pd
import logging
from ..config import path_config
from ..visualization import ConservationPlotter, PhylogeneticPlotter, VariantPlotter

# Configure logger for this module
logger = logging.getLogger(__name__)


def visualize_and_save_trees(tree_files=None, output_dir=None):
    """
    Visualize each Newick tree in the list and save as PNG to the trees output folder.
    Args:
        tree_files (list[Path] or None): List of .nwk tree file paths. If None, will glob all in TREES_OUTPUT_DIR.
        output_dir (Path or None): Where to save PNGs. If None, uses TREES_OUTPUT_DIR.
    
    Returns:
        bool: True if all visualizations succeeded, False if any failed
    """
    success_count = 0
    total_count = 0
    
    try:
        if output_dir is None:
            output_dir = path_config.TREES_OUTPUT_DIR
        
        # Ensure output directory exists with proper error handling
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
        except (OSError, PermissionError) as e:
            logger.error(f"Failed to create output directory {output_dir}: {e}")
            return False

        if tree_files is None:
            try:
                tree_files = list(output_dir.glob("*.nwk"))
            except (OSError, PermissionError) as e:
                logger.error(f"Failed to list tree files in {output_dir}: {e}")
                return False
        
        if not tree_files:
            logger.warning("No tree files found to visualize")
            return True

        for tree_path in tree_files:
            total_count += 1
            try:
                # Check if file exists and is readable
                if not tree_path.exists():
                    logger.error(f"Tree file does not exist: {tree_path}")
                    continue
                
                if not tree_path.is_file():
                    logger.error(f"Path is not a file: {tree_path}")
                    continue
                
                # Read tree with error handling for malformed files
                try:
                    tree = Phylo.read(tree_path, "newick")
                except (ValueError, FileNotFoundError, IOError) as e:
                    logger.error(f"Failed to read tree file {tree_path}: {e}")
                    continue
                except Exception as e:
                    logger.error(f"Unexpected error reading tree file {tree_path}: {e}")
                    continue
                
                # Create matplotlib figure with error handling
                try:
                    fig = plt.figure(figsize=(10, 5))
                    Phylo.draw(tree, do_show=False)
                    
                    png_path = output_dir / f"{tree_path.stem}.png"
                    
                    # Save figure with error handling
                    plt.savefig(png_path, bbox_inches="tight")
                    plt.close(fig)
                    
                    print(f"Saved visualization for {tree_path.name} to {png_path}")
                    success_count += 1
                    
                except (IOError, OSError, PermissionError) as e:
                    logger.error(f"Failed to save tree visualization for {tree_path.name}: {e}")
                    plt.close(fig)  # Ensure figure is closed even on error
                    continue
                except Exception as e:
                    logger.error(f"Unexpected error creating tree visualization for {tree_path.name}: {e}")
                    try:
                        plt.close(fig)
                    except:
                        pass
                    continue
                    
            except Exception as e:
                logger.error(f"Unexpected error processing tree file {tree_path}: {e}")
                continue
        
        if success_count == 0 and total_count > 0:
            logger.error(f"Failed to visualize any of {total_count} tree files")
            return False
        elif success_count < total_count:
            logger.warning(f"Successfully visualized {success_count}/{total_count} tree files")
        
        return success_count > 0 or total_count == 0
        
    except Exception as e:
        logger.error(f"Unexpected error in visualize_and_save_trees: {e}")
        return False


def compute_conservation_scores(msa_file, output_file=None):
    """
    Compute conservation (Shannon entropy) for each column in the MSA and save as CSV.
    Outputs both entropy with and without gaps.
    
    Returns:
        Path or None: Path to output file if successful, None if failed
    """
    try:
        # Validate input file
        msa_file = Path(msa_file)
        if not msa_file.exists():
            logger.error(f"MSA file does not exist: {msa_file}")
            return None
        
        if not msa_file.is_file():
            logger.error(f"MSA path is not a file: {msa_file}")
            return None
        
        # Read alignment with error handling for malformed FASTA files
        try:
            alignment = AlignIO.read(str(msa_file), "fasta")
        except (ValueError, FileNotFoundError, IOError) as e:
            logger.error(f"Failed to read MSA file {msa_file}: {e}")
            return None
        except Exception as e:
            logger.error(f"Unexpected error reading MSA file {msa_file}: {e}")
            return None
        
        # Validate alignment data
        if len(alignment) == 0:
            logger.error(f"MSA file {msa_file} contains no sequences")
            return None
        
        try:
            aln_len = alignment.get_alignment_length()
        except Exception as e:
            logger.error(f"Failed to get alignment length from {msa_file}: {e}")
            return None
        
        if aln_len == 0:
            logger.error(f"MSA file {msa_file} has zero alignment length")
            return None
        
        print(
            f"DEBUG: Alignment has {len(alignment)} sequences, length {aln_len}"
        )
        
        # Print first 5 columns for inspection (with error handling)
        try:
            inspection_cols = min(5, aln_len)
            for i in range(inspection_cols):
                try:
                    column_full = str(alignment[:, i])
                    column_nogap = column_full.replace("-", "")
                    print(f"Col {i+1}: full='{column_full}', nogap='{column_nogap}'")
                    print(f"  unique (full): {set(column_full)}")
                    print(f"  unique (nogap): {set(column_nogap)}")
                except Exception as e:
                    logger.warning(f"Error inspecting column {i+1}: {e}")
                    continue
        except Exception as e:
            logger.warning(f"Error during column inspection: {e}")
        
        scores = []
        failed_positions = 0
        
        for i in range(aln_len):
            try:
                column_full = str(alignment[:, i])  # includes gaps
                column_nogap = column_full.replace("-", "")  # excludes gaps

                # With gaps - handle empty columns and mathematical errors
                total_full = len(column_full)
                if total_full == 0:
                    entropy_full = 0.0
                else:
                    try:
                        unique_chars = set(column_full)
                        if len(unique_chars) == 0:
                            entropy_full = 0.0
                        else:
                            freqs_full = [column_full.count(aa) / total_full for aa in unique_chars]
                            # Handle log2(0) cases
                            entropy_full = -sum(p * np.log2(p) for p in freqs_full if p > 0)
                            entropy_full = abs(entropy_full)
                    except (ValueError, ZeroDivisionError, np.core._exceptions._ArrayMemoryError) as e:
                        logger.warning(f"Error computing entropy with gaps for position {i+1}: {e}")
                        entropy_full = 0.0

                # Without gaps - handle empty columns and mathematical errors
                total_nogap = len(column_nogap)
                if total_nogap == 0:
                    entropy_nogap = 0.0
                    consensus = "-"
                else:
                    try:
                        unique_chars_nogap = set(column_nogap)
                        if len(unique_chars_nogap) == 0:
                            entropy_nogap = 0.0
                            consensus = "-"
                        else:
                            freqs_nogap = [
                                column_nogap.count(aa) / total_nogap for aa in unique_chars_nogap
                            ]
                            # Handle log2(0) cases
                            entropy_nogap = -sum(p * np.log2(p) for p in freqs_nogap if p > 0)
                            entropy_nogap = abs(entropy_nogap)
                            consensus = max(unique_chars_nogap, key=column_nogap.count)
                    except (ValueError, ZeroDivisionError, np.core._exceptions._ArrayMemoryError) as e:
                        logger.warning(f"Error computing entropy without gaps for position {i+1}: {e}")
                        entropy_nogap = 0.0
                        consensus = "-"

                scores.append((i + 1, entropy_full, entropy_nogap, consensus))
                
            except Exception as e:
                logger.warning(f"Error processing alignment position {i+1}: {e}")
                failed_positions += 1
                # Add default values for failed position
                scores.append((i + 1, 0.0, 0.0, "-"))
                continue
        
        if failed_positions > 0:
            logger.warning(f"Failed to process {failed_positions}/{aln_len} alignment positions")
        
        if len(scores) == 0:
            logger.error(f"No conservation scores computed for {msa_file}")
            return None

        # Prepare output file path
        if output_file is None:
            try:
                output_dir = path_config.CONSERVATION_OUTPUT_DIR
                output_dir.mkdir(parents=True, exist_ok=True)
                output_file = output_dir / f"{msa_file.stem}_conservation.csv"
            except (OSError, PermissionError) as e:
                logger.error(f"Failed to create conservation output directory: {e}")
                return None
        
        # Write CSV file with error handling
        try:
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
        except (IOError, OSError, PermissionError) as e:
            logger.error(f"Failed to write conservation scores to {output_file}: {e}")
            return None
        except Exception as e:
            logger.error(f"Unexpected error writing conservation scores to {output_file}: {e}")
            return None
        
        print(f"Saved conservation scores to {output_file}")
        return output_file
        
    except Exception as e:
        logger.error(f"Unexpected error in compute_conservation_scores: {e}")
        return None


def compute_conservation_for_all_msas(msa_dir=None):
    """
    Compute conservation scores for all MSA files in the given directory.
    
    Returns:
        bool: True if all computations succeeded, False if any failed
    """
    success_count = 0
    total_count = 0
    
    try:
        if msa_dir is None:
            msa_dir = path_config.MSA_OUTPUT_DIR
        
        msa_dir = Path(msa_dir)
        if not msa_dir.exists():
            logger.error(f"MSA directory does not exist: {msa_dir}")
            return False
        
        if not msa_dir.is_dir():
            logger.error(f"MSA path is not a directory: {msa_dir}")
            return False
        
        try:
            msa_files = list(msa_dir.glob("*.fasta"))
        except (OSError, PermissionError) as e:
            logger.error(f"Failed to list MSA files in {msa_dir}: {e}")
            return False
        
        if not msa_files:
            logger.warning(f"No FASTA files found in {msa_dir}")
            return True
        
        for msa_file in msa_files:
            total_count += 1
            try:
                result = compute_conservation_scores(msa_file)
                if result is not None:
                    success_count += 1
                else:
                    logger.error(f"Failed to compute conservation scores for {msa_file.name}")
            except Exception as e:
                logger.error(f"Unexpected error processing {msa_file.name}: {e}")
        
        if success_count == 0 and total_count > 0:
            logger.error(f"Failed to compute conservation scores for any of {total_count} MSA files")
            return False
        elif success_count < total_count:
            logger.warning(f"Successfully computed conservation scores for {success_count}/{total_count} MSA files")
        
        return success_count > 0 or total_count == 0
        
    except Exception as e:
        logger.error(f"Unexpected error in compute_conservation_for_all_msas: {e}")
        return False


def plot_conservation_scores(csv_file, output_dir=None):
    """
    Plot Shannon entropy (with and without gaps) from a conservation CSV file.
    Saves the plot as a PNG in the output_dir (defaults to CONSERVATION_OUTPUT_DIR).
    
    Returns:
        Path or None: Path to saved plot if successful, None if failed
    """
    fig = None
    try:
        # Validate input file
        csv_file = Path(csv_file)
        if not csv_file.exists():
            logger.error(f"Conservation CSV file does not exist: {csv_file}")
            return None
        
        if not csv_file.is_file():
            logger.error(f"Conservation CSV path is not a file: {csv_file}")
            return None
        
        # Prepare output directory
        if output_dir is None:
            output_dir = path_config.CONSERVATION_OUTPUT_DIR
        output_dir = Path(output_dir)
        
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
        except (OSError, PermissionError) as e:
            logger.error(f"Failed to create conservation plot output directory {output_dir}: {e}")
            return None
        
        # Read CSV file with error handling
        try:
            df = pd.read_csv(csv_file)
        except (FileNotFoundError, pd.errors.EmptyDataError, pd.errors.ParserError) as e:
            logger.error(f"Failed to read conservation CSV file {csv_file}: {e}")
            return None
        except Exception as e:
            logger.error(f"Unexpected error reading conservation CSV file {csv_file}: {e}")
            return None
        
        # Validate CSV data
        required_columns = ["Position", "ShannonEntropy_WithGaps", "ShannonEntropy_NoGaps"]
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            logger.error(f"Conservation CSV file {csv_file} missing required columns: {missing_columns}")
            return None
        
        if len(df) == 0:
            logger.error(f"Conservation CSV file {csv_file} contains no data")
            return None
        
        # Create plot with error handling
        try:
            fig = plt.figure(figsize=(12, 5))
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
            
            # Save plot with error handling
            plt.savefig(png_path)
            plt.close(fig)
            fig = None  # Mark as closed
            
            print(f"Saved conservation plot to {png_path}")
            return png_path
            
        except (IOError, OSError, PermissionError) as e:
            logger.error(f"Failed to save conservation plot for {csv_file.name}: {e}")
            if fig is not None:
                plt.close(fig)
            return None
        except Exception as e:
            logger.error(f"Unexpected error creating conservation plot for {csv_file.name}: {e}")
            if fig is not None:
                plt.close(fig)
            return None
        
    except Exception as e:
        logger.error(f"Unexpected error in plot_conservation_scores: {e}")
        if fig is not None:
            try:
                plt.close(fig)
            except:
                pass
        return None


def plot_all_conservation_scores(conservation_dir=None):
    """
    Plot conservation scores for all CSVs in the given directory.
    
    Returns:
        bool: True if all plots succeeded, False if any failed
    """
    success_count = 0
    total_count = 0
    
    try:
        if conservation_dir is None:
            conservation_dir = path_config.CONSERVATION_OUTPUT_DIR
        
        conservation_dir = Path(conservation_dir)
        if not conservation_dir.exists():
            logger.error(f"Conservation directory does not exist: {conservation_dir}")
            return False
        
        if not conservation_dir.is_dir():
            logger.error(f"Conservation path is not a directory: {conservation_dir}")
            return False
        
        try:
            csv_files = list(conservation_dir.glob("*_conservation.csv"))
        except (OSError, PermissionError) as e:
            logger.error(f"Failed to list conservation CSV files in {conservation_dir}: {e}")
            return False
        
        if not csv_files:
            logger.warning(f"No conservation CSV files found in {conservation_dir}")
            return True
        
        for csv_file in csv_files:
            total_count += 1
            try:
                result = plot_conservation_scores(csv_file)
                if result is not None:
                    success_count += 1
                else:
                    logger.error(f"Failed to plot conservation scores for {csv_file.name}")
            except Exception as e:
                logger.error(f"Unexpected error plotting conservation scores for {csv_file.name}: {e}")
        
        if success_count == 0 and total_count > 0:
            logger.error(f"Failed to plot conservation scores for any of {total_count} CSV files")
            return False
        elif success_count < total_count:
            logger.warning(f"Successfully plotted conservation scores for {success_count}/{total_count} CSV files")
        
        return success_count > 0 or total_count == 0
        
    except Exception as e:
        logger.error(f"Unexpected error in plot_all_conservation_scores: {e}")
        return False


def plot_variants_on_conservation(conservation_csv, variants_csv, output_dir=None):
    """
    Overlay variant positions on conservation plot.
    
    Returns:
        Path or None: Path to saved plot if successful, None if failed
    """
    fig = None
    try:
        # Validate input files
        conservation_csv = Path(conservation_csv)
        variants_csv = Path(variants_csv)
        
        if not conservation_csv.exists():
            logger.error(f"Conservation CSV file does not exist: {conservation_csv}")
            return None
        
        if not conservation_csv.is_file():
            logger.error(f"Conservation CSV path is not a file: {conservation_csv}")
            return None
        
        if not variants_csv.exists():
            logger.error(f"Variants CSV file does not exist: {variants_csv}")
            return None
        
        if not variants_csv.is_file():
            logger.error(f"Variants CSV path is not a file: {variants_csv}")
            return None
        
        # Read CSV files with error handling
        try:
            consv = pd.read_csv(conservation_csv)
        except (FileNotFoundError, pd.errors.EmptyDataError, pd.errors.ParserError) as e:
            logger.error(f"Failed to read conservation CSV file {conservation_csv}: {e}")
            return None
        except Exception as e:
            logger.error(f"Unexpected error reading conservation CSV file {conservation_csv}: {e}")
            return None
        
        try:
            vars = pd.read_csv(variants_csv)
        except (FileNotFoundError, pd.errors.EmptyDataError, pd.errors.ParserError) as e:
            logger.error(f"Failed to read variants CSV file {variants_csv}: {e}")
            return None
        except Exception as e:
            logger.error(f"Unexpected error reading variants CSV file {variants_csv}: {e}")
            return None
        
        # Validate CSV data
        if len(consv) == 0:
            logger.error(f"Conservation CSV file {conservation_csv} contains no data")
            return None
        
        if len(vars) == 0:
            logger.error(f"Variants CSV file {variants_csv} contains no data")
            return None
        
        if "Position" not in consv.columns:
            logger.error(f"Conservation CSV file {conservation_csv} missing 'Position' column")
            return None
        
        if len(consv.columns) < 2:
            logger.error(f"Conservation CSV file {conservation_csv} has insufficient columns")
            return None
        
        if "position" not in vars.columns:
            logger.error(f"Variants CSV file {variants_csv} missing 'position' column")
            return None

        # Try to extract integer positions from variant CSV
        def parse_pos(x):
            try:
                if "{" in str(x):
                    return int(eval(x)["value"])
                return int(x)
            except Exception:
                return None

        try:
            vars["Position"] = vars["position"].apply(parse_pos)
            # Filter out None values
            vars = vars.dropna(subset=["Position"])
            
            if len(vars) == 0:
                logger.warning(f"No valid variant positions found in {variants_csv}")
                # Continue with empty variant data
                
        except Exception as e:
            logger.error(f"Error parsing variant positions from {variants_csv}: {e}")
            return None
        
        # Prepare output directory
        if output_dir is None:
            output_dir = variants_csv.parent
        output_dir = Path(output_dir)
        
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
        except (OSError, PermissionError) as e:
            logger.error(f"Failed to create variant plot output directory {output_dir}: {e}")
            return None
        
        # Create plot with error handling
        try:
            fig = plt.figure(figsize=(12, 4))
            
            # Plot conservation data
            conservation_column = consv.columns[1]  # Use second column for conservation values
            plt.plot(consv["Position"], consv[conservation_column], label="Conservation")
            
            # Plot variant positions if any exist
            if len(vars) > 0:
                try:
                    conservation_min = consv[conservation_column].min()
                    conservation_max = consv[conservation_column].max()
                    
                    plt.vlines(
                        vars["Position"],
                        ymin=conservation_min,
                        ymax=conservation_max,
                        color="red",
                        alpha=0.5,
                        label="Variants",
                    )
                except Exception as e:
                    logger.warning(f"Error plotting variant lines: {e}")
                    # Continue without variant lines
            
            plt.xlabel("Protein Position")
            plt.ylabel("Conservation (entropy)")
            plt.title(f"Conservation and Variant Positions: {conservation_csv.stem}")
            plt.legend()
            plt.tight_layout()
            
            out_path = output_dir / f"{conservation_csv.stem}_with_variants.png"
            
            # Save plot with error handling
            plt.savefig(out_path)
            plt.close(fig)
            fig = None  # Mark as closed
            
            print(f"Saved conservation+variant plot to {out_path}")
            return out_path
            
        except (IOError, OSError, PermissionError) as e:
            logger.error(f"Failed to save conservation+variant plot: {e}")
            if fig is not None:
                plt.close(fig)
            return None
        except Exception as e:
            logger.error(f"Unexpected error creating conservation+variant plot: {e}")
            if fig is not None:
                plt.close(fig)
            return None
        
    except Exception as e:
        logger.error(f"Unexpected error in plot_variants_on_conservation: {e}")
        if fig is not None:
            try:
                plt.close(fig)
            except:
                pass
        return None


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
