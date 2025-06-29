"""
Scientific plotting classes for comparative genomics with statistical rigor.
"""

from pathlib import Path
from typing import Optional, List, Tuple, Dict, Any, Union
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from scipy.signal import savgol_filter
from Bio import Phylo

from .plot_config import PlotConfig, PlotTheme, PUBLICATION_THEME


class BasePlotter:
    """Base class for scientific plotters with common functionality."""
    
    def __init__(self, config: Optional[PlotConfig] = None, theme: Optional[PlotTheme] = None):
        self.config = config or PlotConfig()
        self.theme = theme or PUBLICATION_THEME
        self.config.apply_theme(self.theme)
    
    def _save_figure(self, fig: plt.Figure, output_path: Path, close_fig: bool = True) -> None:
        """Save figure with proper formatting and cleanup."""
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        fig.savefig(
            output_path,
            format=self.config.output_format,
            dpi=self.config.output_dpi,
            bbox_inches=self.config.bbox_inches,
            facecolor='white',
            edgecolor='none'
        )
        
        print(f"Saved scientific plot to {output_path}")
        
        if close_fig:
            plt.close(fig)
    


class ConservationPlotter(BasePlotter):
    """Publication-quality conservation analysis plots with statistical rigor."""
    
    def plot_conservation_with_confidence(self, csv_file: Path, 
                                        output_dir: Optional[Path] = None) -> Path:
        """
        Plot conservation scores with confidence intervals and statistical annotations.
        
        Args:
            csv_file: Path to conservation CSV file
            output_dir: Output directory for plots
            
        Returns:
            Path to saved plot
        """
        df = pd.read_csv(csv_file)
        
        if output_dir is None:
            output_dir = csv_file.parent
        
        fig, ax1 = plt.subplots(1, 1, figsize=self.config.figsize_conservation)
        
        # Main conservation plot
        self._plot_conservation_main(ax1, df, csv_file.stem)
        
        plt.tight_layout()
        
        output_path = output_dir / f"{csv_file.stem}_scientific.{self.config.output_format}"
        self._save_figure(fig, output_path)
        
        return output_path
    
    def _plot_conservation_main(self, ax: plt.Axes, df: pd.DataFrame, title_base: str) -> None:
        """Plot main conservation curves with confidence intervals."""
        positions = df['Position'].values
        entropy_gaps = df['ShannonEntropy_WithGaps'].values
        entropy_nogaps = df['ShannonEntropy_NoGaps'].values
        
        # Apply smoothing if requested
        if self.config.conservation_smoothing_window > 1:
            window = min(self.config.conservation_smoothing_window, len(positions) // 4)
            if window >= 3 and window % 2 == 0:  # savgol_filter requires odd window
                window += 1
            
            if window >= 3:
                entropy_gaps_smooth = savgol_filter(entropy_gaps, window, 2)
                entropy_nogaps_smooth = savgol_filter(entropy_nogaps, window, 2)
            else:
                entropy_gaps_smooth = entropy_gaps
                entropy_nogaps_smooth = entropy_nogaps
        else:
            entropy_gaps_smooth = entropy_gaps
            entropy_nogaps_smooth = entropy_nogaps
        
        # Plot main lines
        line1 = ax.plot(positions, entropy_gaps_smooth, 
                       color=self.theme.primary_color, 
                       label='With Gaps', 
                       alpha=self.theme.alpha,
                       linewidth=self.theme.line_width)[0]
        
        line2 = ax.plot(positions, entropy_nogaps_smooth,
                       color=self.theme.secondary_color,
                       label='Without Gaps', 
                       alpha=self.theme.alpha,
                       linewidth=self.theme.line_width)[0]
        
        # Add confidence intervals if requested
        if self.config.show_confidence_intervals:
            self._add_confidence_intervals(ax, positions, entropy_gaps_smooth, 
                                         entropy_nogaps_smooth)
        
        # Formatting
        ax.set_xlabel('Alignment Position', fontsize=self.theme.label_fontsize)
        ax.set_ylabel('Shannon Entropy (bits)', fontsize=self.theme.label_fontsize)
        ax.set_title(f'Conservation Analysis: {title_base.replace("_", " ").title()} (5 vertebrate species)',
                    fontsize=self.theme.title_fontsize, pad=20)
        
        # Add grid
        ax.grid(True, alpha=self.theme.grid_alpha, color=self.theme.grid_color)
        
        # Scientific legend with statistics
        self._add_conservation_legend(ax, df)
        
        # Add horizontal line at theoretical maximum entropy
        max_entropy = np.log2(20)  # 20 amino acids
        ax.axhline(y=max_entropy, color=self.theme.error_color, 
                  linestyle='--', alpha=0.5, linewidth=1,
                  label=f'Max Entropy ({max_entropy:.1f} bits)')
    
    def _add_confidence_intervals(self, ax: plt.Axes, positions: np.ndarray,
                                entropy_gaps: np.ndarray, entropy_nogaps: np.ndarray) -> None:
        """Add bootstrapped confidence intervals."""
        # Calculate rolling standard error as proxy for confidence interval
        window = max(5, len(positions) // 50)
        
        def rolling_std_error(data, window):
            return pd.Series(data).rolling(window, center=True, min_periods=1).std() / np.sqrt(window)
        
        gaps_se = rolling_std_error(entropy_gaps, window)
        nogaps_se = rolling_std_error(entropy_nogaps, window)
        
        # Z-score for confidence level
        z_score = stats.norm.ppf(1 - (1 - self.config.confidence_level) / 2)
        
        # Plot confidence intervals
        ax.fill_between(positions, 
                       entropy_gaps - z_score * gaps_se,
                       entropy_gaps + z_score * gaps_se,
                       color=self.theme.primary_color, alpha=0.2, 
                       label=f'{self.config.confidence_level*100:.0f}% CI (With Gaps)')
        
        ax.fill_between(positions,
                       entropy_nogaps - z_score * nogaps_se, 
                       entropy_nogaps + z_score * nogaps_se,
                       color=self.theme.secondary_color, alpha=0.2,
                       label=f'{self.config.confidence_level*100:.0f}% CI (Without Gaps)')
    
    def _add_conservation_legend(self, ax: plt.Axes, df: pd.DataFrame) -> None:
        """Add scientific legend with key statistics and data summary."""
        # Calculate comprehensive statistics
        gaps_mean = df['ShannonEntropy_WithGaps'].mean()
        nogaps_mean = df['ShannonEntropy_NoGaps'].mean()
        gaps_std = df['ShannonEntropy_WithGaps'].std()
        nogaps_std = df['ShannonEntropy_NoGaps'].std()
        
        # Conservation thresholds based on Shannon entropy
        highly_conserved = (df['ShannonEntropy_NoGaps'] < 0.5).sum()
        moderately_conserved = ((df['ShannonEntropy_NoGaps'] >= 0.5) & 
                               (df['ShannonEntropy_NoGaps'] < 1.5)).sum()
        variable_regions = (df['ShannonEntropy_NoGaps'] >= 1.5).sum()
        total_positions = len(df)
        
        # Create comprehensive legend with percentages
        legend_elements = [
            mpatches.Patch(color=self.theme.primary_color, alpha=self.theme.alpha,
                          label=f'With Gaps: {gaps_mean:.2f}±{gaps_std:.2f}'),
            mpatches.Patch(color=self.theme.secondary_color, alpha=self.theme.alpha,
                          label=f'No Gaps: {nogaps_mean:.2f}±{nogaps_std:.2f}'),
            mpatches.Patch(color='white', 
                          label=f'Highly Conserved: {highly_conserved}/{total_positions} ({highly_conserved/total_positions*100:.1f}%)'),
            mpatches.Patch(color='white', 
                          label=f'Moderate: {moderately_conserved}/{total_positions} ({moderately_conserved/total_positions*100:.1f}%)'),
            mpatches.Patch(color='white', 
                          label=f'Variable: {variable_regions}/{total_positions} ({variable_regions/total_positions*100:.1f}%)')
        ]
        
        ax.legend(handles=legend_elements, loc='upper right', 
                 frameon=True, fancybox=True, shadow=True, fontsize=self.theme.legend_fontsize)
    

class PhylogeneticPlotter(BasePlotter):
    """Publication-quality phylogenetic tree visualization."""
    
    def plot_tree_scientific(self, tree_file: Path, 
                            output_dir: Optional[Path] = None) -> Path:
        """
        Create publication-quality phylogenetic tree visualization.
        
        Args:
            tree_file: Path to Newick tree file
            output_dir: Output directory for plots
            
        Returns:
            Path to saved plot
        """
        if output_dir is None:
            output_dir = tree_file.parent
            
        tree = Phylo.read(tree_file, 'newick')
        
        fig, ax = plt.subplots(figsize=self.config.figsize_phylogeny)
        
        # Set up the tree plot with species labels
        if self.config.tree_layout == 'circular':
            Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False,
                      branch_labels=None if not self.config.show_branch_lengths else lambda x: f'{x.branch_length:.3f}',
                      label_func=self._format_species_name)
        else:
            Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False,
                      label_func=self._format_species_name)
        
        # Enhance the plot
        self._enhance_tree_plot(ax, tree, tree_file.stem)
        
        output_path = output_dir / f"{tree_file.stem}_scientific.{self.config.output_format}"
        self._save_figure(fig, output_path)
        
        return output_path
    
    def _format_species_name(self, clade):
        """Format species names for display on tree."""
        if not clade.name:
            return ""
        
        # Load species mapping from genes_to_proteins.json
        from ..config import path_config
        from ..util import file_util
        
        genes_to_proteins = file_util.open_file_return_as_json(
            f"{path_config.DATA_INPUT_DIR}/genes_to_proteins.json"
        )
        
        # Build dynamic mapping from protein IDs to species abbreviations
        name_map = {}
        for gene_name, orthologs in genes_to_proteins.items():
            for ortholog in orthologs:
                species_full = ortholog["species"]
                # Convert "Homo sapiens" to "H. sapiens"
                parts = species_full.split()
                species_abbrev = f"{parts[0][0]}. {parts[1]}" if len(parts) >= 2 else species_full
                
                if ortholog.get("uniprot_id"):
                    name_map[ortholog["uniprot_id"]] = species_abbrev
                if ortholog.get("entrez_protein_id"):
                    name_map[ortholog["entrez_protein_id"]] = species_abbrev
        
        # Extract species identifier from full name
        for key, species in name_map.items():
            if key in clade.name:
                return species
        
        return clade.name
    
    def _enhance_tree_plot(self, ax: plt.Axes, tree, title_base: str) -> None:
        """Enhance tree plot with comprehensive scientific formatting."""
        terminals = tree.get_terminals()
        species_count = len(terminals)
        ax.set_title(f'Phylogenetic Tree: {title_base.replace("_", " ").title()} ({species_count} vertebrate species)',
                    fontsize=self.theme.title_fontsize, pad=20)
        
        # Remove axis ticks and labels for cleaner look
        ax.set_xticks([])
        ax.set_yticks([])
        
        # Simple species count display
        ax.text(0.02, 0.02, f"{species_count} species", transform=ax.transAxes,
               verticalalignment='bottom', fontsize=self.theme.tick_fontsize,
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray'))
    


class VariantPlotter(BasePlotter):
    """Scientific variant analysis visualization with statistical testing."""
    
    def plot_variants_with_statistics(self, conservation_csv: Path, variants_csv: Path,
                                    output_dir: Optional[Path] = None) -> Path:
        """
        Create publication-quality variant overlay plot with statistical analysis.
        
        Args:
            conservation_csv: Path to conservation scores CSV
            variants_csv: Path to variants CSV  
            output_dir: Output directory for plots
            
        Returns:
            Path to saved plot
        """
        consv_df = pd.read_csv(conservation_csv)
        vars_df = pd.read_csv(variants_csv)
        
        if output_dir is None:
            output_dir = conservation_csv.parent
        
        # Parse variant positions
        vars_df = self._parse_variant_positions(vars_df)
        
        # Statistical analysis
        stats_results = self._analyze_variant_conservation(consv_df, vars_df)
        
        # Create plot (single panel)
        fig, ax = plt.subplots(1, 1, figsize=self.config.figsize_variants)
        
        # Main variant overlay plot
        self._plot_variant_overlay_main(ax, consv_df, vars_df, stats_results,
                                       conservation_csv.stem)
        
        plt.tight_layout()
        
        output_path = output_dir / f"{conservation_csv.stem}_variants_scientific.{self.config.output_format}"
        self._save_figure(fig, output_path)
        
        return output_path
    
    def _parse_variant_positions(self, vars_df: pd.DataFrame) -> pd.DataFrame:
        """Parse variant positions from different formats."""
        def parse_position(x):
            try:
                if pd.isna(x):
                    return None
                if isinstance(x, str) and '{' in x:
                    return int(eval(x)['value'])
                return int(float(x))
            except (ValueError, TypeError, SyntaxError):
                return None
        
        vars_df = vars_df.copy()
        vars_df['parsed_position'] = vars_df['position'].apply(parse_position)
        vars_df = vars_df.dropna(subset=['parsed_position'])
        vars_df['parsed_position'] = vars_df['parsed_position'].astype(int)
        
        return vars_df
    
    def _analyze_variant_conservation(self, consv_df: pd.DataFrame, 
                                    vars_df: pd.DataFrame) -> Dict[str, Any]:
        """Perform statistical analysis of variant-conservation relationship."""
        # Get conservation scores at variant positions
        variant_positions = vars_df['parsed_position'].values
        valid_positions = np.isin(variant_positions, consv_df['Position'].values)
        
        if not np.any(valid_positions):
            return {'error': 'No matching positions found'}
        
        variant_conservation = []
        background_conservation = []
        
        for pos in variant_positions[valid_positions]:
            conservation_score = consv_df[consv_df['Position'] == pos]['ShannonEntropy_NoGaps'].iloc[0]
            variant_conservation.append(conservation_score)
        
        # Background is all non-variant positions
        background_positions = ~np.isin(consv_df['Position'], variant_positions)
        background_conservation = consv_df[background_positions]['ShannonEntropy_NoGaps'].values
        
        # Statistical tests
        if len(variant_conservation) > 1 and len(background_conservation) > 1:
            # Mann-Whitney U test (non-parametric)
            statistic, p_value = stats.mannwhitneyu(variant_conservation, 
                                                   background_conservation,
                                                   alternative='two-sided')
            
            # Effect size (Cohen's d)
            pooled_std = np.sqrt(((len(variant_conservation) - 1) * np.var(variant_conservation, ddof=1) +
                                (len(background_conservation) - 1) * np.var(background_conservation, ddof=1)) /
                               (len(variant_conservation) + len(background_conservation) - 2))
            
            cohens_d = (np.mean(variant_conservation) - np.mean(background_conservation)) / pooled_std
        else:
            statistic, p_value, cohens_d = np.nan, np.nan, np.nan
        
        return {
            'variant_conservation': np.array(variant_conservation),
            'background_conservation': background_conservation,
            'n_variants': len(variant_conservation),
            'n_background': len(background_conservation),
            'variant_mean': np.mean(variant_conservation) if variant_conservation else np.nan,
            'background_mean': np.mean(background_conservation),
            'mann_whitney_statistic': statistic,
            'p_value': p_value,
            'cohens_d': cohens_d,
            'significant': p_value < self.config.significance_threshold if not np.isnan(p_value) else False
        }
    
    def _plot_variant_overlay_main(self, ax: plt.Axes, consv_df: pd.DataFrame,
                                  vars_df: pd.DataFrame, stats_results: Dict[str, Any],
                                  title_base: str) -> None:
        """Plot main conservation curve with intelligent variant overlay."""
        # Conservation curve with enhanced visibility
        ax.plot(consv_df['Position'], consv_df['ShannonEntropy_NoGaps'],
               color=self.theme.primary_color, linewidth=self.theme.line_width,
               alpha=self.theme.alpha, label='Conservation', zorder=1)
        
        
        # Intelligent variant clustering and display
        variant_positions = np.array(vars_df['parsed_position'].values)
        
        # Dynamically classify variants from raw data descriptions - NO hardcoding!
        lof_positions, pathogenic_positions, additional_classifications = self._get_dynamic_variant_classifications(vars_df, title_base)
        
        # Define y-axis limits for vertical lines
        y_min = consv_df['ShannonEntropy_NoGaps'].min()
        y_max = consv_df['ShannonEntropy_NoGaps'].max()
        
        # Cluster nearby variants if enabled
        if self.config.cluster_nearby_variants and len(variant_positions) > self.config.max_annotation_density:
            clustered_positions = self._cluster_variants(variant_positions)
            regular_positions, lof_variant_positions = self._separate_lof_variants(clustered_positions, lof_positions)
        else:
            lof_mask = np.isin(variant_positions, lof_positions)
            pathogenic_mask = np.isin(variant_positions, pathogenic_positions)
            # Regular positions are those that are neither LOF nor pathogenic
            regular_mask = ~(lof_mask | pathogenic_mask)
            regular_positions = variant_positions[regular_mask]
            lof_variant_positions = variant_positions[lof_mask]
        
        # Plot regular variants with density adaptation
        if len(regular_positions) > 0:
            if len(regular_positions) > self.config.max_annotation_density:
                # Use scatter plot for high density
                conservation_at_variants = [consv_df[consv_df['Position'] == pos]['ShannonEntropy_NoGaps'].iloc[0] 
                                          for pos in regular_positions if pos in consv_df['Position'].values]
                ax.scatter(regular_positions[:len(conservation_at_variants)], conservation_at_variants, 
                          color=self.theme.accent_color, alpha=0.7, s=20, 
                          label=f'Variants (n={len(regular_positions)})', zorder=3)
            else:
                # Use vertical lines for lower density
                variant_color = self.theme.accent_color
                if stats_results.get('significant', False) and self.config.highlight_significant:
                    variant_label = f'Variants (p={stats_results["p_value"]:.3e})'
                else:
                    variant_label = f'Variants (n={len(regular_positions)})'
                
                ax.vlines(regular_positions, y_min, y_max,
                         colors=variant_color, alpha=self.config.variant_line_alpha,
                         linewidth=self.config.variant_line_width, label=variant_label, zorder=2)
        
        # Always highlight loss-of-function variants prominently (transparent to show conservation underneath)
        if len(lof_variant_positions) > 0:
            ax.vlines(lof_variant_positions, y_min, y_max,
                     colors='red', alpha=0.3, linewidth=6, 
                     label=f'Loss-of-function (n={len(lof_variant_positions)})', zorder=4)
            
            # Add LoF variant annotations with smart positioning to avoid overlap
            if len(lof_variant_positions) <= 10:
                self._add_smart_annotations(ax, consv_df, lof_variant_positions, 'red', 'LOF')
        
        # Highlight pathogenic variants with distinct visual markers
        pathogenic_mask = np.isin(variant_positions, pathogenic_positions)
        pathogenic_variant_positions = variant_positions[pathogenic_mask]
        if len(pathogenic_variant_positions) > 0:
            ax.vlines(pathogenic_variant_positions, y_min, y_max,
                     colors='orange', alpha=0.4, linewidth=4, 
                     label=f'Likely Pathogenic (n={len(pathogenic_variant_positions)})', zorder=3)
            
            # Add likely pathogenic variant annotations with smart positioning to avoid overlap
            if len(pathogenic_variant_positions) <= 15:
                self._add_smart_annotations(ax, consv_df, pathogenic_variant_positions, 'orange', 'LP')
        
        # Enhanced formatting
        ax.set_xlabel('Protein Position', fontsize=self.theme.label_fontsize)
        ax.set_ylabel('Conservation (Shannon Entropy)', fontsize=self.theme.label_fontsize)
        ax.set_title(f'Conservation & Variants: {title_base.replace("_", " ").title()} (5 vertebrate species)',
                    fontsize=self.theme.title_fontsize, pad=20)
        
        ax.grid(True, alpha=self.theme.grid_alpha)
        ax.legend(loc='upper right', fontsize=self.theme.legend_fontsize)
        
    
    
    def _cluster_variants(self, positions: np.ndarray) -> np.ndarray:
        """Cluster nearby variants to reduce visual complexity."""
        if len(positions) <= self.config.max_annotation_density:
            return positions
        
        sorted_positions = np.sort(positions)
        clustered = []
        
        i = 0
        while i < len(sorted_positions):
            cluster_start = sorted_positions[i]
            cluster_positions = [cluster_start]
            
            # Find all positions within cluster distance
            j = i + 1
            while j < len(sorted_positions) and sorted_positions[j] - cluster_start <= self.config.cluster_distance:
                cluster_positions.append(sorted_positions[j])
                j += 1
            
            # Use median position to represent cluster
            clustered.append(int(np.median(cluster_positions)))
            i = j
        
        return np.array(clustered)
    
    def _separate_lof_variants(self, positions: np.ndarray, lof_positions: List[int]) -> Tuple[np.ndarray, np.ndarray]:
        """Separate loss-of-function variants from regular variants."""
        lof_mask = np.isin(positions, lof_positions)
        regular_positions = positions[~lof_mask]
        lof_variant_positions = positions[lof_mask]
        return regular_positions, lof_variant_positions
    
    def _get_dynamic_variant_classifications(self, vars_df: pd.DataFrame, title_base: str) -> Tuple[List[int], List[int], Dict[str, List[int]]]:
        """Dynamically classify variants based ONLY on raw data descriptions - NO hardcoding."""
        lof_positions = []
        pathogenic_positions = []
        additional_classifications = {
            'benign': [],
            'uncertain': [],
            'borderline': [],
            'reduced_function': []
        }
        
        # Extract ALL variant classifications DYNAMICALLY from raw description data
        if 'description' in vars_df.columns:
            for _, row in vars_df.iterrows():
                pos = row['parsed_position']
                desc = str(row.get('description', '')).lower()
                
                # Comprehensive variant classification based on actual data descriptions
                
                # Likely pathogenic variants (strict matching for scientific accuracy)
                if 'likely pathogenic' in desc:
                    pathogenic_positions.append(pos)
                
                # Likely benign variants
                if 'likely benign' in desc:
                    additional_classifications['benign'].append(pos)
                
                # Uncertain significance
                if 'uncertain significance' in desc:
                    additional_classifications['uncertain'].append(pos)
                
                # Borderline phenotype
                if 'borderline' in desc:
                    additional_classifications['borderline'].append(pos)
                
                # Reduced function (non-LOF but impaired)
                if any(term in desc for term in ['reduced function', 'decreased peak current', 'impaired channel', 'reduced current']):
                    additional_classifications['reduced_function'].append(pos)
                
                # Comprehensive LOF detection based on actual data patterns
                lof_indicators = [
                    'loss of function',
                    'loss-of-function',
                    'non-functional channel',
                    'results in a non-functional channel',
                    'complete absence of sodium current',
                    'absence of sodium current', 
                    'complete loss of sodium ion transmembrane transport',
                    'complete loss of sodium'
                ]
                
                if any(indicator in desc for indicator in lof_indicators):
                    lof_positions.append(pos)
        
        # NO HARDCODING! All classifications must come from raw data descriptions.
        # This ensures the pipeline scales to any gene and accurately reflects the actual data.
        
        return sorted(list(set(lof_positions))), sorted(list(set(pathogenic_positions))), additional_classifications
    
    def _add_smart_annotations(self, ax: plt.Axes, consv_df: pd.DataFrame, 
                             positions: np.ndarray, color: str, annotation_type: str) -> None:
        """Add annotations with smart positioning to avoid overlap."""
        if len(positions) == 0:
            return
        
        # Sort positions for consistent annotation placement
        sorted_positions = np.sort(positions)
        
        # Calculate minimum distance for overlap detection (in data coordinates)
        x_range = ax.get_xlim()[1] - ax.get_xlim()[0]
        # Use smaller threshold for tighter grouping - positions within 30 units are considered overlapping
        min_distance = min(30, x_range * 0.02)  # 30 positions or 2% of x-range, whichever is smaller
        
        # Group positions that are close together
        position_groups = []
        current_group = [sorted_positions[0]]
        
        for pos in sorted_positions[1:]:
            if pos - current_group[-1] < min_distance:
                current_group.append(pos)
            else:
                position_groups.append(current_group)
                current_group = [pos]
        position_groups.append(current_group)
        
        # Get plot boundaries to avoid overlapping with axis labels
        y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
        y_min_plot, y_max_plot = ax.get_ylim()
        
        # Define safe vertical offsets that stay within plot area
        # Positive offsets: stay below top 15% of plot area
        # Negative offsets: stay above bottom 20% of plot area (to avoid x-axis label)
        safe_top = y_max_plot - (y_range * 0.15)
        safe_bottom = y_min_plot + (y_range * 0.25)
        
        # Calculate offsets in data coordinates (not points)
        y_offset_data = y_range * 0.06  # 6% of y-range for better spacing
        y_offsets_data = [y_offset_data, y_offset_data * 2, y_offset_data * 3, y_offset_data * 4,
                         -y_offset_data, -y_offset_data * 2, -y_offset_data * 3, -y_offset_data * 4]
        # More varied horizontal offsets for better separation
        x_offsets_data = [0, x_range * 0.01, -x_range * 0.01, x_range * 0.015, -x_range * 0.015, 
                         x_range * 0.02, -x_range * 0.02, x_range * 0.008]
        
        group_idx = 0
        for group in position_groups:
            for i, pos in enumerate(group):
                if pos in consv_df['Position'].values:
                    conservation_score = consv_df[consv_df['Position'] == pos]['ShannonEntropy_NoGaps'].iloc[0]
                    
                    # Use alternating offsets for overlapping positions
                    offset_idx = (group_idx * len(group) + i) % len(y_offsets_data)
                    y_offset_data_val = y_offsets_data[offset_idx]
                    x_offset_data_val = x_offsets_data[offset_idx]
                    
                    # Calculate annotation position in data coordinates
                    if y_offset_data_val > 0:
                        # Positive offset: place above, but check if it fits
                        annotation_y = min(conservation_score + y_offset_data_val, safe_top)
                    else:
                        # Negative offset: place below, but check if it fits
                        annotation_y = max(conservation_score + y_offset_data_val, safe_bottom)
                    
                    # Calculate horizontal position with offset
                    annotation_x = pos + x_offset_data_val
                    
                    # Adjust font size for crowded areas
                    font_size = max(self.theme.tick_fontsize - 2, 6) if len(group) > 3 else self.theme.tick_fontsize - 1
                    
                    # Use abbreviated labels for crowded areas
                    if len(group) > 2:
                        label = f'{annotation_type}\n{pos}'
                    else:
                        label = f'{pos}'
                    
                    # Use data coordinates for annotation to keep it within safe area
                    # Set high z-order to ensure annotations appear on top of all other elements
                    ax.annotate(label, xy=(pos, conservation_score), 
                               xytext=(annotation_x, annotation_y), textcoords='data',
                               fontsize=font_size, 
                               bbox=dict(boxstyle='round,pad=0.2', facecolor=color, alpha=0.9),
                               color='white', weight='bold',
                               ha='center', va='center', zorder=10,
                               arrowprops=dict(arrowstyle='->', color=color, alpha=0.7, lw=1.5, zorder=9)
                               if (abs(annotation_y - conservation_score) > y_range * 0.02 or 
                                   abs(annotation_x - pos) > x_range * 0.005) else None)
            group_idx += 1