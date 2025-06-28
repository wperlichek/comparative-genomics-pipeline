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
import warnings

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
    
    def _add_statistical_annotations(self, ax: plt.Axes, data: pd.DataFrame, 
                                   x_col: str, y_col: str) -> None:
        """Add basic statistical annotations to plots."""
        n_points = len(data)
        mean_val = data[y_col].mean()
        std_val = data[y_col].std()
        
        # Add text box with statistics
        stats_text = f'N = {n_points}\nMean ± SD = {mean_val:.3f} ± {std_val:.3f}'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
                verticalalignment='top', bbox=dict(boxstyle='round', 
                facecolor='white', alpha=0.8, edgecolor='gray'))


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
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=self.config.figsize_conservation,
                                      height_ratios=[3, 1])
        
        # Main conservation plot
        self._plot_conservation_main(ax1, df, csv_file.stem)
        
        # Conservation distribution histogram
        self._plot_conservation_distribution(ax2, df)
        
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
        ax.set_title(f'Conservation Analysis: {title_base.replace("_", " ").title()}',
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
        """Add scientific legend with key statistics."""
        # Calculate key statistics
        gaps_mean = df['ShannonEntropy_WithGaps'].mean()
        nogaps_mean = df['ShannonEntropy_NoGaps'].mean()
        highly_conserved = (df['ShannonEntropy_NoGaps'] < 0.5).sum()
        variable_regions = (df['ShannonEntropy_NoGaps'] > 2.0).sum()
        
        # Create legend with statistics
        legend_elements = [
            mpatches.Patch(color=self.theme.primary_color, alpha=self.theme.alpha,
                          label=f'With Gaps (μ={gaps_mean:.2f})'),
            mpatches.Patch(color=self.theme.secondary_color, alpha=self.theme.alpha,
                          label=f'Without Gaps (μ={nogaps_mean:.2f})'),
            mpatches.Patch(color='white', label=f'Highly Conserved: {highly_conserved} sites'),
            mpatches.Patch(color='white', label=f'Variable Regions: {variable_regions} sites')
        ]
        
        ax.legend(handles=legend_elements, loc='upper right', 
                 frameon=True, fancybox=True, shadow=True)
    
    def _plot_conservation_distribution(self, ax: plt.Axes, df: pd.DataFrame) -> None:
        """Plot conservation score distribution."""
        entropy_nogaps = df['ShannonEntropy_NoGaps'].values
        
        # Histogram
        n_bins = min(30, len(entropy_nogaps) // 5)
        ax.hist(entropy_nogaps, bins=n_bins, color=self.theme.secondary_color,
               alpha=0.7, density=True, edgecolor='white')
        
        # Add distribution curve if we have enough data
        if len(entropy_nogaps) > 10:
            kde_x = np.linspace(entropy_nogaps.min(), entropy_nogaps.max(), 100)
            try:
                kde = stats.gaussian_kde(entropy_nogaps)
                ax.plot(kde_x, kde(kde_x), color=self.theme.error_color, 
                       linewidth=2, label='Kernel Density')
            except np.linalg.LinAlgError:
                pass  # Skip KDE if singular matrix
        
        ax.set_xlabel('Shannon Entropy (bits)', fontsize=self.theme.label_fontsize)
        ax.set_ylabel('Density', fontsize=self.theme.label_fontsize) 
        ax.set_title('Conservation Score Distribution', fontsize=self.theme.label_fontsize)
        ax.grid(True, alpha=self.theme.grid_alpha)


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
        
        # Set up the tree plot
        if self.config.tree_layout == 'circular':
            Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False,
                      branch_labels=None if not self.config.show_branch_lengths else lambda x: f'{x.branch_length:.3f}')
        else:
            Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False)
        
        # Enhance the plot
        self._enhance_tree_plot(ax, tree, tree_file.stem)
        
        output_path = output_dir / f"{tree_file.stem}_scientific.{self.config.output_format}"
        self._save_figure(fig, output_path)
        
        return output_path
    
    def _enhance_tree_plot(self, ax: plt.Axes, tree, title_base: str) -> None:
        """Enhance tree plot with scientific formatting."""
        ax.set_title(f'Phylogenetic Tree: {title_base.replace("_", " ").title()}',
                    fontsize=self.theme.title_fontsize, pad=20)
        
        # Remove axis ticks and labels for cleaner look
        ax.set_xticks([])
        ax.set_yticks([])
        
        # Add scale information if branch lengths are meaningful
        if self.config.show_branch_lengths:
            # Add scale bar (this is a simplified implementation)
            tree_info = self._get_tree_statistics(tree)
            scale_text = f"Tree Statistics:\n"
            scale_text += f"Total Length: {tree_info['total_length']:.4f}\n"
            scale_text += f"Max Distance: {tree_info['max_distance']:.4f}\n"
            scale_text += f"Taxa: {tree_info['taxa_count']}"
            
            ax.text(0.02, 0.02, scale_text, transform=ax.transAxes,
                   verticalalignment='bottom', fontsize=self.theme.tick_fontsize,
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    def _get_tree_statistics(self, tree) -> Dict[str, float]:
        """Calculate basic tree statistics."""
        terminals = tree.get_terminals()
        total_length = tree.total_branch_length()
        
        # Calculate maximum root-to-tip distance
        max_distance = 0
        for terminal in terminals:
            distance = tree.distance(tree.root, terminal)
            max_distance = max(max_distance, distance)
        
        return {
            'total_length': total_length or 0,
            'max_distance': max_distance,
            'taxa_count': len(terminals)
        }


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
        
        # Create plot
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=self.config.figsize_variants,
                                      height_ratios=[3, 1])
        
        # Main variant overlay plot
        self._plot_variant_overlay_main(ax1, consv_df, vars_df, stats_results,
                                       conservation_csv.stem)
        
        # Variant conservation distribution
        self._plot_variant_conservation_dist(ax2, consv_df, vars_df, stats_results)
        
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
        """Plot main conservation curve with variant overlay."""
        # Conservation curve
        ax.plot(consv_df['Position'], consv_df['ShannonEntropy_NoGaps'],
               color=self.theme.primary_color, linewidth=self.theme.line_width,
               alpha=self.theme.alpha, label='Conservation')
        
        # Variant positions
        variant_positions = vars_df['parsed_position'].values
        y_min, y_max = ax.get_ylim()
        
        # Identify loss-of-function variants
        lof_positions = [177, 227, 393, 939, 959, 1289]  # From grep analysis
        lof_mask = np.isin(variant_positions, lof_positions)
        regular_positions = variant_positions[~lof_mask]
        lof_variant_positions = variant_positions[lof_mask]
        
        # Plot regular variants
        if len(regular_positions) > 0:
            if stats_results.get('significant', False) and self.config.highlight_significant:
                variant_color = self.theme.accent_color
                variant_label = f'Variants (p={stats_results["p_value"]:.3e})'
            else:
                variant_color = self.theme.accent_color
                variant_label = f'Variants (n={len(regular_positions)})'
            
            ax.vlines(regular_positions, y_min, y_max,
                     colors=variant_color, alpha=self.config.variant_line_alpha,
                     linewidth=self.config.variant_line_width, label=variant_label)
        
        # Highlight loss-of-function variants in red
        if len(lof_variant_positions) > 0:
            ax.vlines(lof_variant_positions, y_min, y_max,
                     colors='red', alpha=0.8,
                     linewidth=3, label=f'Loss-of-function (n={len(lof_variant_positions)})')
        
        # Formatting
        ax.set_xlabel('Protein Position', fontsize=self.theme.label_fontsize)
        ax.set_ylabel('Conservation (Shannon Entropy)', fontsize=self.theme.label_fontsize)
        ax.set_title(f'Conservation & Variants: {title_base.replace("_", " ").title()}',
                    fontsize=self.theme.title_fontsize, pad=20)
        
        ax.grid(True, alpha=self.theme.grid_alpha)
        ax.legend(loc='upper right')
        
        # Add statistics text box
        if not stats_results.get('error'):
            self._add_statistics_textbox(ax, stats_results)
    
    def _plot_variant_conservation_dist(self, ax: plt.Axes, consv_df: pd.DataFrame,
                                      vars_df: pd.DataFrame, stats_results: Dict[str, Any]) -> None:
        """Plot conservation score distributions for variants vs background."""
        if stats_results.get('error'):
            ax.text(0.5, 0.5, 'Insufficient data for distribution analysis',
                   transform=ax.transAxes, ha='center', va='center')
            return
        
        variant_scores = stats_results['variant_conservation']
        background_scores = stats_results['background_conservation']
        
        # Create overlaid histograms
        bins = np.linspace(min(np.min(variant_scores), np.min(background_scores)),
                          max(np.max(variant_scores), np.max(background_scores)), 20)
        
        ax.hist(background_scores, bins=bins, alpha=0.7, density=True,
               color=self.theme.primary_color, label='Background', edgecolor='white')
        ax.hist(variant_scores, bins=bins, alpha=0.7, density=True,
               color=self.theme.accent_color, label='Variants', edgecolor='white')
        
        # Add mean lines
        ax.axvline(stats_results['background_mean'], color=self.theme.primary_color,
                  linestyle='--', linewidth=2, alpha=0.8)
        ax.axvline(stats_results['variant_mean'], color=self.theme.accent_color,
                  linestyle='--', linewidth=2, alpha=0.8)
        
        ax.set_xlabel('Conservation Score', fontsize=self.theme.label_fontsize)
        ax.set_ylabel('Density', fontsize=self.theme.label_fontsize)
        ax.set_title('Conservation Distribution Comparison', fontsize=self.theme.label_fontsize)
        ax.legend()
        ax.grid(True, alpha=self.theme.grid_alpha)
    
    def _add_statistics_textbox(self, ax: plt.Axes, stats_results: Dict[str, Any]) -> None:
        """Add statistical results textbox to plot."""
        if np.isnan(stats_results['p_value']):
            stats_text = "Insufficient data for statistical testing"
        else:
            p_val = stats_results['p_value']
            cohens_d = stats_results['cohens_d']
            
            stats_text = f"Statistical Analysis:\n"
            stats_text += f"Mann-Whitney U test\n"
            stats_text += f"p-value: {p_val:.3e}\n"
            stats_text += f"Effect size (d): {cohens_d:.3f}\n"
            
            if stats_results['significant']:
                stats_text += "Result: Significant difference"
            else:
                stats_text += "Result: No significant difference"
        
        ax.text(0.98, 0.02, stats_text, transform=ax.transAxes,
               verticalalignment='bottom', horizontalalignment='right',
               fontsize=self.theme.tick_fontsize,
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray'))