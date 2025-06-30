"""
Configuration classes for scientific plotting with publication-quality defaults.
"""

from dataclasses import dataclass
from typing import Dict, Tuple, Optional, Any
import matplotlib.pyplot as plt
import matplotlib as mpl


@dataclass
class PlotTheme:
    """Scientific publication theme configuration."""
    
    # Color schemes
    primary_color: str = '#2E86AB'
    secondary_color: str = '#A23B72'  
    accent_color: str = '#F18F01'
    error_color: str = '#C73E1D'
    grid_color: str = '#E5E5E5'
    text_color: str = '#2D3748'
    
    # Typography
    font_family: str = 'DejaVu Sans'
    title_fontsize: int = 14
    label_fontsize: int = 12
    tick_fontsize: int = 10
    legend_fontsize: int = 10
    
    # Layout
    figure_dpi: int = 300
    line_width: float = 1.5
    marker_size: float = 4.0
    alpha: float = 0.8
    
    # Grid and spines
    show_grid: bool = True
    grid_alpha: float = 0.3
    spine_width: float = 0.8


@dataclass 
class PlotConfig:
    """Configuration for specific plot types."""
    
    # Figure dimensions (inches) - optimized for information density
    figsize_conservation: Tuple[int, int] = (14, 8)
    figsize_phylogeny: Tuple[int, int] = (12, 10) 
    figsize_variants: Tuple[int, int] = (16, 8)
    
    # Conservation plot specific
    conservation_smoothing_window: int = 5
    show_confidence_intervals: bool = True
    confidence_level: float = 0.95
    show_data_summary: bool = True
    max_annotation_density: int = 50  # Maximum annotations before summarizing
    
    # Phylogeny specific
    tree_layout: str = 'rectangular'  # 'rectangular', 'circular', 'radial'
    show_branch_lengths: bool = True
    show_internal_labels: bool = False
    
    # Variant overlay specific
    variant_line_alpha: float = 0.6
    variant_line_width: float = 1.0
    highlight_significant: bool = True
    significance_threshold: float = 0.05
    show_variant_summary: bool = True
    cluster_nearby_variants: bool = True
    cluster_distance: int = 10  # Group variants within this distance
    
    # Variant classification colors
    variant_colors: Dict[str, str] = None
    
    # Output settings
    output_format: str = 'png'
    output_dpi: int = 300
    bbox_inches: str = 'tight'
    
    def __post_init__(self):
        if self.variant_colors is None:
            self.variant_colors = {
                'pathogenic': '#C73E1D',
                'likely_pathogenic': '#F18F01', 
                'uncertain': '#808080',
                'likely_benign': '#A8D5BA',
                'benign': '#2E86AB',
                'other': '#9B9B9B'
            }
    
    def apply_theme(self, theme: PlotTheme) -> None:
        """Apply scientific theme to matplotlib."""
        plt.style.use('default')  # Reset to default first
        
        # Set global parameters
        mpl.rcParams.update({
            'font.family': theme.font_family,
            'font.size': theme.label_fontsize,
            'axes.titlesize': theme.title_fontsize,
            'axes.labelsize': theme.label_fontsize,
            'xtick.labelsize': theme.tick_fontsize,
            'ytick.labelsize': theme.tick_fontsize,
            'legend.fontsize': theme.legend_fontsize,
            'figure.dpi': theme.figure_dpi,
            'savefig.dpi': self.output_dpi,
            'savefig.bbox': self.bbox_inches,
            'lines.linewidth': theme.line_width,
            'lines.markersize': theme.marker_size,
            'axes.linewidth': theme.spine_width,
            'axes.spines.top': False,
            'axes.spines.right': False,
            'axes.grid': theme.show_grid,
            'grid.alpha': theme.grid_alpha,
            'grid.color': theme.grid_color,
            'text.color': theme.text_color,
            'axes.labelcolor': theme.text_color,
            'xtick.color': theme.text_color,
            'ytick.color': theme.text_color,
        })


# ClinVar plotting constants
CLINICAL_SIGNIFICANCE_MAPPING = {
    'pathogenic': 'Pathogenic',
    'likely_pathogenic': 'Likely Pathogenic',
    'benign': 'Benign',
    'likely_benign': 'Likely Benign',
    'uncertain': 'Uncertain',
    'vus': 'Uncertain',
    'other': 'Other'
}

PLOT_POSITIONING = {
    'example_text_x': 0.02,
    'example_text_y': 0.98,
    'bar_label_offset': 0.5,
    'example_box_padding': 0.3
}


# Predefined themes
PUBLICATION_THEME = PlotTheme()

