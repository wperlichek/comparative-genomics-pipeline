"""
Scientific visualization module for comparative genomics pipeline.

This module provides publication-quality plotting functions with statistical rigor
for genomic conservation analysis, phylogenetic visualization, and variant analysis.
"""

from .scientific_plots import ConservationPlotter, PhylogeneticPlotter, VariantPlotter
from .plot_config import PlotConfig, PlotTheme

__all__ = [
    'ConservationPlotter',
    'PhylogeneticPlotter', 
    'VariantPlotter',
    'PlotConfig',
    'PlotTheme'
]