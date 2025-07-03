"""
Conservation-Domain Mapping Visualization

Creates publication-quality plots showing conservation scores mapped to protein domains,
revealing evolutionary constraint patterns across functional regions.
"""

import json
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pathlib import Path
from typing import Dict, List, Optional
import numpy as np

from ..config.path_config import CONSERVATION_OUTPUT_DIR, DATA_OUTPUT_DIR


class ConservationDomainPlotter:
    """Visualizes conservation scores overlaid with protein domain architecture."""
    
    def __init__(self):
        self.conservation_dir = CONSERVATION_OUTPUT_DIR
        self.structures_dir = DATA_OUTPUT_DIR / "structures"
        
        # Domain type color mapping for visual consistency
        self.domain_colors = {
            'Transmembrane': '#FF6B6B',
            'Intramembrane': '#4ECDC4', 
            'Topological domain': '#45B7D1',
            'Repeat': '#96CEB4',
            'Domain': '#FFEAA7',
            'Region': '#DDA0DD',
        }
        
    def load_conservation_data(self, gene_name: str) -> pd.DataFrame:
        """Load conservation scores from CSV file."""
        conservation_file = self.conservation_dir / f"{gene_name}_conservation.csv"
        if not conservation_file.exists():
            raise FileNotFoundError(f"Conservation file not found: {conservation_file}")
        return pd.read_csv(conservation_file)
    
    def load_domain_data(self, uniprot_id: str) -> Dict:
        """Load domain annotations from JSON file."""
        domain_file = self.structures_dir / f"{uniprot_id}_domains.json"
        if not domain_file.exists():
            raise FileNotFoundError(f"Domain file not found: {domain_file}")
        with open(domain_file) as f:
            return json.load(f)
    
    def calculate_domain_conservation(self, conservation_df: pd.DataFrame, domains: List[Dict]) -> List[Dict]:
        """Calculate average conservation score for each domain."""
        domain_stats = []
        
        for domain in domains:
            start = domain['start']
            end = domain['end']
            
            # Filter conservation data for this domain region
            domain_region = conservation_df[
                (conservation_df['Position'] >= start) & 
                (conservation_df['Position'] <= end)
            ]
            
            if not domain_region.empty:
                avg_entropy = domain_region['ShannonEntropy_NoGaps'].mean()
                min_entropy = domain_region['ShannonEntropy_NoGaps'].min()
                max_entropy = domain_region['ShannonEntropy_NoGaps'].max()
                
                domain_stats.append({
                    **domain,
                    'avg_conservation': avg_entropy,
                    'min_conservation': min_entropy,
                    'max_conservation': max_entropy,
                    'conservation_range': max_entropy - min_entropy
                })
        
        return domain_stats
    
    def create_conservation_domain_plot(
        self, 
        gene_name: str, 
        uniprot_id: str,
        output_file: Optional[str] = None,
        conservation_threshold: float = 0.5
    ) -> str:
        """
        Create combined conservation-domain visualization.
        
        Args:
            gene_name: Gene name (e.g., 'SCN1A')
            uniprot_id: UniProt ID (e.g., 'P35498')
            output_file: Custom output filename
            conservation_threshold: Entropy threshold for highlighting high conservation
        
        Returns:
            Path to saved plot
        """
        # Load data
        conservation_df = self.load_conservation_data(gene_name)
        domain_data = self.load_domain_data(uniprot_id)
        domains = domain_data['domains']
        
        # Calculate domain-level conservation statistics
        domain_stats = self.calculate_domain_conservation(conservation_df, domains)
        
        # Create figure with two subplots
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 10), height_ratios=[2, 1])
        fig.suptitle(f'{gene_name} Conservation Scores Mapped to Protein Domains', 
                     fontsize=16, fontweight='bold')
        
        # Top plot: Conservation scores
        positions = conservation_df['Position']
        entropy_scores = conservation_df['ShannonEntropy_NoGaps']
        
        ax1.plot(positions, entropy_scores, linewidth=1.5, color='#2C3E50', alpha=0.8)
        ax1.fill_between(positions, entropy_scores, alpha=0.3, color='#3498DB')
        
        # Highlight highly conserved regions (low entropy)
        highly_conserved = conservation_df[conservation_df['ShannonEntropy_NoGaps'] <= conservation_threshold]
        if not highly_conserved.empty:
            ax1.scatter(highly_conserved['Position'], highly_conserved['ShannonEntropy_NoGaps'], 
                       color='#E74C3C', s=20, alpha=0.7, zorder=5, label=f'Highly conserved (≤{conservation_threshold})')
        
        ax1.set_xlabel('Amino Acid Position')
        ax1.set_ylabel('Shannon Entropy (Conservation Score)')
        ax1.set_title('Conservation Profile Across Protein Sequence')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        # Bottom plot: Domain architecture
        protein_length = domain_data['protein_length']
        
        # Group domains by type for organized display
        domain_types = {}
        for domain_stat in domain_stats:
            domain_type = domain_stat['type']
            if domain_type not in domain_types:
                domain_types[domain_type] = []
            domain_types[domain_type].append(domain_stat)
        
        # Plot domains by type on separate tracks
        y_position = 0
        track_height = 0.8
        track_spacing = 1.0
        
        for domain_type, type_domains in domain_types.items():
            color = self.domain_colors.get(domain_type, '#95A5A6')
            
            for domain in type_domains:
                start = domain['start']
                end = domain['end']
                length = domain['length']
                avg_conservation = domain['avg_conservation']
                
                # Color intensity based on conservation (darker = more conserved)
                alpha = 1.0 - min(avg_conservation / 2.0, 0.8)  # Scale alpha by conservation
                
                # Draw domain rectangle
                rect = patches.Rectangle(
                    (start, y_position), length, track_height,
                    facecolor=color, alpha=alpha, edgecolor='black', linewidth=0.5
                )
                ax2.add_patch(rect)
                
                # Add domain label if space permits
                if length > 30:  # Only label if domain is wide enough
                    label_text = domain['description'][:15] + ('...' if len(domain['description']) > 15 else '')
                    ax2.text(start + length/2, y_position + track_height/2, label_text,
                           ha='center', va='center', fontsize=8, weight='bold')
            
            # Add track label
            ax2.text(-50, y_position + track_height/2, domain_type,
                    ha='right', va='center', fontsize=10, weight='bold')
            
            y_position += track_spacing
        
        ax2.set_xlim(0, protein_length)
        ax2.set_ylim(-0.5, y_position)
        ax2.set_xlabel('Amino Acid Position')
        ax2.set_ylabel('Domain Type')
        ax2.set_title('Protein Domain Architecture (Color Intensity = Conservation Level)')
        
        # Remove y-axis ticks for cleaner appearance
        ax2.set_yticks([])
        
        # Add conservation color bar explanation
        ax2.text(protein_length * 0.02, y_position - 0.3, 
                'Domain Color Intensity: Darker = More Conserved',
                fontsize=10, style='italic', bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.7))
        
        plt.tight_layout()
        
        # Save plot
        if output_file is None:
            output_file = f"{gene_name}_{uniprot_id}_conservation_domains.png"
        
        output_path = self.structures_dir / output_file
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        return str(output_path)
    
    def generate_domain_conservation_report(self, gene_name: str, uniprot_id: str) -> Dict:
        """
        Generate statistical report of conservation by domain type.
        
        Returns:
            Dictionary with conservation statistics by domain type
        """
        conservation_df = self.load_conservation_data(gene_name)
        domain_data = self.load_domain_data(uniprot_id)
        domain_stats = self.calculate_domain_conservation(conservation_df, domain_data['domains'])
        
        # Group by domain type and calculate statistics
        type_stats = {}
        for domain in domain_stats:
            domain_type = domain['type']
            if domain_type not in type_stats:
                type_stats[domain_type] = {
                    'count': 0,
                    'avg_conservation_scores': [],
                    'total_length': 0
                }
            
            type_stats[domain_type]['count'] += 1
            type_stats[domain_type]['avg_conservation_scores'].append(domain['avg_conservation'])
            type_stats[domain_type]['total_length'] += domain['length']
        
        # Calculate summary statistics
        summary_stats = {}
        for domain_type, stats in type_stats.items():
            scores = stats['avg_conservation_scores']
            summary_stats[domain_type] = {
                'domain_count': stats['count'],
                'total_length': stats['total_length'],
                'mean_conservation': np.mean(scores),
                'std_conservation': np.std(scores),
                'min_conservation': np.min(scores),
                'max_conservation': np.max(scores)
            }
        
        return summary_stats


def create_conservation_domain_plot(gene_name: str, uniprot_id: str) -> str:
    """Convenience function to create conservation-domain plot."""
    plotter = ConservationDomainPlotter()
    return plotter.create_conservation_domain_plot(gene_name, uniprot_id)


if __name__ == "__main__":
    # Example usage
    plotter = ConservationDomainPlotter()
    
    # Generate plot for SCN1A
    output_path = plotter.create_conservation_domain_plot("SCN1A", "P35498")
    print(f"Conservation-domain plot saved to: {output_path}")
    
    # Generate conservation report
    report = plotter.generate_domain_conservation_report("SCN1A", "P35498")
    print("\nConservation by Domain Type:")
    for domain_type, stats in report.items():
        print(f"{domain_type}: {stats['mean_conservation']:.3f} ± {stats['std_conservation']:.3f}")