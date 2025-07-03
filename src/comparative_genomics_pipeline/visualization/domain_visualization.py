import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pathlib import Path
from typing import Dict, List, Any, Optional
import logging

logger = logging.getLogger(__name__)

def visualize_protein_domains(domain_data: Dict[str, Any], out_dir: Optional[Path] = None) -> Optional[Path]:
    """
    Create a clean visualization of protein domains showing their positions along the protein sequence.
    
    Args:
        domain_data: Dictionary containing protein length and domain information
        out_dir: Output directory for the plot (defaults to structures directory)
    
    Returns:
        Path to the saved plot file, or None on failure
    """
    if out_dir is None:
        from ..config import path_config
        out_dir = path_config.DATA_OUTPUT_DIR / "structures"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    accession = domain_data.get('accession')
    protein_length = domain_data.get('protein_length', 0)
    domains = domain_data.get('domains', [])
    
    if not accession or not protein_length:
        logger.error("Invalid domain data: missing accession or protein_length")
        return None
    
    # Set up the plot with more vertical space
    fig, ax = plt.subplots(figsize=(16, 10))
    
    # Color map for different domain types
    type_colors = {
        'Domain': '#FF6B6B',
        'Transmembrane': '#4ECDC4', 
        'Repeat': '#45B7D1',
        'Topological domain': '#FFA07A',
        'Region': '#98D8C8',
        'Intramembrane': '#F7DC6F'
    }
    
    # Group domains by type for better organization
    domain_groups = {}
    for domain in domains:
        domain_type = domain['type']
        if domain_type not in domain_groups:
            domain_groups[domain_type] = []
        domain_groups[domain_type].append(domain)
    
    # Define track positions for each domain type (higher number = higher on plot)
    track_positions = {
        'Repeat': 8,
        'Domain': 7,
        'Transmembrane': 6,
        'Intramembrane': 5,
        'Topological domain': 4,
        'Region': 3
    }
    
    # Draw protein backbone at the bottom
    backbone_y = 1
    backbone_height = 0.3
    ax.add_patch(patches.Rectangle((0, backbone_y), protein_length, backbone_height, 
                                   facecolor='lightgray', edgecolor='black', linewidth=2))
    ax.text(protein_length/2, backbone_y + backbone_height/2, f'SCN1A Protein ({protein_length} aa)', 
            ha='center', va='center', fontweight='bold', fontsize=12)
    
    # Draw domains organized by type
    for domain_type, domains_of_type in domain_groups.items():
        if domain_type not in track_positions:
            continue
            
        y_pos = track_positions[domain_type]
        color = type_colors.get(domain_type, '#999999')
        domain_height = 0.4
        
        for domain in domains_of_type:
            start = domain['start']
            end = domain['end']
            length = end - start + 1
            description = domain.get('description', '')
            
            # Draw domain rectangle
            ax.add_patch(patches.Rectangle((start, y_pos), length, domain_height,
                                           facecolor=color, edgecolor='black', alpha=0.8, linewidth=1))
            
            # Add labels only for significant domains or those with space
            should_label = (
                length > 50 or  # Large domains
                domain_type in ['Domain', 'Repeat'] or  # Important domain types
                (domain_type == 'Transmembrane' and 'S1' in description)  # First transmembrane of each repeat
            )
            
            if should_label:
                # Determine label text
                if domain_type == 'Repeat':
                    label_text = f"Repeat {description}"
                elif domain_type == 'Domain':
                    label_text = f"{description} Domain"
                elif domain_type == 'Transmembrane' and 'repeat' in description:
                    # Extract repeat info for major transmembrane segments
                    label_text = description.split(';')[1].strip() if ';' in description else description
                    if 'Name=' in label_text:
                        label_text = label_text.replace('Name=', '')
                else:
                    label_text = description[:15] + '...' if len(description) > 15 else description
                
                # Position text
                text_x = start + length/2
                text_y = y_pos + domain_height/2
                
                # Choose text size and rotation based on domain size
                if length > 200:
                    fontsize = 9
                    rotation = 0
                elif length > 100:
                    fontsize = 8
                    rotation = 0
                else:
                    fontsize = 7
                    rotation = 45
                
                ax.text(text_x, text_y, label_text, ha='center', va='center', 
                       fontsize=fontsize, rotation=rotation, fontweight='bold')
    
    # Add track labels on the left
    for domain_type, y_pos in track_positions.items():
        if domain_type in domain_groups:
            count = len(domain_groups[domain_type])
            ax.text(-50, y_pos + 0.2, f"{domain_type}\n({count})", ha='right', va='center', 
                   fontweight='bold', fontsize=10, 
                   bbox=dict(boxstyle="round,pad=0.3", facecolor=type_colors.get(domain_type, '#999999'), alpha=0.3))
    
    # Set plot properties
    ax.set_xlim(-100, protein_length + 100)
    ax.set_ylim(0, 10)
    ax.set_xlabel('Amino Acid Position', fontsize=12, fontweight='bold')
    ax.set_title(f'SCN1A (P35498) Protein Domain Architecture', fontsize=16, fontweight='bold', pad=20)
    
    # Remove y-axis ticks and labels
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    # Add major position markers
    major_ticks = range(0, protein_length + 1, 500)
    ax.set_xticks(major_ticks)
    ax.grid(True, axis='x', alpha=0.3, linestyle='--')
    
    # Create legend
    legend_elements = []
    for domain_type in ['Repeat', 'Domain', 'Transmembrane', 'Intramembrane', 'Topological domain', 'Region']:
        if domain_type in domain_groups:
            color = type_colors.get(domain_type, '#999999')
            count = len(domain_groups[domain_type])
            legend_elements.append(patches.Patch(facecolor=color, label=f"{domain_type} ({count})"))
    
    if legend_elements:
        ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=10)
    
    plt.tight_layout()
    
    # Save the plot
    out_path = out_dir / f"{accession}_domains.png"
    try:
        plt.savefig(out_path, dpi=300, bbox_inches='tight')
        plt.close()
        logger.info(f"Saved domain visualization for {accession} to {out_path}")
        return out_path
    except Exception as e:
        logger.error(f"Error saving domain visualization for {accession}: {e}")
        plt.close()
        return None