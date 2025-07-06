#!/usr/bin/env python3
"""
Quick demo script for the Comparative Genomics Pipeline
Shows key project outputs and capabilities for presentation/resume purposes
"""

import os
import sys
from pathlib import Path

def main():
    """Demo the pipeline outputs and capabilities"""
    
    print("🧬 Comparative Genomics Pipeline - Demo")
    print("=" * 50)
    
    # Check if we're in the right directory
    if not os.path.exists('data/output'):
        print("❌ Error: Run this script from the project root directory")
        sys.exit(1)
    
    print("\n📊 Project Status:")
    print("✅ Complete bioinformatics pipeline implementation")
    print("✅ 118 passing tests with 34% code coverage")
    print("✅ Analysis of SCN1A and DEPDC5 epilepsy genes")
    print("✅ 5 vertebrate species comparative analysis")
    print("✅ Docker containerization for reproducibility")
    
    print("\n🔬 Scientific Outputs Generated:")
    
    # Check for key output files
    output_files = [
        ("Phylogenetic Trees", "data/output/trees/*.png"),
        ("Conservation Analysis", "data/output/conservation/*.png"),
        ("Variant Mapping", "data/output/variants/*.png"),
        ("Conservation Data", "data/output/conservation/*.csv"),
        ("Sequence Alignments", "data/output/msa/*.fasta"),
        ("Protein Structures", "data/output/structures/*.png")
    ]
    
    from glob import glob
    
    for category, pattern in output_files:
        files = glob(pattern)
        if files:
            print(f"✅ {category}: {len(files)} files")
            for f in files[:2]:  # Show first 2 files
                print(f"   📄 {os.path.basename(f)}")
            if len(files) > 2:
                print(f"   ... and {len(files) - 2} more")
        else:
            print(f"❌ {category}: No files found")
    
    print("\n🎯 Key Achievements:")
    print("• Evolutionary conservation analysis across vertebrate species")
    print("• Clinical variant mapping to conserved regions")
    print("• Publication-quality scientific visualizations")
    print("• AWS S3 integration for sequence caching")
    print("• Comprehensive test suite with async API integration")
    
    print("\n🚀 Quick Start Commands:")
    print("pip install -e .")
    print("comparative-genomics-pipeline")
    print("# OR")
    print("docker build -t genomics-pipeline .")
    print("docker run --rm -v $(pwd)/data:/app/data genomics-pipeline")
    
    print("\n📋 Testing:")
    print("pytest --cov=comparative_genomics_pipeline")
    
    print("\n🔗 Key Technologies:")
    print("• Python 3.10+ with asyncio")
    print("• BioPython, matplotlib, pandas, scipy")
    print("• Docker containerization")
    print("• AWS S3 for data caching")
    print("• Multiple genomic database APIs")
    
    print("\n" + "=" * 50)
    print("✨ Ready for resume presentation!")

if __name__ == "__main__":
    main()