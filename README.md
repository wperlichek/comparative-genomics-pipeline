# Comparative Genomics Pipeline

A Python pipeline for analyzing evolutionary conservation of epilepsy-associated genes across vertebrate species, focusing on SCN1A (Dravet syndrome) and DEPDC5 (focal epilepsy).

**Current Status:** Work in progress. Core functionality implemented, testing and validation ongoing.

**Pipeline:** Ortholog retrieval → Multiple sequence alignment → Phylogenetic analysis → Conservation scoring → Human variant mapping

**Primary Genes:** SCN1A (voltage-gated sodium channel, Dravet syndrome) and DEPDC5 (mTOR pathway regulator, focal epilepsy)
**Species:** Human, mouse, macaque, chicken, great tit

## Implementation

**Data Sources:**
- UniProtKB: Protein sequences, variant annotations
- NCBI Entrez: Orthologous sequences
- EBI Clustal Omega: MSA and phylogenetic trees
- RCSB PDB: Protein structures

**Technologies:**
- Python 3.9+ with asyncio for concurrent API calls
- Biopython for sequence analysis
- Docker for containerization
- NumPy/matplotlib for data processing and visualization


## Current Functionality

**Working:**
- Complete 8-step async pipeline with proper orchestration
- Multi-database API integration (UniProt, NCBI, EBI Clustal Omega, ClinVar, PDB)
- AWS S3 caching for UniProt sequences with graceful fallback
- Multiple sequence alignment and phylogenetic tree generation
- Shannon entropy conservation scoring with statistical analysis
- Clinical variant mapping and visualization for both genes
- Publication-quality scientific plots with matplotlib
- CLI interface and Docker containerization
- Comprehensive input validation and file I/O error handling

**Limitations:**
- Research-scale scope (5 vertebrate species, 2 primary genes)
- Sequential gene processing (not parallelized)
- Basic API retry logic without exponential backoff or circuit breakers
- Partial S3 caching (sequences only, variants/NCBI data not cached)
- No CLI argument parsing (fixed configuration only)
- Limited progress feedback for long-running operations
- Output validation minimal compared to input validation

## 🧠 AI-Assisted Development

This project uses Claude AI to help with boilerplate code generation, plot refinements, documentation, and even entire feature generation. I believe AI makes it so we can focus more on the research rather than the code — experienced programmers can shepherd and review code and orchestrate an app rather than spend time on coding details that are sometimes important to dive into but often not.

## 🛠️ Installation & Usage

**Requirements:** Python 3.10+, Docker (recommended for reproducibility), and git.

### 1. Clone the repository
```bash
git clone https://github.com/wperlichek/comparative-genomics-pipeline.git
cd comparative-genomics-pipeline
```

### 2. Install Python dependencies (local development)
```bash
pip install -e .
```

### 3. Run the pipeline (local)
```bash
comparative-genomics-pipeline
```

### 4. Build and run with Docker (recommended)
```bash
docker build -t genomics-pipeline .
docker run --rm -v $(pwd)/data:/app/data genomics-pipeline
```
- The `-v $(pwd)/data:/app/data` flag mounts your local `data/` directory for input/output.

### 5. Clear previous results
```bash
rm -rf ./data/output/*
```

### 6. Run tests
```bash
pytest
```

---

- Input configuration: Edit `data/input/genes_to_proteins.json` to specify which genes and species to analyze.
- Results (plots, trees, CSVs) are saved in `data/output/` subfolders.
- For troubleshooting Docker, see the README and `docs/requirements_notes.md` for platform-specific tips.
- For cloud/AWS usage, see the DevOps section in the documentation.

## Development Roadmap

**Immediate Priorities:**
- Comprehensive testing suite (pytest)
- Robust error handling and API retry logic
- Configuration management (YAML-based)
- Command-line interface with proper argument parsing

**Planned Extensions:**
- Enhanced analysis of SCN1A and DEPDC5 variants
- Protein domain annotation integration
- Statistical significance testing for conservation scores
- Clinical variant database integration (ClinVar, OMIM)

## Current Analysis Outputs 🔄

*These outputs are generated from genomic databases with S3 caching for protein sequences*

### SCN1A

<img src="data/output/trees/SCN1A_scientific.png" width="600">

**Protein Phylogenetic Tree:** Evolutionary relationships confirming expected species divergence patterns.

<img src="data/output/conservation/SCN1A_conservation_scientific.png" width="800">

**Conservation Analysis:** SCN1A evolutionary conservation across 5 vertebrate species, with SCN1A showing 90.2% of positions highly conserved.

<img src="data/output/variants/SCN1A_conservation_variants_scientific.png" width="800">

**Conservation & Variants:** SCN1A conservation analysis with loss-of-function and likely pathogenic variants overlaid, showing variant distribution across conserved and variable regions.

### DEPDC5

<img src="data/output/trees/DEPDC5_scientific.png" width="600">

**Protein Phylogenetic Tree:** Evolutionary relationships confirming expected species divergence patterns.

<img src="data/output/conservation/DEPDC5_conservation_scientific.png" width="800">

**Conservation Analysis:** DEPDC5 evolutionary conservation across 5 vertebrate species.

<img src="data/output/variants/DEPDC5_conservation_variants_scientific.png" width="800">

**Conservation & Variants:** DEPDC5 conservation analysis with likely pathogenic variants overlaid, showing how mutations cluster in specific functional regions.

<img src="data/output/variants/clinvar_variants_comparison.png" width="800">

**ClinVar Variant Analysis:** Raw ClinVar data showing pathogenic and likely pathogenic variants for SCN1A and DEPDC5.

## Research Focus

**SCN1A (Dravet Syndrome):** Voltage-gated sodium channel predominantly expressed in GABAergic interneurons. Loss-of-function variants impair interneuron firing, reducing GABA release and causing network disinhibition leading to hyperexcitability and seizures.

**DEPDC5 (Focal Epilepsy):** DEP domain-containing protein 5, part of the GATOR1 complex that regulates mTOR signaling. Mutations cause focal cortical dysplasia and familial focal epilepsy with variable foci, affecting neural development and excitability through mTOR pathway dysregulation.

*For research motivations, see [docs/early_hypothesis_research.md](docs/early_hypothesis_research.md).*

*For ongoing research findings and insights, see [docs/research_log.md](docs/research_log.md).*

*For academic papers and literature references, see [docs/reference_papers.md](docs/reference_papers.md).*

---

*Self-directed computational biology project. Ongoing development and validation.*
