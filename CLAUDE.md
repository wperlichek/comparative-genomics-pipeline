# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

**Installation & Execution:**
```bash
# Local development
pip install -e .
comparative-genomics-pipeline

# Docker (recommended for reproducibility)
docker build -t genomics-pipeline .
docker run --rm -v $(pwd)/data:/app/data genomics-pipeline

# Clear previous results
rm -rf ./data/output/*
```

**No testing framework is currently implemented** - this is a known limitation that needs to be addressed.

## Architecture Overview

This is a bioinformatics pipeline for comparative genomics analysis, specifically focused on epilepsy-related ion channel genes (primarily SCN1A) across vertebrate species. The pipeline performs ortholog retrieval, multiple sequence alignment, phylogenetic analysis, conservation scoring, and human variant mapping.

**Core Workflow:**
1. Ortholog Collection (UniProt/NCBI APIs)
2. Multiple Sequence Alignment (EBI Clustal Omega)
3. Phylogenetic Tree Generation
4. Conservation Analysis (Shannon entropy)
5. Variant Mapping
6. Scientific Visualization
7. Structure Retrieval (PDB)

**Tech Stack:** Python 3.10+ with asyncio, BioPython, httpx, matplotlib, pandas, scipy, Docker

## Code Organization

**Entry Point:** `src/comparative_genomics_pipeline/__main__.py` - Main async pipeline orchestration

**Key Directories:**
- `client/` - External API integrations (EBI, NCBI, UniProt, PDB)
- `service/` - Business logic layer (biopython_service.py)
- `config/` - Configuration management (paths, logging, constants)
- `visualization/` - Scientific plotting with publication-quality output
- `util/` - File I/O utilities

**Data Flow:** 
- Input: `data/input/genes_to_proteins.json` (species and UniProt ID mappings)
- Output: `data/output/` (FASTA, CSV, PNG, markdown results)

## Development Patterns

**Async/Await Architecture:** Heavy use of asyncio for concurrent API calls to external databases. The main pipeline in `__main__.py` orchestrates multiple async clients.

**Client-Service Pattern:** 
- Clients handle external API communications with proper error handling
- Services contain bioinformatics computations and business logic
- Clear separation between data retrieval and processing

**Configuration-Driven Design:** Centralized configuration in `config/` directory with structured input via JSON rather than command-line arguments.

## Current Limitations

- **No testing infrastructure** (no pytest, no unit tests, empty quality/ directory)
- **Limited error handling** for API failures and network issues
- **Hardcoded configuration** for specific gene/species combinations
- **No automated validation** of bioinformatics outputs
- **Proof-of-concept scale** (currently 5 species, 1 primary gene)

## External Dependencies

**APIs:** UniProtKB, NCBI Entrez, EBI Clustal Omega, RCSB PDB - requires internet connectivity

**Data Sources:** Pipeline fetches live data from multiple genomic databases, results are cached locally in file system structure

## Key Files for Context

- `data/input/genes_to_proteins.json` - Species and gene configuration
- `src/comparative_genomics_pipeline/config/path_config.py` - File path management
- `README.md` - Project overview with current results and roadmap
- `docs/requirements.md` - Detailed technical requirements