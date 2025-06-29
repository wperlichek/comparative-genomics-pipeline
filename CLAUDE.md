# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## CRITICAL: Scientific Accuracy Mandate

**ABSOLUTE PRIORITY:** Scientific accuracy is the highest priority when working with this bioinformatics pipeline. Never make biological or scientific claims that are not directly supported by the data or analysis results. Always:

1. **Stick to the data:** Only describe what the pipeline actually measures/computes, not biological interpretations
2. **Avoid inferences:** Do not infer protein domains, functional regions, or biological mechanisms unless explicitly annotated in the data
3. **Distinguish analysis from biology:** Clearly separate computational results (e.g., "highly conserved regions") from biological claims (e.g., "transmembrane domains")
4. **Verify annotations:** Check that claimed features (domains, motifs, functional sites) are actually present in the pipeline output, not assumed from biological knowledge
5. **When uncertain:** Always err on the side of being more conservative and factual rather than more interpretive

This is a research tool - scientific integrity is non-negotiable.

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

**Testing Framework:**
```bash
# Install dev dependencies
pip install -e ".[dev]"

# Run all tests
pytest

# Run specific test categories
pytest -m unit          # Unit tests only
pytest -m integration   # Integration tests only
pytest -m "not slow"    # Skip slow tests

# Run with coverage
pytest --cov=comparative_genomics_pipeline --cov-report=html

# Run specific test file
pytest tests/unit/test_uniprot_client.py
```

## Architecture Overview

This is a bioinformatics pipeline for comparative genomics analysis, specifically focused on epilepsy-related genes SCN1A (Dravet syndrome) and DEPDC5 (focal epilepsy) across vertebrate species. The pipeline performs ortholog retrieval, multiple sequence alignment, phylogenetic analysis, conservation scoring, and human variant mapping.

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

- **Limited test coverage** (basic unit and integration tests implemented, but needs expansion)
- **Limited error handling** for API failures and network issues
- **Hardcoded configuration** for specific gene/species combinations
- **No automated validation** of bioinformatics outputs
- **Proof-of-concept scale** (currently 5 species, 2 primary genes: SCN1A and DEPDC5)

## External Dependencies

**APIs:** UniProtKB, NCBI Entrez, EBI Clustal Omega, RCSB PDB - requires internet connectivity

**Data Sources:** Pipeline fetches live data from multiple genomic databases, results are cached locally in file system structure

## Key Files for Context

- `data/input/genes_to_proteins.json` - Species and gene configuration
- `src/comparative_genomics_pipeline/config/path_config.py` - File path management
- `README.md` - Project overview with current results and roadmap
- `docs/requirements.md` - Detailed technical requirements