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

**API Format Preference:** Always prefer JSON over XML for API calls unless JSON is not available for a specific endpoint.

## Current Implementation Status

**Project Maturity Level: Mid-stage prototype (60-70% complete)**

The pipeline is a functional bioinformatics research tool that has successfully generated scientific results for SCN1A and DEPDC5 epilepsy genes across 5 vertebrate species. Key implementation achievements:

**✅ Fully Implemented:**
- Complete 8-step async pipeline orchestration
- All external API integrations (UniProt, NCBI, EBI, ClinVar, PDB)
- AWS S3 caching for UniProt FASTA sequences with graceful fallback
- Multiple sequence alignment and phylogenetic tree generation
- Conservation scoring with Shannon entropy analysis
- Clinical variant mapping and visualization
- Scientific plotting with basic statistical analysis
- Docker containerization for reproducibility

**⚠️ Partially Implemented:**
- Test coverage at 18% (basic unit tests for core clients)
- Error handling with retry logic (basic implementation)
- Performance optimization (sequential gene processing)

**❌ Not Yet Implemented:**
- Automated CI/CD pipeline
- Advanced error recovery systems

## Current Limitations & Technical Debt

**Performance Bottlenecks:**
- **EBI Clustal Omega dependency**: External service wait times for alignment jobs
- **Sequential processing**: Genes processed one at a time rather than in parallel
- **Partial API caching**: FASTA sequences cached in S3, but variants and NCBI data still hit APIs repeatedly

**Infrastructure Gaps:**
- **Limited test coverage**: 18% coverage, significant gaps in service layer testing
- **Basic error handling**: Minimal retry logic, no circuit breakers
- **Partial S3 integration**: Only UniProt sequences cached, variants and NCBI data still uncached

**Research Scale Constraints:**
- **Species limitation**: Currently optimized for 5 vertebrate species
- **Gene scope**: Focused on 2 primary epilepsy genes (SCN1A, DEPDC5)
- **Storage approach**: Local file system only, no distributed storage

## External Dependencies

**APIs:** UniProtKB, NCBI Entrez, EBI Clustal Omega, RCSB PDB - requires internet connectivity

**Data Sources:** Pipeline fetches live data from multiple genomic databases, results are cached locally in file system structure

## Key Files for Context

- `data/input/genes_to_proteins.json` - Species and gene configuration
- `src/comparative_genomics_pipeline/config/path_config.py` - File path management
- `README.md` - Project overview with current results and roadmap
- `docs/requirements.md` - Detailed technical requirements
- **`docs/early_hypothesis_research.md` - Core research questions and hypothesis driving all analysis work**
- `docs/reference_papers.md` - Key research papers and citations supporting the analysis
- `docs/research_log.md` - Development progress and research notes

## Working Preferences & Collaboration Style

**Research Context Awareness:** ALWAYS consult `docs/early_hypothesis_research.md` before performing any analysis or making recommendations. All work should be framed within the context of the core research questions and hypothesis outlined in that document.

**Concise Communication:** Keep explanations brief and focused. Avoid unnecessary preamble or postamble unless detail is specifically requested.

**Scientific Precision:** When updating documentation or plots, ensure scientific accuracy in terminology. For example:
- Specify "protein phylogenetic trees" rather than generic "phylogenetic trees"
- Distinguish between different variant types (loss-of-function vs. likely pathogenic)
- Maintain consistency between documentation and code output

**Iterative Refinement:** Expect requests for consistency across multiple components:
- README documentation should match plot titles and descriptions
- Plot titles should reflect the underlying methodology (e.g., protein-based analysis)
- Update both documentation and visualization code when making conceptual changes

**Parallel Structure:** When adding content for multiple genes (SCN1A, DEPDC5), maintain parallel formatting and structure in documentation, using appropriate subheadings for organization.

**Implementation Details:** Be prepared to trace through the codebase to understand implementation specifics (e.g., whether trees use protein vs. DNA sequences) when making documentation updates.