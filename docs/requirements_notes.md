# Comparative Genomics Pipeline for Epilepsy-Related Ion Channels

## Project Overview

This pipeline enables modular, reproducible comparative genomics analysis of epilepsy-associated ion channel genes across a broad range of species. It integrates public gene/protein datasets, multiple sequence alignment, phylogenetic tree building, conservation scoring, human variant mapping, and visualization. The project is under active development, with a focus on expanding species coverage, improving analysis, and supporting scientific interpretation.

## Current Capabilities

- Fetches orthologous protein sequences for selected genes (**SCN1A**, **DEPDC5**) from UniProt and NCBI for human, mouse, and monkey (with easy extension to more species).
- Performs multiple sequence alignment (MSA) using EBI Clustal Omega API.
- Builds phylogenetic trees from unaligned FASTA files.
- Computes conservation scores (Shannon entropy, consensus) for each MSA and outputs results as CSV and PNG.
- Fetches and maps human disease variants from UniProt onto conservation plots.
- Fetches 3D protein structures (PDB) for selected genes.
- Modular design: easy to add new genes, species, or analysis steps.
- Markdown documentation and progress logs for scientific transparency.

## Planned/Upcoming Features

- Generalize variant/conservation overlay for all genes and species.
- Add/expand unit tests, code linting, and static analysis for code quality.
- Expand and polish documentation for public/academic use.
- Add support for more model and non-model organisms (e.g., reptiles, amphibians, fish, invertebrates).
- Integrate 3D structure visualization and mapping of variants onto structures.
- Enable cloud/DevOps support (e.g., AWS, Docker, Nextflow) for scalable, reproducible runs.
- Summarize research findings, plan further research, and write scientific conclusions.

## Requirements

- Python 3.10+
- pandas
- biopython
- matplotlib
- httpx
- (optional) jupyter, seaborn, scipy, requests, tqdm

Install with:

```bash
pip install -e .
```

## Usage Notes

- Not yet battle-tested for public use—expect quirks and rough edges. More details, tips, and polish coming soon.
- To clear all pipeline output data:

```bash
rm -f ./data/output/*/*
```

- Input species/genes are currently configured in `data/input/genes_to_proteins.json`. To focus on a subset (e.g., human, mouse, monkey), edit this file accordingly.
- All data and results are under active analysis; interpretation is ongoing.

## Scientific Context

The pipeline focuses on comparative genomics of epilepsy-associated genes, currently **SCN1A** and **DEPDC5**, across mammals (human, mouse, monkey). This enables insight into evolutionary conservation and divergence, and may help identify when seizure susceptibility first emerged in vertebrate evolution.

---

## 🧠 Phase 1: Define Scope & Collect Data

### ✅ Step 1. Define Research Scope

**Ion Channel Gene:**
A gene that encodes a protein forming an ion channel — a pore in a cell membrane that selectively allows charged ions (such as sodium, potassium, calcium, or chloride) to pass in and out of cells. Ion channels are essential for generating electrical signals in nerve and muscle cells.

**Neuronal Excitability:**
The ability of a neuron to fire an electrical signal (action potential) in response to stimulation. It reflects how easily a neuron "gets excited" to send messages. Neuronal excitability depends on the balance of ion flow across the neuron’s membrane through ion channels.

**Conserved Region:**
A segment of a protein or gene sequence that remains highly similar or identical across different species over evolutionary time, indicating its important functional or structural role.

**Divergent Region:**
A segment of a protein or gene sequence that shows significant variation between species, often reflecting adaptations, less critical function, or evolutionary change.

> *Early question:*
> If someone is epileptic, is there a mutation in a conserved or divergent region, or both, and why?

---

## Epilepsy-Related Genes Under Study

| Gene Symbol | Protein Type                                   | UniProt Accession | Notes                                                             |
| ----------- | ---------------------------------------------- | ----------------- | ----------------------------------------------------------------- |
| **SCN1A**   | Voltage-gated sodium channel α subunit         | `P35498`          | Key gene linked to Dravet syndrome and other epilepsy forms.      |
| **DEPDC5**  | GATOR1 complex subunit                         |                   | Mutations cause familial focal epilepsy; mTOR pathway regulator.  |

---

### ✅ Step 2. Collect Orthologous Protein Sequences

> *Question:*
> I need to select species other than humans to compare with. How do I decide which species to select? Should I just select randomly?

> *Answer:*
> Do **not** select randomly — we want to select for biologically meaningful coverage across the evolutionary tree.

**Guidelines for species selection:**

* Span vertebrates and invertebrates
* Vary in nervous system complexity
* Have good sequence availability

#### Example group-level selection (future expansion):

* Mammals
* Birds
* Reptiles
* Amphibians
* Bony Fish
* Insects
* Nematodes
* Cephalopods
* Jawless Vertebrates

---

**Orthologous Genes:**
Genes that evolved from a common ancestral gene by speciation and generally perform the same function.

> *Question:*
> With the diversity of groups we've selected (insects, nematodes, etc.), are we sure that the epileptic genes we've selected will be found in all of them, e.g., will it be an orthologous gene in this group?

> *Answer:*
> No, it is possible these genes won't be present in some groups. However, knowing which groups lack these genes can provide insight into the research.

**Keep in mind when selecting species:**

1. **Start with well-studied model organisms:**

   * Well-annotated genomes and proteomes available
   * Widely used in research with abundant data and tools
   * Representative of the taxonomic group and relevant to neuroscience/evolution

2. **Check ortholog availability:**

   * Use NCBI Gene, Ensembl, OrthoDB to confirm gene orthologs exist
   * Avoid species with missing or poorly annotated orthologs

3. **Consider evolutionary distance and biological relevance:**

   * Include both close relatives (fine differences) and distant relatives (to find conserved regions)
   * Select species with interesting traits linked to your research question

4. **Verify data quality:**

   * Prefer species with high-quality reference genomes and protein data
   * Avoid species with incomplete or low-quality sequence data

---

### Final selection (species level):

* Sourced from [NCBI](https://www.ncbi.nlm.nih.gov/)
* Choosing proteins for later Multiple Sequence Alignment (MSA) — proteins are more "obviously" informative than DNA sequences for this research project

---

## Mammals Under Study

### Human (*Homo sapiens*)

#### SCN1A
* **UniProt accession:** P35498
* **Entrez gene ID:** 6323
* **Entrez protein accession:** NP_001159435.1
* **Ensembl gene ID:** ENSG00000144285
* **Ensembl protein ID:** ENSP00000375401

#### DEPDC5
* **UniProt accession:** Q8TB45
* **Entrez gene ID:** 23344
* **Entrez protein accession:** NP_001138350.1
* **Ensembl gene ID:** ENSG00000198715
* **Ensembl protein ID:** ENSP00000369497

### Mouse (*Mus musculus*)

#### SCN1A
* **UniProt accession:** Q8VHK0
* **Entrez gene ID:** 20264
* **Entrez protein accession:** NP_001240789.1
* **Ensembl gene ID:** ENSMUSG00000075301
* **Ensembl protein ID:** ENSMUSP00000097412

#### DEPDC5
* **UniProt accession:** Q8BGT7
* **Entrez gene ID:** 233876
* **Entrez protein accession:** NP_001074923.1
* **Ensembl gene ID:** ENSMUSG00000037780
* **Ensembl protein ID:** ENSMUSP00000039313

### Monkey (*Macaca mulatta*)

#### SCN1A
* **UniProt accession:** F7GZC2
* **Entrez gene ID:** 100426356
* **Entrez protein accession:** XP_001104495.2
* **Ensembl gene ID:** ENSMMUG00000019913
* **Ensembl protein ID:** ENSMMUP00000026013

#### DEPDC5
* **UniProt accession:** F7GJQ2
* **Entrez gene ID:** 105595066
* **Entrez protein accession:** XP_015002345.1
* **Ensembl gene ID:** ENSMMUG00000020713
* **Ensembl protein ID:** ENSMMUP00000027113

---

## Project Goals & Aims
- Build a modular, reproducible comparative genomics pipeline for epilepsy-related genes (**SCN1A**, **DEPDC5**)
- Integrate public gene/protein datasets, MSA, phylogenetic tree building, conservation analysis, and human variant mapping
- Support analysis across a broad range of species (currently: human, mouse, monkey; easily extendable)
- Enable robust, scalable, and cloud-ready workflows (Docker, AWS, Nextflow)
- Output standard formats (FASTA, CSV, PNG, PDB, markdown) for scientific analysis and reproducibility
- Document all steps, code, and scientific reasoning for transparency and collaboration
- Plan for future expansion: more species, 3D structure mapping, cloud/DevOps, and deeper scientific interpretation

## Technical Requirements & Cloud/DevOps (Planned & In Progress)

### Core Software & Libraries
- Python >=3.10
- pandas
- biopython
- matplotlib
- httpx
- (optional) jupyter, seaborn, scipy, requests, tqdm

### Data & File Structure
- Input: `data/input/genes_to_proteins.json` (configurable for any species set)
- Output: All results (FASTA, MSA, trees, CSV, PNG, PDB) in `data/output/` subfolders
- Results and analysis: Markdown and images in `data/results/`

### Cloud/DevOps & AWS (Planned/Prototype)
- Pipeline designed for modularity and reproducibility, with future support for:
  - **AWS S3**: Store input/output data, results, and logs in S3 buckets for scalable, cloud-based workflows
  - **AWS Batch or EC2**: Run compute-intensive steps (e.g., MSA, tree building) on cloud VMs or batch jobs
  - **AWS Lambda**: Trigger pipeline steps or automate data movement
  - **Docker**: Containerize the pipeline for reproducible local/cloud execution
  - **Nextflow**: Orchestrate multi-step workflows, enable parallelization, and support cloud backends (including AWS Batch)
  - **CloudFormation/Terraform**: Infrastructure as code for reproducible AWS setup
- Credentials and configuration managed via environment variables or config files (never hardcoded)
- Logging and error handling designed for cloud-scale runs

### Example AWS Usage (Planned)
- Upload input FASTA/MSA files to S3
- Trigger pipeline run via Lambda or Nextflow on AWS Batch
- Store all outputs (CSV, PNG, PDB, logs) in S3 for downstream analysis
- Use Docker images for all compute steps to ensure reproducibility

---

Clustal Omega (Multiple Sequence Alignment)
- See https://www.ebi.ac.uk/jdispatcher/msa/clustalo for example
- Pass concatenated FASTA sequence
- Web API: https://www.ebi.ac.uk/jdispatcher/docs/webservices/#openapi

Phylogenetic Tree Generation
- PHYLIP or FASTA often preferred (Clustal Omega can produce these)

---

## 2025-06-27: Pipeline Focus Update

- For now, the analysis is focused on two genes (**SCN1A** and **DEPDC5**) and three mammalian species (human, mouse, monkey) to ensure data quality and realistic comparative genomics. This temporary scaling down is for robust validation and to simplify initial testing. Broader species and gene sets can be reintroduced after this minimal set is fully validated.

- Update: Only **SCN1A** and **DEPDC5** are present in `data/input/genes_to_proteins.json`, with entries for Homo sapiens and Mus musculus. Add Macaca mulatta (monkey) as soon as valid orthologs are confirmed.

- Remove checklist items for reptiles and amphibians for now, as the dataset is mammals-only and focused on **SCN1A** and **DEPDC5**.

## Ortholog ID Confirmation Checklist

For each gene and species in `data/input/genes_to_proteins.json`, confirm:
- [ ] Species name is correct and matches the intended model organism
- [ ] UniProt ID is valid and fetches the correct protein
- [ ] Entrez protein ID (if present) is valid and fetches the correct protein
- [ ] Sequence length and description match expectations for the gene
- [ ] No duplicate or missing entries for any gene/species pair

Record any issues or corrections needed below each entry.