# Comparative Genomics Pipeline for Epilepsy-Related Ion Channels

## 🧠 Phase 1: Define Scope & Collect Data

### ✅ Step 1. Define Research Scope

* Focus genes: e.g., `SCN1A` (sodium channel), `KCNQ2` (potassium channel).
* Research Question:
  *"How have key epilepsy-related ion channel genes evolved across species, and are there conserved or divergent regions that may relate to neuronal excitability and epilepsy risk?"*

### ✅ Step 2. Collect Orthologous Protein Sequences

* Select \~5–10 species: e.g., **human, mouse, rat, zebrafish, Drosophila, C. elegans**.
* Use **Biopython** or **UniProt REST API** to fetch protein sequences.
* Save sequences as **FASTA**, with consistent naming.
* Handle missing or failed API responses gracefully.

---

## 🧬 Phase 2: Build the Python-Based Computational Pipeline

### ✅ Step 3. Align Sequences (MSA)

* Use external tool (e.g., **MAFFT**, **Clustal Omega**).
* Python script should:

  * Call MSA tool via `subprocess`.
  * Pass FASTA as input.
  * Save alignment result in standard format (e.g., **CLUSTAL**, **FASTA**).

### ✅ Step 4. Build Phylogenetic Tree

* Use external tool (e.g., **IQ-TREE**, **FastTree**).
* Python script should:

  * Call tree-building tool using aligned sequences.
  * Save output as **Newick** format.

### ✅ Step 5. Analyze Conservation

* Compute **conservation scores** across aligned sequences:

  * Use entropy or simple consensus-based scoring.
* Output: table of conservation scores by alignment position.

### ✅ Step 6. Map Known Human Variants

* Search **ClinVar**, **UniProt**, or other public databases for epilepsy-causing variants.
* Map variants to the protein sequence.
* Visualize:

  * Alignment with conserved regions highlighted.
  * Variant positions overlaid.

### ✅ Step 7. Optional: 3D Protein Mapping (Stretch Goal)

* Use **PDB API** to fetch available protein structure (if any).
* Identify structural domains where conserved/variant regions fall.

---

## ☁️ Phase 3: AWS Cloud Integration

### ✅ Step 8. Use AWS S3 for Storage

* Create versioned **S3 buckets**:

  * `raw_sequences/`, `alignments/`, `trees/`, `results/`
* Upload and retrieve data programmatically.

### ✅ Step 9. Containerize Tools with Docker

* Create **Dockerfiles** for:

  * Sequence analysis scripts.
  * MSA tool.
  * Tree-building tool.

### ✅ Step 10. Run on AWS Batch

* Create **AWS Batch jobs** for each pipeline stage.
* Configure **compute environments** (e.g., EC2-based).
* Dynamically allocate resources (stretch goal).

### ✅ Step 11. Orchestrate with Nextflow (Stretch)

* Write a **Nextflow pipeline** to run all stages:

  * Each step: Docker container → reads from S3 → writes to S3.
  * Track params, logs, and provenance.

---

## 🛠️ Phase 4: Infra-as-Code & Deployment

### ✅ Step 12. Write CloudFormation Template

* Template provisions:

  * S3 buckets
  * IAM roles/policies
  * AWS Batch compute environments
  * Job queues
* Ensure teardown & redeployability.

---

## 📆 Phase 5: Deliverables

### ✅ Step 13. GitHub Repo Structure

* Organize:

  * `/scripts/` – Python code
  * `/docker/` – Dockerfiles
  * `/nextflow/` – Workflow DSL (optional)
  * `/cloudformation/` – Infra-as-code
  * `/results/` – Example outputs (plots, trees)
  * `README.md`
  * `requirements.txt`

### ✅ Step 14. Documentation

* In `README.md`:

  * State research question and biological motivation
  * Describe pipeline architecture
  * Include setup instructions for AWS env
  * Explain how to run pipeline and interpret results
  * Link to public data sources

### ✅ Step 15. Results Summary

* Tables: conservation scores per gene/position.
* Plots:

  * Conservation vs. position
  * Phylogenetic tree
  * Annotated MSA with variant overlay
* Write a short interpretation of findings.