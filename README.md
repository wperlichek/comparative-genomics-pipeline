## 🧬 Epilepsy Ion Channel Comparative Genomics Pipeline

Explore how epilepsy-related ion channel genes have evolved across species using foundational comparative genomics tools.

> 🧪 *Built as part of my transition into computational biology, this pipeline integrates public gene/protein datasets with practical bioinformatics techniques to enable early-stage, cross-species analysis of epilepsy-associated genes.*

### 🚀 Pipeline Flow  
`Fetch orthologs → Align sequences → Build phylogenetic trees → Score conservation → Map human variants → (Optional) Visualize in 3D → Store & run on AWS`

📁 **Requirements:** [docs/requirements.md](docs/requirements.md)

### Current Status

6/24/25 —  
Collected orthologous protein sequences for 5 key genes using NCBI and UniProt APIs; saved in `data/output/orthologs`.
Ran multiple sequence alignments (MSAs); results are in `data/output/msa`.  
Next up: figuring out how to actually interpret these alignments—reading up on best practices for MSA analysis before moving on to phylogenetic tree building.