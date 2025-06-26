# Comparative Genomics Pipeline Requirements

> This pipeline is designed for comparative genomics of epilepsy-related ion channel genes, but is modular and extensible for other gene families and research questions.

## Data Sources & Outputs
- **Data Sources/APIs:**
    - UniProtKB: protein sequences, human variant annotations
    - NCBI Entrez: additional protein orthologs
    - EBI Clustal Omega: multiple sequence alignment, phylogenetic tree generation
    - RCSB PDB: 3D protein structure files
- **Outputs:**
    - FASTA (orthologs, MSA), CSV (conservation scores), PNG (plots), PDB (structures), markdown (results/notes)
    - All outputs organized in `data/output/` and `data/results/` subfolders

## Project Goals & Aims
- Build a modular, reproducible comparative genomics pipeline for epilepsy-related ion channel genes
- Integrate public gene/protein datasets, MSA, phylogenetic tree building, conservation analysis, and human variant mapping
- Support analysis across a broad range of species (currently: human, rabbit, chicken, lamprey; easily extendable)
- Enable robust, scalable, and cloud-ready workflows (Docker, AWS, Nextflow)
- Output standard formats (FASTA, CSV, PNG, PDB, markdown) for scientific analysis and reproducibility
- Document all steps, code, and scientific reasoning for transparency and collaboration
- Plan for future expansion: more species, 3D structure mapping, cloud/DevOps, and deeper scientific interpretation

## Core Requirements
- Python >=3.10
- pandas
- biopython
- matplotlib
- httpx

## Optional/Recommended
- jupyter
- seaborn
- scipy
- requests
- tqdm

## Pipeline Steps
1. Fetch orthologous protein sequences for selected genes and species
2. Align sequences using EBI Clustal Omega (MSA)
3. Build phylogenetic trees from unaligned FASTA files
4. Compute conservation scores (Shannon entropy, consensus) for each MSA
5. Map and visualize human disease variants on conservation plots
6. Fetch and save 3D protein structures (PDB) for selected genes
7. (Optional) Visualize in 3D
8. (Optional) Store & run on AWS/cloud
9. Clean and unit test all code
10. Code analysis (linting, type checks, static analysis)
11. Expand and polish documentation
12. Summarize research findings, plan further research, and write conclusions

## Installation

```bash
pip install -e .
```

## Usage Notes
- Not yet battle-tested for public use—expect quirks and rough edges. More details, tips, and polish coming soon.
- To clear all pipeline output data:
  ```bash
  rm -f ./data/output/*/*
  ```
- Input species/genes are currently configured in `data/input/genes_to_proteins.json`. To focus on a subset (e.g., human, rabbit, chicken, lamprey), edit this file accordingly.
- All data and results are under active analysis; interpretation is ongoing.

## Project Status
- Modular pipeline for ortholog fetching, MSA, tree building, conservation scoring, variant mapping, and visualization is complete for human, rabbit, chicken, and lamprey.
- Documentation, code quality, and scientific analysis are being expanded.
- Future plans: more species, 3D structure mapping, cloud/DevOps support, and deeper scientific interpretation.

## Cloud & Containerization Roadmap (AWS, Docker, Nextflow)

### 1. Containerize the Pipeline with Docker
- [ ] Write a `Dockerfile`:
    - Use `python:3.10-slim` as the base image
    - Copy all source code and requirements into the image
    - Run `pip install -e .` during build
    - Set the entrypoint to your main script (e.g., `comparative-genomics-pipeline`)
- [ ] Build and test locally:
    - `docker build -t genomics-pipeline .`
    - `docker run -v $(pwd)/data:/app/data genomics-pipeline`

### 2. Push Docker Image to AWS ECR
- [ ] Create an ECR repository in AWS
- [ ] Authenticate Docker to ECR (`aws ecr get-login-password`)
- [ ] Tag and push your image to ECR
    - `docker tag genomics-pipeline:latest <your-account-id>.dkr.ecr.<region>.amazonaws.com/genomics-pipeline:latest`
    - `docker push <your-account-id>.dkr.ecr.<region>.amazonaws.com/genomics-pipeline:latest`

### 3. Set Up AWS S3 for Data Storage
- [ ] Create S3 buckets for input, output, and logs
- [ ] Upload your input files (e.g., FASTA, JSON) to S3
- [ ] Plan to sync results back to S3 after each run

### 4. Configure AWS Batch for Scalable Compute
- [ ] Create a compute environment and job queue in AWS Batch
- [ ] Define a job definition using your ECR image
- [ ] Set up IAM roles for Batch and S3 access
- [ ] Submit jobs to Batch for each pipeline step (MSA, tree, conservation, etc.)
- [ ] Use environment variables for S3 paths and config

### 5. Orchestrate with Nextflow (Optional, for Multi-step Pipelines)
- [ ] Write a `main.nf` Nextflow script to define pipeline steps as processes
- [ ] Configure each process to use your Docker image and S3 for input/output
- [ ] Write a `nextflow.config` for AWS Batch profile
- [ ] Run with: `nextflow run main.nf -profile awsbatch`

### 6. Infrastructure as Code (CloudFormation/Terraform)
- [ ] Write templates to define:
    - S3 buckets
    - ECR repository
    - Batch compute environment, job queue, job definitions
    - IAM roles and policies
- [ ] Store templates in an `infra/` or `aws/` folder in your repo

### 7. Best Practices & Automation
- [ ] Never hardcode credentials; use IAM roles or environment variables
- [ ] Log all pipeline steps and errors to S3 or CloudWatch
- [ ] Use versioned S3 buckets for reproducibility
- [ ] Document all setup, deployment, and teardown steps in your repo
