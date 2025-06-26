# Comparative Genomics Pipeline for Epilepsy-Related Ion Channels

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

## Top 5 Epilepsy-Related (Human) Genes to Study

| Gene Symbol | Protein Type                                   | UniProt Accession | Notes                                                             |
| ----------- | ---------------------------------------------- | ----------------- | ----------------------------------------------------------------- |
| **SCN1A**   | Voltage-gated sodium channel α subunit         | `P35498`          | Key gene linked to Dravet syndrome and other epilepsy forms.      |
| **KCNQ2**   | Voltage-gated potassium channel subunit        | `O43526`          | Involved in neonatal epilepsy; regulates neuronal excitability.   |
| **SCN2A**   | Voltage-gated sodium channel α subunit         | `Q99250`          | Associated with early-onset epilepsy and developmental issues.    |
| **CACNA1H** | T-type voltage-gated calcium channel           | `O95180`          | Implicated in generalized and childhood absence seizures.         |
| **GABRA1**  | GABA\_A receptor subunit α1 (chloride channel) | `P14867`          | Important for inhibitory signaling; mutations linked to epilepsy. |


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

#### Example group-level selection:

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

## Mammal: Rabbit (*Oryctolagus cuniculus*)

### UniProt data for orthologs:

#### SCN1A

* **UniProt accession:** G1SSP8
* **Entrez gene ID:** 100009591
* **Entrez protein accession:** XP\_051704506.1
* **Ensembl gene ID:** ENSOCUG00000007266.4
* **Ensembl protein ID:** ENSOCUP00000006283.2

#### KCNQ2

* **Entrez gene ID:** 127485776
* **Entrez protein accession:** XP\_051685089.1
* **NCBI fallback:** https://www.ncbi.nlm.nih.gov/protein/XP_051685089.1 

#### SCN2A

* **UniProt accession:** A0A5F9C773
* **Entrez gene ID:** 100351658
* **Entrez protein accession:** XP\_008256914.1
* **Ensembl gene ID:** ENSOCUG00000005533.4
* **Ensembl protein ID:** ENSOCUP00000025500.1; ENSOCUP00000029697.1

#### CACNA1H

* **Entrez gene ID:** 127486418
* **Entrez protein:** XP_051686345.1
* **NCBI fallback:** https://www.ncbi.nlm.nih.gov/gene/127486418

#### GABRA1

* **UniProt**: G1TAX4
* **Entrez gene**: 100357047
* **Entrez protein**: XP_002710402.2
* **Ensembl**: ENSOCUG00000016109.4
* **Ensembl protein**: ENSOCUP00000013846.3

---

### Reasoning for selecting Rabbit:

* ☑ Species has a publicly available genome and proteome

  * \~94,000 protein sequences, \~22,000 protein-coding genes
  * Source: [https://www.ncbi.nlm.nih.gov/datasets/genome/GCF\_964237555.1/](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_964237555.1/)
* ☑ High-quality genome annotation

  * Level: Chromosome (Top tier)
  * Annotation: NCBI RefSeq (Gold standard)
  * Genes: 38,822
  * Protein-coding: 21,895 (robust coverage and annotation depth)
* ☑ Species is commonly used in research (e.g., neuroscience, evolution)

  * PubMed hits: 228 (as of 6/23/25) for "Oryctolagus cuniculus"
* ☑ Species represents its group well (classic model organism)

  * Not an outlier; divergent enough from human without being too divergent
* ☑ Ortholog presence confirmed via NCBI Gene, Ensembl, OrthoDB

  * Confirmed ortholog presence for rabbit
* ☑ Avoids missing or poorly annotated orthologs
* ☑ Species has interesting traits linked to the research question

  * Literature shows rabbits experience spontaneous seizures from diverse causes
  * Veterinary reviews suggest rabbits are suitable for seizure research, with measurable episodes and often benign recovery

---

## Mammal: Chicken ("Gallus gallus")

## Reasoning: Chicken is the most model big organism.

### UniProt data for orthologs:

#### SCN1A

* **UniProt accession:** A0A8V1A1D3
* **Entrez gene ID:** SCN1A
* **Entrez protein accession:** XP\_004942836.4
* **Ensembl gene ID:** ENSGALG00010026770.1
* **Ensembl protein ID:** ENSGALP00010039668.1

#### KCNQ2

* **UniProt accession:** A0A8V0ZLH7
* **Entrez gene ID:** KCNQ2
* **Entrez protein accession:** XP\_040506889.1
* **Ensembl gene ID:** ENSGALG00010023174.1
* **Ensembl protein ID:** ENSGALP00010034224.1

#### SCN2A

* **UniProt accession:** A0A8V1A5Z8
* **Entrez gene ID:** SCN2A
* **Entrez protein accession:** XP\_040559277.1
* **Ensembl gene ID:** ENSGALG00010025597.1
* **Ensembl protein ID:** ENSGALP00010038628.1

#### CACNA1H

* **UniProt accession:** A0A8V0Z849
* **Entrez gene ID:** CACNA1H
* **Entrez protein accession:** XP\_015149910.3
* **Ensembl gene ID:** ENSGALG00010017430.1
* **Ensembl protein ID:** ENSGALP00010024632.1

#### GABRA1

* **UniProt**: A0A8V0YC88
* **Entrez gene**: GABRA1
* **Entrez protein**: XP\_046756226.1
* **Ensembl**: ENSGALG00010009841
* **Ensembl protein**: ENSGALP00010013552

---

## Expanded Group-Level Selection

* Mammals: Human (*Homo sapiens*) — all five genes present, well-studied
* Birds: Chicken (*Gallus gallus*) — classic model, all five genes annotated
* Reptiles: Green anole (*Anolis carolinensis*) — robust UniProt/NCBI/Ensembl coverage
* Amphibians: Western clawed frog (*Xenopus tropicalis*) — classic model, all five genes annotated
* Bony Fish: Zebrafish (*Danio rerio*) — gold-standard, all five genes present
* Insects: Fruit fly (*Drosophila melanogaster*) — all five genes have clear orthologs
* Nematodes: *Caenorhabditis elegans* — all five genes have orthologs
* Cephalopods: California two-spot octopus (*Octopus bimaculoides*) — genome/proteome available, some ion channel orthologs present
* Jawless Vertebrates: Sea lamprey (*Petromyzon marinus*) — NCBI/Ensembl annotated proteome, some ion channel orthologs present. The lamprey has been extensively studied because its relatively simple brain is thought in many respects to reflect the brain structure of early vertebrate ancestors.

---

Clustal Omega (Multiple Sequence Alignment)
- See https://www.ebi.ac.uk/jdispatcher/msa/clustalo for example
- Pass concataneted FASTA sequence
- Web API: https://www.ebi.ac.uk/jdispatcher/docs/webservices/#openapi

Phylogenic Tree Generation
- PHYLIP or FASTA often preferred (Clustal Omega can produce these)