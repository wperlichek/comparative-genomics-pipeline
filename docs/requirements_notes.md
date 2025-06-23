# Comparative Genomics Pipeline for Epilepsy-Related Ion Channels

## 🧠 Phase 1: Define Scope & Collect Data

### ✅ Step 1. Define Research Scope

**Ion Channel Gene:**
A gene that encodes a protein forming an ion channel — a pore in a cell membrane that selectively allows charged ions (such as sodium, potassium, calcium, or chloride) to pass in and out of cells. Ion channels are essential for generating electrical signals in nerve and muscle cells.

**Neuronal Excitability**
The ability of a neuron to fire an electrical signal (action potential) in response to stimulation. It reflects how easily a neuron "gets excited" to send messages. Neuronal excitability depends on the balance of ion flow across the neuron’s membrane through ion channels.

**Converved Region**
A segment of a protein or gene sequence that remains highly similar or identical across different species over evolutionary time, indicating its important functional or structural role.

**Divergent Region**
A segment of a protein or gene sequence that shows significant variation between species, often reflecting adaptations, less critical function, or evolutionary change.

_Early question_: If someone is epileptic, is there a mutation in a conserved or divergent region, or both, and why?

## Top 5 Epilepsy-Related Genes to Study

| Gene Symbol | Protein Type                             | Notes                                                                                 |
|-------------|----------------------------------------|---------------------------------------------------------------------------------------|
| **SCN1A**   | Voltage-gated sodium channel α subunit | Key gene linked to Dravet syndrome and other epilepsy forms.                           |
| **KCNQ2**   | Voltage-gated potassium channel subunit| Involved in neonatal epilepsy; regulates neuronal excitability.                       |
| **SCN2A**   | Voltage-gated sodium channel α subunit | Associated with early-onset epilepsy and developmental issues.                        |
| **CACNA1H** | T-type voltage-gated calcium channel   | Implicated in generalized and childhood absence seizures.                            |
| **GABRA1**  | GABA_A receptor subunit α1 (chloride channel) | Important for inhibitory signaling; mutations linked to epilepsy.                    |

### ✅ Step 2. Collect Orthologous Protein Sequences

_Question_: I need to select species other than humans to compare with. How do I decide which species to select? Should I just select randomly?

_Answer_: Do not select randomly - we want to select for biologically meaningful coverage across the evolutionary tree. 

Some guidelines for selection:
- Span vertebrates and invertebrates 
- Vary in nervous system complexity 
- Have good sequence availability

#### Example selection (group level):

    -Mammals
    -Birds
    -Reptiles
    -Amphibians
    -Bony Fish
    -Insects
    -Nematodes
    -Cephalopods
    -Jawless Vertebrates

**Orthologous Genes**
Genes that evolved from a common ancestral gene by speciation and generally perform the same function.

_Question_: With the diversity of groups we've selected (insects, nematodes, etc.), are we sure
that the epileptic genes we've selected will be found in all of them, e.g. will it be a orthologous gene in this group?

_Answer_: No, it is possible these genes won't be present in some of the groups selected. However,
knowing which groups have these gene(s) absent can provide insight for the research too.

When selecting, keep this in mind:

1. **Start with well-studied model organisms:**

   * Well-annotated genomes and proteomes available
   * Widely used in research with abundant data and tools
   * Representative of the taxonomic group and relevant to neuroscience/evolution

2. **Check ortholog availability:**

   * Use NCBI Gene, Ensembl, OrthoDB to confirm gene orthologs exist
   * Avoid species with missing or poorly annotated orthologs

3. **Consider evolutionary distance and biological relevance:**

   * Include both close relatives (fine differences) and distant relatives (conserved regions)
   * Select species with interesting traits linked to your research question

4. **Verify data quality:**

   * Prefer species with high-quality reference genomes and protein data
   * Avoid species with incomplete or low-quality sequence data

#### Final selection (species level):

- Sourced from https://www.ncbi.nlm.nih.gov/
- Choosing proteins for later MSA - proteins are more "obviously" informative than DNA sequences
  for this research project

Mammal: Rabbit (Oryctolagus cuniculus)
Reasoning:
- [x] Species has a publicly available genome and proteome - https://www.ncbi.nlm.nih.gov/datasets/taxonomy/9986/
- [ ] Species is commonly used in research (e.g. neuroscience, evolution)
- [ ] Species represents its group well (e.g. classic model organism)
