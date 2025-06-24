# MSA Observations – [Comparative-Genomics-Pipeline]
**Date:** 2025-06-24  

# Ortholog Collection Summary  
**Genes Selected:** SCN1A, KCNQ2, SCN2A, CACNA1H, GABRA1  
**Species:** Homo sapiens, Oryctolagus cuniculus, Gallus gallus  
**Notes:** All genes have at least 3 orthologous protein sequences collected. In some cases, UniProt or Entrez IDs are missing; these will be updated as needed.

---

## Gene: SCN1A
| Species               | UniProt ID | Entrez Protein ID     |
|----------------------|------------|------------------------|
| Homo sapiens          | P35498     | *(missing)*            |
| Oryctolagus cuniculus | G1SSP8     | *(missing)*            |
| Gallus gallus         | A0A8V1A1D3 | XP_004942836.4         |

## Gene: KCNQ2
| Species               | UniProt ID | Entrez Protein ID     |
|----------------------|------------|------------------------|
| Homo sapiens          | O43526     | *(missing)*            |
| Oryctolagus cuniculus | *(missing)*| XP_051685089.1         |
| Gallus gallus         | A0A8V0ZLH7 | XP_040506889.1         |

## Gene: SCN2A
| Species               | UniProt ID | Entrez Protein ID     |
|----------------------|------------|------------------------|
| Homo sapiens          | Q99250     | *(missing)*            |
| Oryctolagus cuniculus | A0A5F9C773 | XP_008256914.1         |
| Gallus gallus         | A0A8V1A5Z8 | XP_040559277.1         |

## Gene: CACNA1H
| Species               | UniProt ID | Entrez Protein ID     |
|----------------------|------------|------------------------|
| Homo sapiens          | O95180     | *(missing)*            |
| Oryctolagus cuniculus | *(missing)*| XP_051686345.1         |
| Gallus gallus         | A0A8V0Z849 | XP_015149910.3         |

## Gene: GABRA1
| Species               | UniProt ID | Entrez Protein ID     |
|----------------------|------------|------------------------|
| Homo sapiens          | P14867     | *(missing)*            |
| Oryctolagus cuniculus | G1TAX4     | XP_002710402.2         |
| Gallus gallus         | A0A8V0YC88 | XP_046756226.1         |

---

## 1. Alignment Summary
- **Alignment tool:** Clustal Omega version X.X  
- **Input format:** FASTA MSA  
- **Output format(s) used/viewed:** FASTA / PHYLIP / Clustal  

---

## 2. Visual Inspection Notes
| Region (positions) | Conservation Level | Observations / Comments |
|--------------------|--------------------|--------------------------|
| 1–40               | High               | Completely conserved across all species; may indicate functional core |
| 41–60              | Low-to-moderate    | Gaps in Species B and C; possible insertion in Species A |
| 61–90              | Moderate           | Several substitutions, mostly conservative (polar ↔ polar) |
| 91–110             | Poor               | Many gaps & mismatches—consider trimming |
| 111–end            | High               | Terminal region well-conserved |

---

## 3. Domain / Motif Overlay (if applicable)
- **Domain:** (e.g., transmembrane helix, catalytic site) spans positions X–Y  
- **Motif pattern:** (e.g., GxGxxG) at positions P–Q; conserved in Species A/B only  

---

## 4. Suspicious Regions
- **Region positions:** 41–50  
- **Issue:** Large insertion in Species C; unclear if it’s a real structural feature or alignment artifact  
- **Action needed:** Re-check alignment with local tool or reference genomic sequence  

---

## 5. Preliminary Decisions
- ✅ Retain well-conserved regions (1–40, 111–end) for phylogenetic analysis  
- ⚠️ Mark region (91–110) for **potential trimming**  
- ⏳ Defer trimming until after generating preliminary tree  

---

## 6. Questions / Next Steps
1. Does the insertion in Species C correspond to a known exon/intron boundary?  
2. Should I test trimming with `trimAl` and compare trees?  
3. Which phylogenetic tool (IQ-TREE vs. FastTree vs. MEGA) should I try first?

---

## 7. Archive & References
- **MSA file:** `project_name_alignment.fasta`  
- **Screenshots:**  
  - `screen1_conserved.png` showing positions 1–40  
  - `screen2_indels.png` highlighting gaps  
- **Tools used:**  
  - MSAViewer (web)  
  - Jalview (local)  

---

