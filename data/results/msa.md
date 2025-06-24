# MSA Observations – [Project_Name]
**Date:** YYYY-MM-DD  
**Species / Sequences Analyzed:**  
- Species A – accession XYZ  
- Species B – accession ABC  
- Species C – accession DEF  

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

