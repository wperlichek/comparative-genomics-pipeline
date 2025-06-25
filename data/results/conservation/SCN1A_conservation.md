# Conservation Score Analysis Worksheet

**Note:** These are _preliminary_ results to help interpret conservation across aligned sequences.

---

## Basic Info

1. **What gene(s) is this analysis for?**
   SCN1A

2. **Which species are included in the alignment?**
   Home sapiens, Gallus gallus, and Oryctolagus cuniculus

3. **How were conservation scores calculated?**
   Shannon entropy per alignment column, consensus residue - see `compute_conservation_scores()` in biopython_service.py

---

## Conservation Results

4. **Which regions are most highly conserved?**
   *~81-210
   *~1254-1383

   There are others but these stand out

5. **Which regions are most variable?**
   *~1900-End
   *~1384-1414
   *~1026-1085
   *~276-327

6. **Are there any notable patterns in conservation (e.g., domains, motifs, or functional sites)?**
   TBD

7. **Do conserved regions correspond to known functional or structural domains?**
   TBD

8. **Are there any regions with missing data or gaps that affect conservation scores?**
   *680-690

---

## Interpretation

9. **What might explain the observed conservation/variability?**
   Some regions are used for core neural functionality and required. 
   Some variable regions don't contain required functionality
   or indicate loss of functionality (like in Gallus gallus 680-690)

10. **How might these results inform further analysis or experiments?**
    Focus on regions which are known to be assoicated with seizures. 
    Then see if they are conservered.
    This will require understanding epilepsy in humans, rabbits, and chickens.

---

## Visualization

- [Optional: Insert or link to plots, heatmaps, or color-coded alignments showing conservation scores.]
