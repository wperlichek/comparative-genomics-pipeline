# Phylogenetic Tree Analysis Worksheet

**Note:** These are *preliminary* results intended to build confidence in the pipeline and are not final conclusions.

### Phylogenetic Tree Analysis

1. **What major clades or groupings are observed in the tree?**
   Grouping: Human and rabbit are grouped together. Chicken usually its own group. This makes sense as
   human/rabbit are mammals and chickens are birds. 

2. **Are there any unexpected relationships or placements in the tree?**
   No. We see GABRA1 SCN2A and SCN1A have rabbit/human close and chicken far horizontally, which is intuitive. 
   CACNA1H and KCNQ2 the chicken is closer to human/rabbit horizontally. 
   On a high level, we would expect some proteins to diverge and some to be conserved depending on evolutionary constraints.

3. **How well supported are the branches (if support values are available)?**
   N/A Clustal Omega does not have branch support. Need to use IQ-Tree if we wanted this. 

4. **Do the tree’s groupings match known evolutionary relationships?**
   Yes, mostly. The chicken is a bird that diverged from mammals so we
   expect to see distance from rabbit and human.

5. **What might explain any discrepancies or surprises in the tree?**
   The KCNQ2.png suggests this protein is highly conserved across human,
   rabbit, and chicken. Other genes show more obvious evolutionary
   distance.

6. **How do the expanded groups (reptiles, amphibians, bony fish, insects, nematodes, cephalopods, jawless vertebrates) affect the tree’s structure and interpretation?**
   - Are new major clades or unexpected groupings observed with the inclusion of these additional species?

   Taking SCNA1 as example: all new groups are still about same group length, suggesting as before that SCNA1 is highly conserved and important. 

   Taking GARBA1 as example, it's interesting that human, rabbit, and chicken now appear close together compared to distance in the other groups.

   My key observation here is that although chicken and human might have appear divergent when we only had Human, rabbit, chicken plotted, when we start plotting much different groups, we see how similar human, chicken, and rabbit are in all protein structures we're analyzing.

   - Do the relationships among vertebrates and invertebrates match established evolutionary history?

   

   - Are there any surprising placements or long branches that suggest rapid evolution or gene loss in certain lineages?

7. **What changes, if any, are observed in the findings or interpretations compared to the tree with only mammals and birds?**
   - Are the mammal/bird groupings still as distinct, or do new patterns emerge?

   - Do any of the new groups cluster closer to mammals, birds, or form their own distinct branches?

   - How does the inclusion of distant groups (e.g., insects, nematodes) highlight conserved versus divergent regions of the protein?
