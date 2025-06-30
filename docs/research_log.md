# Research Log

## 2025-06-28
Confirmed that loss of function variants occurred in highly conserved regions, which was expected and validates the accuracy of the data pipeline.

**Question:** After asking what "loss of function" means and learning that SCN1A LoF variants reduce sodium current, I expected this would cause hypoexcitability, but Dravet syndrome causes hyperexcitability and seizures. Why does REDUCING sodium current promote hyperexcitability?

**Answer:** SCN1A is preferentially expressed in GABAergic interneurons. LoF variants impair interneuron firing, reducing GABA release and disinhibiting the network. Excitatory neurons use multiple sodium channel subtypes and are less SCN1A-dependent.

**Reference:** [Science Translational Medicine - Interneuron-specific dual-AAV SCN1A gene replacement corrects epileptic phenotypes in mouse models of Dravet syndrome](https://www.science.org/doi/10.1126/scitranslmed.adn5603)

Updated conservation-variant plot to separate LoF variants (red) from all variants (orange) in the bottom histogram for better visualization of where critical variants cluster in the conservation distribution.

**Question:** Are variants more likely to occur in regions that are highly conserved because of nature or because we more strictly analyze variants in sequences that are highly conserved and it's a part of their definition?

**Answer:** This is primarily a methodological artifact rather than biological reality. The apparent enrichment of variants in conserved regions occurs because: (1) We preferentially analyze protein-coding and functionally important regions that are inherently more conserved, (2) These regions receive more clinical scrutiny leading to better variant detection, and (3) Biologically, variants in truly conserved regions are actually less frequent due to purifying selection, but when present are more likely to be pathogenic and thus reported. The SCN1A analysis validates this - LoF variants cluster in conserved functional domains where they have measurable phenotypic impact.

## 2025-06-29
**Plot Simplifications:** Removed density distribution panel from conservation plots for cleaner scientific presentation. Single-panel design maintains all essential conservation statistics in legend while improving readability.

**Research Focus Update:** Pipeline now configured for dual-gene analysis focusing exclusively on SCN1A (Dravet syndrome, voltage-gated sodium channel) and DEPDC5 (focal epilepsy, mTOR pathway regulator). Removed references to other epilepsy genes to concentrate analysis on these two distinct but complementary mechanisms: SCN1A affecting interneuron excitability through sodium channel dysfunction, and DEPDC5 affecting cortical development through mTOR pathway dysregulation.

**DEPDC5 data added:**
Early data from DEPDC5 looks promising but raises questions on how best to display comparisons with SCN1A, since no LoF variants have been identified in DEPDC5 and the focus will be on likely pathogenic variants instead. The plan is to map loss-of-function, likely pathogenic, and pathogenic variants in both SCN1A and DEPDC5 to capture key candidates. Raw variant data for both genes will be re-verified and vetted to ensure accuracy.

## 2025-06-30
Starting to dig into the ClinVar data to see how it compares with the variant info I got from UniProt. Want to cross-reference these two sources since they might have different coverage or different ways of classifying the same variants. ClinVar tends to be more clinically focused while UniProt has broader variant annotations, so comparing them should give a better picture of what variants are actually relevant for patients vs. what's just been observed in the lab.