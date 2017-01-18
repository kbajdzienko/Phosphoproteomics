#### Wishes:

- [x] script to extract query table or peptide table from mascot file
- [x] calculate confidence score (similar to delta-score)
- [x] 99 report quality script -> add functions to report quality metrics and plot histograms from tidy_df's
- [x] distribution (frequency, not density) of single, double, triple and quadruple phosphorylated peptides (taken from df)
- [x] histograms -> rewrite from counts to density
- [x] Profile plot -- all peaks by Accession.
- [x] Pareto and Auto scaling + intensity profile (like in Metaboanalyst)
- [x] Reshape the quality check table

- [ ] plots: fuzzy clustering of dynamic phosphorylation profiles
- [ ] plot: PCA (with loadings)
- [ ] plot: Volcano plot (p-val/fold change)
- [x] plot: ann_ID intensity profile, choice of which ann_ID are plotted
- [x] ions table -> join All_Accessions to peakData
- [ ] for future: import function for MaxQuant output

#### FIX:
- [x] sitesMerge - different number of ann_IDs in peakData and annIntData - (was no problem)
- [x] ADD check for uniqueness of sequence-accession-mass-score combination for merging df and mascot in make_annID
- [x] Round Neutral mass in make_annID before joining df and mascot

#### Suggestions:
- [ ] Missing values statistics
- [x] Improve read.mascot performance for big files
- [x] Improve make_conf performance
- [ ] Simplify filter_unique_peptides since Accessions in All_accessions are arranged in the same way  
- [ ] Update comments
