#### Wishes:

- [x] script to extract query table or peptide table from mascot file
- [x] calculate confidence score (similar to delta-score)
- [x] 99 report quality script -> add functions to report quality metrics and plot histograms from tidy_df's
- [x] distribution (frequency, not density) of single, double, triple and quadruple phosphorylated peptides (taken from df)
- [x] histograms -> rewrite from counts to density
- [x] Profile plot -- all peaks by Accession.
- [x] Pareto and Auto scaling + intensity profile (like in Metaboanalyst)
- [x] Reshape the quality check table
- [x] plot: PCA (with loadings)
- [x] plot: Volcano plot (p-val/fold change)
- [x] plots: fuzzy clustering of dynamic phosphorylation profiles


- [ ] 2wayANOVA test (randomized blocks beacause all replicates are independent?)
- [ ] plot: PLS-DA
- [ ] fuzzy clustering - export list of accessions for each cluster (export, csv, col: siteID, cluster)
- [ ] for applicable functions: add parameter to choose: intData, andIntData
- [ ] change "ann" to "psite"
- [ ] Add whiskers for time profile plot as an option

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
- [x] Ability to make DF with subset of original data (list of patterns in sample names?, e.g. - AZD/GLU-015; GLU-015)
- [ ] for future: import function for MaxQuant output
- [ ] fucking Gene Ontology