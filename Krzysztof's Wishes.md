#### Wishes:

- [x] script to extract query table or peptide table from mascot file
- [x] calculate confidence score (similar to delta-score)
- [x] 99 report quality script -> add functions to report quality metrics and plot histograms from tidy_df's
- [x] distribution (frequency, not density) of single, double, triple and quadruple phosphorylated peptides (taken from df)
- [ ] plots: clustering of dynamic phosphorylation profiles, Functional Categories of Regulated Phosphoproteins (please check F1000 :Fig3 and Fig4, Global, in vivo, and site-specific phosphorylation dynamics in signaling networks.)
- [x] histograms -> rewrite from counts to density
- [x] Profile plot -- all peaks by Accession.
- [ ] PCA (with loadings)
- [ ] Volcano plot (p-val/fold change)

- [ ] for future: import profile for MaxQuant output to PDA workflow

#### FIX:
- [ ] sitesMerge - different number of ann_IDs in peakData and annIntData - HOW
- [ ] ADD check for uniqueness of sequence-accession-mass-score combination for merging df and mascot in make_annID
- [ ] Round Neutral mass in make_annID before joining df and mascot

#### Suggestions:
- [ ] Missing values statistics

#### Disposed Wishes:
- [ ] From tidy df -> identified unique proteins within each run -- median and CV
- [ ] calculate mascot delta score and normalized mascot delta score