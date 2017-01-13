#Load all needed functions and librarys
source(00 - sources.utility.R)

# Conversion of the file into the list of easy-to-use dataframes
# All zeros are replaced with NA
df <- tidy_PQI("QI peptide output", abundance = "Raw")

# Creata ann_ID in format "AT3G47070.1_T65_T66"
# Accession, phosphorylated amino acid, number
df <- make_annID(df, mascot_file = "mascot peptide and query output")

# Remove peptides with more than one Accession based on ion PQI output table
df <- filter_pep_unique(df, ions_file = "QI peptide ions output")

# Normalization by median intensity of all the peaks in sample
df <- normMedian(df)

# Filter peaks with too many NA.
# Keep a peak if in at least one group there is less NA, then the set threshold
df <- filter_NA(df, threshold = 0.5)  

# Filter peaks with Score less than the set threshold
df <- filter_score(df, score_threshold = 5)

# Fill missing values with PCA-based methods (PPCA or BPCA)
df <- fillNA(df, method = c("ppca", "bpca"))
    
# Log-transform data
df <- logTransform(df, base = 2) 

# Split not cleaved peptides into single phosphorilation sites and sum intensitites for the same ann_ID (protein/site)
# Save intensity data in annIntData table
# Save annotation data in annData table
df <- sitesMerge(df)


######PLOTING######
#Report abundance distibution profile
int_profile_median(df)
# Plot histograms of a distribution of coeficient of variation within replicate group
plot_list_hist_cv(file = "CV_df.pdf", df)
#Quality report from mascot output
QC_stat("mascot_file")
#Quality control report from filtered dataframe (df)
QC_stat(df)