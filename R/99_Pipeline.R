#do the thing function
#Function to automate data processing pipeline
# df <- dtt(mascot_file = "", 
#           ions_file = "",
#           peptide_file = "", 
#           min_score = 10)

dtt <- function(df, peptide_file, mascot_file, ions_file, 
                abundance = "Raw", NA_allowed = 0.5, 
                min_score = 10, NAimputation = "ppca", base = 2) {
  
  df <- tidy_PQI(file = peptide_file, abundance = abundance)
  
  # Creata ann_ID in format "AT3G47070.1_T65_T66"
  # Accession, phosphorylated amino acid, number
  df <- make_annID(df, mascot_file = mascot_file)
  
  # Remove peptides with more than one Accession based on ion PQI output table
  #WARNING: This function can return a warning: NAs introduced by coercion. This is because peptides
  #after charge deconvolution in QI can be left with no score
  df <- filter_pep_unique(df, ions_file = ions_file)
  
  # Normalization by median intensity of all the peaks in sample
  df <- normMedian(df)
  
  # Filter peaks with too many NA.
  # Keep a peak if in at least one group there is less NA, then the set threshold
  df <- filter_NA(df, threshold = NA_allowed)  
  
  # Filter peaks with Score less than the set threshold
  df <- filter_score(df, score_threshold = min_score)
  
  # Fill missing values with PCA-based methods (PPCA or BPCA)
  df <- fillNA(df, method = NAimputation)
  
  # Log-transform data
  df <- logTransform(df, base = base) 
  
  # Split not cleaved peptides into single phosphorilation sites and sum intensitites for the same ann_ID (protein/site)
  # Save intensity data in annIntData table
  # Save annotation data in annData table
  df <- sitesMerge(df)
  
  return(df)
  
}