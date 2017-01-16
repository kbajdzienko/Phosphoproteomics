#do the thing function
#Function to automate data processing pipeline
# df <- dtt(mascot_file = "", 
#           ions_file = "",
#           peptide_file = "", 
#           min_score = 10)



dtt_pre <- function(df, peptide_file, mascot_file, 
                    ions_file, 
                abundance = "Raw", NA_allowed = 0.5, base = 2) {
  
  df <- tidy_PQI(file = peptide_file, abundance = abundance)
  
  df <- make_annID(df, mascot_file = mascot_file)
  
  df <- filter_pep_unique(df, ions_file = ions_file)
  
  df <- filter_NA(df, threshold = NA_allowed)  
  
  df <- logTransform(df, base = base) 
  
  return(df)
  
}





dtt <- function(df, peptide_file, mascot_file, 
                ions_file, 
                abundance = "Raw", NA_allowed = 0.5, 
                NAimputation = "ppca", base = 2) {
  
  df <- tidy_PQI(file = peptide_file, abundance = abundance)
  
  df <- make_annID(df, mascot_file = mascot_file)
  
  df <- filter_pep_unique(df, ions_file = ions_file)

  df <- normMedian(df)
  
  df <- filter_NA(df, threshold = NA_allowed)  
  
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