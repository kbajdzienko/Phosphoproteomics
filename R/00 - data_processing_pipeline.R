
dtt_import <- function(df, peptide_file, mascot_file, 
                      abundance = "Raw") {
  
  df <- tidy_PQI(file = peptide_file, abundance = abundance)
  
  df <- make_annID(df, mascot_file = mascot_file)
  
  df <- make_conf(df, mascot_file)
  
  
  
 return(df)
  
}





dtt_clean <- function(
                df,  
                ions_file, 
                NA_allowed = 0.5, 
                NAimputation = "ppca", 
                base = 2) {
  
  df <- filter_pep_unique(df, ions_file = ions_file)
  
  df <- filter_NA(df, threshold = NA_allowed)  
  
  df <- normMedian(df)
  
  df <- fillNA(df, method = NAimputation)
  
  df <- logTransform(df, base = base) 
  
  df <- normScale(df, method = "auto")
  
  df <- sitesMerge(df)
  
  return(df)
  
}