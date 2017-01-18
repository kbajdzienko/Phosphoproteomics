
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
                NAimputation = "ppca",
                score = 29
                ) {
  
  filter_phos_conf(df, score_threshold = 70)
  
  filter_score(df, score_threshold = score)
  
  df <- filter_NA(df, threshold = 0.7)  
  
  df <- normMedian(df)
  
  df <- fillNA(df, method = NAimputation)

  df <- sitesMerge(df)
  
  return(df)
  
}