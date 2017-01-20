
dtt_import <- function(df, peptide_file, mascot_file, ions_file,
                      abundance = "Raw") {
  
  df <- tidy_PQI(file = peptide_file, abundance = abundance)
  
  df <- make_annID(df, mascot_file = mascot_file)
  
  df <- make_conf(df, mascot_file = mascot_file)
  
  df <- make_all_acc(df, ions_file = ions_file)
  
  
  
 return(df)
  
}





dtt_clean <- function(
                df,  
                NAimputation = "ppca",
                score = 29
                ) {
  
  df <- filter_phos_conf(df, score_threshold = 70)
  
  df <- filter_score(df, score_threshold = score)
  
  df <- filter_NA(df, threshold = 0.7)  
  
  df <- normMedian(df)
  
  df <- fillNA(df, method = NAimputation)

  df <- sitesMerge(df)
  
  return(df)
  
}