
dtt_import <- function(peptide_file, mascot_file, ions_file,
                      abundance = "Raw") {
  
  df <- tidy_PQI(peptide_file, abundance)
  
  df <- make_annID(df, mascot_file)
  
  df <- make_conf(df, mascot_file)
  
  df <- make_all_acc(df, ions_file)
  
  
  
  
  
 return(df)
  
}


dtt_clean <- function(
                df,  
                NAimputation = "ppca",
                score = 25
                ) {
  
  
  df <- filter_phos_conf(df, score_threshold = 50)
  
  df <- filter_score(df, score=25)
  
  df <- filter_NA_Mann(df)  
  
  df <- normMedian(df)
  
  df <- fillNA(df, "ppca")

  df <- sitesMerge(df)
  
  return(df)
  
}
