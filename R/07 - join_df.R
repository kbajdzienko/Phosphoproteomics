join_df <- function(df1, df2) {
  if(is.null(df1$annData) | is.null(df2$annData))
    stop("Merge sites first to get annData")
  sampleData <- bind_rows(df1$sampleData, df2$sampleData)
  annData <- inner_join(df1$annData, df2$annData, by="ann_ID")
  annIntData <- bind_rows(semi_join(df1$annIntData, annData),
                          semi_join(df2$annIntData, annData))
    
  df <- list(sampleData = sampleData,
             annData = annData,
             annIntData = annIntData)
  return(df)
}

join_df_proteins <- function(df1, df2) {
  
  sampleData <- bind_rows(df1$sampleData, df2$sampleData)
  peakData <- inner_join(df1$peakData, df2$peakData, by="peak_ID")
  intData <- bind_rows(semi_join(df1$intData, peakData),
                       semi_join(df2$intData, peakData))
  
  df <- list(sampleData = sampleData,
             peakData = peakData,
             intData = intData)
  return(df)}