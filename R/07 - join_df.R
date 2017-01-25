join_df <- function(df1, df2) {
  if(is.null(df1$annData) | is.null(df2$annData))
    stop("Merge sites first to get annData")
  sampleData <- bind_rows(df1$sampleData, df2$sampleData)
  annData <- inner_join(df1$annData, df2$annData)
  annIntData <- bind_rows(semi_join(df1$annIntData, annData),
                          semi_join(df2$annIntData, annData))
  df <- list(sampleData = sampleData,
             annData = annData,
             annIntData = annIntData)
  return(df)
}
