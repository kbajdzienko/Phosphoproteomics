# Function to merge all the dataframes into a single table
# and write the output into the file
report_PQI <- function (df, file) {

  output <-
    df$peakData %>%
    #left_join(df$seqData) %>%
    #left_join(df$annData) %>%
    left_join(intTable(df))

  write.csv(output, file, row.names = F)

  return(output)
}
