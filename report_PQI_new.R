library(plyr)
library(dplyr)
library(tidyr)

# Function to get a classic intensity table with peaks in rows, samples in columns 
intTable <- function (df) {
  intTable <- 
    df$intData %>%
    spread(sample_ID, intensity)
  return(intTable)
}

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