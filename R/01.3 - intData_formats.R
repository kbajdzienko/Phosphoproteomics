# Function to extract a classic intensity table with peaks in rows, samples in columns
# Output is a table with first column peak_ID
intTable <- function (df) {
  intTable <-
    df$intData %>%
    #mutate(sample_ID = stringr::str_extract(sample_ID, "(\\d{1,}$)")) %>%
    spread(sample_ID, intensity)
  return(intTable)
}

# Function to extract a classic intensity matrix with peaks in rows, samples in columns
# Output is a matrix with rows called
intMatrix <- function (df) {
  intMatrix <- intTable(df)
  intMatrix <-
    intMatrix %>%
    select(-peak_ID) %>%
    as.matrix() %>%
    'rownames<-'(intMatrix$peak_ID)
  return(intMatrix)
}

# Reverse function to convert product of intTable or intMatrix into tidy format
intData <- function(intTable) {
  if (is.data.frame(intTable)) {
    intTable <-
      tbl_df(intTable) %>%
      gather(sample_ID, intensity, -peak_ID)
  } else if (is.matrix(intTable)) {
    intTable <-
      tbl_df(intTable) %>%
      mutate(peak_ID = rownames(intTable)) %>%
      gather(sample_ID, intensity, -peak_ID)
  } else {
    stop("The intensity table is not a data frame neither a matrix")
  }
}
