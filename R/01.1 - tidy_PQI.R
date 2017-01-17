# Conversion of the file into the list of easy-to-use dataframes
# All zeros are replaced with NA

tidy_PQI <- function(file, abundance = "Raw") {

  match.arg(abundance, c("Normalized", "Raw"))
  df <- read.csv(file,
                 header = F, sep = ",",
                 stringsAsFactors = F,
                 colClasses = "character")

  df <- tbl_df(df)

  # Get indexes of columns where intensity table starts
  col.norm <- grep("Normalized abund", df[1, ])
  col.raw <- grep("Raw abund", df[1, ])
  col.spect <- grep("Spectral coun", df[1, ])

  # Select either normalized or raw columns
  if (abundance == "Normalized") {
    col.int <- col.norm:(col.raw-1)
  } else if (abundance == "Raw") {
    col.int <- col.raw:(col.spect-1)
  }

  # Change the future variable names in line 3
  df[3, 1:12] <- c("peak_ID", "Ions_used", "Ions", "Decon_ions", "Charge",
                   "RT", "Neutral_mass", "Score", "Sequence", "Modifications",
                   "Accession", "Description")

  # Data about samples
  sampleData <-
    df %>%
    slice(2:3) %>%
    select(col.int) %>%
    t() %>%
    tbl_df() %>%
    setNames(c("group", "sample_ID"))

  # Expand groupnames
  for (i in 2:nrow(sampleData)) {
    if (sampleData[i, "group"] == "")
      sampleData[i, "group"] <- sampleData[i - 1, "group"]
  }

  #Rename samples: add group name
  sampleData <-
    sampleData %>%
    mutate(sample_ID = paste(substr(group, 1, 3),
                             gsub("\\D", "", group),
                             sample_ID,
                             sep = "_"))


  # Funciton, adding _1, _2, _3, ... indexes to duplicated peak IDs
  rename_dupl <- function(ids) {
    dupl <- duplicated(ids) | duplicated(ids, fromLast = T)
    id.dupl <- ids[dupl]
    for (i in 1:length(id.dupl)) {
      j <- id.dupl == id.dupl[i]
      if (sum(j) > 1) id.dupl[j] <- paste(id.dupl[j], 1:sum(j), sep = "_")
    }
    ids[dupl] <- id.dupl
    return(ids)
  }

  # Intensity table
  intData <-
    df %>%
    select(1, col.int) %>%
    setNames(c("peak_ID", sampleData$sample_ID)) %>%
    slice(-(1:3)) %>%
    distinct() %>%
    mutate(peak_ID = rename_dupl(.$peak_ID)) %>%
    gather("sample_ID", "intensity", -peak_ID) %>%
    mutate_at(vars(intensity), as.numeric) %>%
    mutate_at(vars(intensity), zero.to.na) %>%
    arrange(peak_ID)

  # Peak information table
  peakData <-
    df %>%
    select(1, 4, 6:12) %>%
    setNames(.[3, ]) %>%
    slice(-(1:3)) %>%
    distinct() %>%
    mutate(peak_ID = rename_dupl(.$peak_ID)) %>%
    mutate_at(vars(RT, Neutral_mass, Score), as.numeric) %>%
    arrange(peak_ID)

  return(list(sampleData = sampleData,
              peakData = peakData,
              intData = intData))
}
