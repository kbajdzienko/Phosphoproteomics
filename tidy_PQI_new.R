for (pckg in c("plyr", "dplyr", "tidyr"))  {
  if (!require(pckg, character.only = T)) install.packages(pckg)
}

library(plyr)
library(dplyr)
library(tidyr)

source("report_PQI_new.R")
source("make_annID_new.R")

# This funciton is basically just conversion of the file into the set of easy-to-use dataframes

tidy_PQI <- function(file) {
  
  df <- read.csv(file,
                 header = F, sep = ",",
                 stringsAsFactors = F,
                 colClasses = "character")
  
  df <- tbl_df(df)
  
  # Get indexes of columns where intensity table starts
  col.norm <- grep("Normalized abund", df[1, ])
  col.raw <- grep("Raw abund", df[1, ])
  col.spect <- grep("Spectral coun", df[1, ])
  
  # Change the future variable names in line 3
  df[3, 1:12] <- c("peak_ID", "Ions_used", "Ions", "Decon_ions", "Charge",
                   "RT", "Neutral_mass", "Score", "Sequence", "Modifications",
                   "Accession", "Description")
  
  # Data about samples
  sampleData <- 
    df %>%
    slice(2:3) %>%
    select(col.raw:(col.spect-1)) %>%
    t() %>%
    tbl_df() %>%
    setNames(c("group", "sample_ID"))
  
  # Expand groupnames
  for (i in 2:nrow(sampleData)) {
    if (sampleData[i, "group"] == "") 
      sampleData[i, "group"] <- sampleData[i - 1, "group"]
  }
  
  #Rename samples
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
    select(1, col.raw:(col.spect-1)) %>%
    setNames(c("peak_ID", sampleData$sample_ID)) %>%
    slice(-(1:3)) %>%
    distinct() %>%
    mutate(peak_ID = rename_dupl(.$peak_ID)) %>%
    mutate_each(funs(as.numeric), -peak_ID) %>%
    gather("sample_ID", "intensity", -peak_ID) %>%
    arrange(peak_ID)
  
  # Peak information table
  peakData <- 
    df %>%
    select(1, 6:12) %>%
    setNames(.[3, ]) %>%
    slice(-(1:3)) %>%
    distinct() %>%
    mutate(peak_ID = rename_dupl(.$peak_ID)) %>%
    mutate_each(funs(as.numeric), RT, Neutral_mass) %>%
    mutate_each(funs(as.numeric), Score) %>%
    arrange(RT)
  
  return(list(sampleData = sampleData,
              peakData = peakData,
              intData = intData))
}

filter_NA <- function(df, threshold = 0.5) {
  
  shitpeaks <- 
    df$intData %>%
    left_join(df$sampleData) %>%
    mutate(intensity = zero.to.na(intensity)) %>%
    group_by(peak_ID, group) %>%
    summarize(NA_ratio = sum(is.na(intensity))/n()) %>%
    group_by(peak_ID) %>%
    filter(!any(NA_ratio < threshold)) %>%
    select(peak_ID) %>%
    distinct()
  
  df$peakData <- anti_join(df$peakData, shitpeaks)
  df$intData <- anti_join(df$intData, shitpeaks)
  return(df)
}

zero.to.na <- function(x) {
  x[x==0] <- NA
  return(x)
}


# Use like that:
#df <- filter_NA(df)
# OR
#df <- filter_NA(df, threshold = 0.2)


filter_score <- function(df, score_threshold = 5) {
  
  shitpeaks <- 
    df$intData %>%
    left_join(df$peakData) %>%
    group_by(peak_ID) %>%
    filter(!any(Score > score_threshold)) %>%
    select(peak_ID) %>%
    distinct()
  
  df$peakData <- anti_join(df$peakData, shitpeaks)
  df$intData <- anti_join(df$intData, shitpeaks)
  return(df)
}