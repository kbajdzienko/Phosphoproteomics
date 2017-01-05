
# This funciton is basically just conversion of the file into the set of easy-to-use dataframes

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

  # Change the future variable names in line 3
  df[3, 1:12] <- c("peak_ID", "Ions_used", "Ions", "Decon_ions", "Charge",
                   "RT", "Neutral_mass", "Score", "Sequence", "Modifications",
                   "Accession", "Description")


  if (abundance == "Normalized") {
    col.int <- col.norm:(col.raw-1)
  } else if (abundance == "Raw") {
    col.int <- col.raw:(col.spect-1)
  }

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
    select(1, 6:12) %>%
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
    df$peakData %>%
    filter(Score < score_threshold) %>%
    select(peak_ID)

  df$peakData <- anti_join(df$peakData, shitpeaks)
  df$intData <- anti_join(df$intData, shitpeaks)
  return(df)
}


#Calculate coeficient of variation within replicates of each group
peakCV <- function(df) {
  df_cv <-
    df$intData %>%
    mutate(intensity = zero.to.na(intensity)) %>%
    left_join(df$sampleData) %>%
    group_by(peak_ID, group) %>%
    summarize(CV = sd(intensity, na.rm = T)/mean(intensity, na.rm = T)) %>%
    spread(group, CV)
  return(df_cv)
}

#Sum intensities of duplicate phos site entries in each sample

annCV <- function(df) {
  df_cv_sum2 <-
    df$intData %>%
    left_join(select(df$peakData, peak_ID, ann_ID)) %>%
    group_by(ann_ID, sample_ID) %>%
    summarize_at(vars(intensity), sum, na.rm = T) %>%
    mutate(intensity = zero.to.na(intensity)) %>%
    left_join(df$sampleData) %>%
    group_by(ann_ID, group) %>%
    summarize(CV = sd(intensity, na.rm = T)/mean(intensity, na.rm = T)) %>%
    spread(group, CV)
  return(df_cv_sum2)
}

# Function to plot histograms of several data frames' columns to pdf file
plot_list_hist <- function(df_list, file) {

  # Get the number of columns in dataframes
  n <- unique(sapply(df_list, ncol))
  if (length(n) != 1) stop("Number of columns in data frames is different")

  # Check if all the names are the same between data frames
  name_test <-
    sapply(df_list, function(df) names(df)[-1]) %>%
    apply(1, function(group_ID) length(unique(group_ID)) == 1) %>%
    all()
  if (!name_test) stop("Group names are different in data frames")


  # Open pdf to plot
  pdf(file)
  par(mfrow=c(3,2),
      oma = c(0, 0, 2, 0))

  for (i in 2:n) {
    for (j in 1:length(df_list)) {
      CV <- df_list[[j]][[i]]
      hist(CV,
           xlim = c(0, 2.5),
           breaks = seq(0, 2.5, length.out = 19),
           xlab = "CV",
           main = names(df_list)[j])
      abline(v = median(CV, na.rm = T), col = "red")
    }
    mtext(names(df_list[[1]][i]), outer = TRUE)
    par(mfrow=c(3,2))
  }
  dev.off()
}


# Function to plot CV histograms to pdf file
plot_list_hist_cv <- function(..., file = "CV_report.pdf") {
  df_names <- as.character(substitute(...()))
  df_list <- c(sapply(list(...), function(df) list(peakCV(df), annCV(df))))
  names(df_list) <- c(sapply(df_names,
                             paste0, c("_cv", "_cv_sum")))
  plot_list_hist(df_list, file)
}

# Sum intensities for the same protein, the same modification
sumANN <- function(df) {
  df$intData <-
    df$intData %>%
    left_join(select(df$peakData, peak_ID, ann_ID)) %>%
    group_by(ann_ID, sample_ID) %>%
    summarize_at(vars(intensity), sum, na.rm = T) %>%
    mutate(intensity = zero.to.na(intensity)) %>%
    ungroup()
  return(df)
}



