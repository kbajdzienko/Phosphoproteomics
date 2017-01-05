# Remove peptides with more than one Accession based on ion PQI output table

filter_pep_unique <- function (df, ions_file) {
  ions <- read.csv(ions_file,
                   skip = 2,
                   header = T, sep = ",",
                   stringsAsFactors = F,
                   colClasses = "character")

  ions <-
    tbl_df(ions) %>%
    select(1:2, 5, 8:12) %>%
    setNames(c("ion_ID", "RT", "Neutral_mass",
               "Score", "Sequence", "Modifications",
               "Accession", "All_accessions")) %>%
    mutate_at(vars(RT:Score), as.numeric)

  dupl_ions <- filter(ions, grepl(";", ions$All_accessions))

  # Test if the list of accessions is always the same for a certain sequence
  is.diff <- function(all_acc) {
    is.diff <-
      all_acc %>%
      strsplit(";") %>%
      outer(., ., Vectorize(setequal)) %>%
      all() %>%
      !.
    return(is.diff)
  }

  #Apply is.diff function to every sequence and filter the different ones
  dupl_ions_diff <-
    dupl_ions %>%
    group_by(Sequence) %>%
    filter(is.diff(All_accessions))
  if (nrow(dupl_ions_diff) != 0) {
    warning("Some sequences have different lists of possible accessions:\n",
            paste0(unique(dupl_ions$Sequence), "\n"))
  }

  # Exclude non-unique from data
  df_unique <- df
  dupl_ions <- semi_join(df$peakData, by = "Sequence")
  df_unique$peakData <- anti_join(df$peakData, dupl_ions, by = "Sequence")
  df_unique$intData <- anti_join(df$intData, dupl_ions, by = "Sequence")

  # Get the list of proteins not having a unique peptide
  prot.nonuniq <-
    setdiff(unique(df$peakData$Accession),
            unique(df_unique$peakData$Accession))

  if (length(prot.nonuniq) != 0)
    message("Proteins having no unique peptides:\n",
            paste0(prot.nonuniq, "\n"))

  return(df_unique)

  }
