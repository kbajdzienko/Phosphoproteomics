# Create ann_ID in format "AT3G47070.1_T65_T66"
# Accession, phosphorilated amino acid, number

make_annID <- function(df, mascot_file) {

  # Read mascot and rename some variables for future joining
  mascot <-
    read.mascot(mascot_file) %>%
    rename(Accession = prot_acc,
           Sequence = pep_seq,
           Score = pep_score,
           Neutral_mass = pep_exp_mr) %>%
    distinct()

  # Function removing NAs by considering I/L as the same
  mascot <- cure_IL_NA(df, mascot)

  # Remove everything but Sequence, Accession, start and end position
  # and number of missed cleaveges and
  # merge original annotation data table with mascot
  peakData <-
    df$peakData %>%
    left_join(select(mascot, Accession, Sequence, pep_miss, pep_start, pep_end))

  # Separate modification columns into single modifications
  # Select only phospho
  # Remove everything but phosphate position number in sequence
  phosph_pos <-
    peakData$Modifications %>%
    strsplit("\\|") %>%
    lapply(function(x) as.integer(gsub("\\D", "", x[grepl("Phospho", x)])))

  # Function looking up aminoacid by position in peptide and merging it
  # with position number in protein
  aaposMerge <- function(sequence, pep_start, phosph_pos) {
    paste0(unlist(strsplit(sequence, ""))[phosph_pos],
           pep_start + phosph_pos,
           collapse = "_")
  }

  # Merge accession and aminoacid+position
  ann_ID <- paste(peakData$Accession,
                  mapply(aaposMerge, peakData$Sequence, peakData$pep_start, phosph_pos),
                  sep = "_")

  # Add new ann_ID to peakData table
  df$peakData <-
    peakData %>%
    mutate(ann_ID = ann_ID) %>%
    select(peak_ID:Score, ann_ID, everything())

  return(df)
}


# Function, adding to mascot missing Sequence/Accessions in case of I/L bug
cure_IL_NA <- function (df, mascot){
  # Funciton fixing NAs appearance
  IL_seq_na <-
    df$peakData %>%
    anti_join(mascot, by = c("Sequence", "Accession")) %>%
    distinct(Sequence, Accession)

  # Form regular expression instead of sequence with optional I/L positions
  IL_seq_na$ILSequence <-
    IL_seq_na$Sequence %>%
    strsplit("") %>%
    lapply(gsub, pattern = "[LI]", replacement = "[LI]") %>%
    lapply(function(x) paste0(c("^", x, "$"), collapse = "")) %>%
    unlist()

  pep_starts <- list()
  for (i in 1:nrow(IL_seq_na)) {
    add_data <-
      mascot %>%
      filter(Accession == IL_seq_na$Accession[i] & grepl(IL_seq_na$ILSequence[i], Sequence)) %>%
      distinct(pep_start, pep_end, pep_miss)
    if (nrow(add_data) > 1) {
      warning("Multiple I/L matches found in Mascot for ",
              IL_seq_na$Sequence[i], " ", IL_seq_na$Accession[i])
      add_data <-
        data.frame(pep_start = NA, pep_end = NA, pep_miss = NA) %>%
        mutate_all(as.integer) %>%
        tbl_df()
    } else if (nrow(add_data) == 0) {
      warning("No I/L matches found in Mascot for ",
              IL_seq_na$Sequence[i], " ", IL_seq_na$Accession[i])
      add_data <-
        data.frame(pep_start = NA, pep_end = NA, pep_miss = NA) %>%
        mutate_all(as.integer) %>%
        tbl_df()
    }
    pep_starts[[i]] <- add_data
  }

  mascot <-
    IL_seq_na %>%
    bind_cols(Reduce(bind_rows, pep_starts)) %>%
    select(-ILSequence) %>%
    bind_rows(mascot, .)
  return(mascot)
}
