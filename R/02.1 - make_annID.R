# Create ann_ID in format "AT3G47070.1_T65_T66"
# Accession, phosphorilated amino acid, number

make_annID <- function(df, mascot_file) {

  # Read mascot and rename some variables for future joining
  # Remove everything but Sequence, Accession, start and end position
  # and number of missed cleaveges
  mascot <-
    read.mascot(mascot_file) %>%
    rename(Accession = prot_acc,
           Sequence = pep_seq,
           Score = pep_score,
           Neutral_mass = pep_exp_mr) %>%
    select(Accession, Sequence, pep_miss, pep_start, pep_end) %>%
    distinct()

  # Peaks not found in Mascot for future curation
  miss_peakData <-
    df$peakData %>%
    anti_join(mascot, by = c("Sequence", "Accession"))

  # Merge original annotation data table with mascot
  # Missed peaks are left out (use inner_join)
  peakData <-
    df$peakData %>%
    inner_join(mascot)

  # Curation of NAs by considering I/L as the same
  peakData <- bind_rows(peakData, cure_IL_NA(miss_peakData, mascot))

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
cure_IL_NA <- function (miss_peakData, mascot) {

  # Get sequences/accessions with no match
  IL_seq <- distinct(miss_peakData, Sequence, Accession)

  # Instead of sequence form regular expression with optional I/L positions
  IL_seq$ILSequence <-
    strsplit(IL_seq$Sequence, "") %>%
    lapply(gsub, pattern = "[LI]", replacement = "[LI]") %>%
    sapply(function(x) paste0(c("^", x, "$"), collapse = ""))

  # Create future data frame for joining IL_seq_na
  IL_mascot <- NULL

  for (i in 1:nrow(IL_seq)) {
    acc. <- IL_seq$Accession[i]
    seq. <- IL_seq$Sequence[i]
    ilseq. <- IL_seq$ILSequence[i]

    # Get mascot records with Sequence matching that regular expression and
    # and corresponding Accession

    add_data <-
      mascot %>%
      filter((Accession == acc.) & grepl(ilseq., Sequence)) %>%
      distinct(Sequence, pep_start, pep_end, pep_miss)

    # If not a single sequence, replace all numbers it with NA
    if (nrow(add_data) != 1) {
      if (nrow(add_data) > 1)
        warning("Multiple I/L matches found in Mascot for ", seq., " ", acc.)
      if (nrow(add_data) == 0)
        warning("No I/L matches found in Mascot for ", seq., " ", acc.)

      add_data <-
        data.frame(Sequence = seq.,
                   pep_start = NA, pep_end = NA, pep_miss = NA) %>%
        mutate_at(vars(-Sequence), as.integer) %>%
        tbl_df()
    }

    # In any case add either NAs or a single record to the table for future joining
    IL_mascot <- bind_rows(IL_mascot, add_data)
  }

  # Join
  IL_mascot <-
    IL_mascot %>%
    rename(mascotSequence = Sequence) %>%
    bind_cols(IL_seq, .)

  # Join with miss_peakData and take I/L sequences from mascot data frame
  filled_peakData <-
    left_join(miss_peakData, IL_mascot) %>%
    mutate(Sequence = mascotSequence) %>%
    select(-mascotSequence, -ILSequence)

  return(filled_peakData)
}
