# Creata ann_ID in format "AT3G47070.1_T65_T66"
# Accession, phosphorilated amino acid, number

make_annID <- function(df, mascot_file) {

  # Reads mascot file, takes peptide table and joins pep_var_mod_conf from query table
  mascot_conf <- function(mascot_file) {
    dupl <- function(x) duplicated(x) | duplicated(x, fromLast = T)

    query <- read.mascot(mascot_file, "query")
    mascot <- read.mascot(mascot_file, "pep")
    query$pep_var_mod_conf[!dupl(query$query_number)] <- 100
    mascot <- left_join(mascot, query)
    return(mascot)
  }

  # Remove everything but Sequence, Accession, start and end position
  # and number of missed cleaveges from mascot output file
  mascot <-
    mascot_conf(mascot_file) %>%
    select(prot_acc, pep_seq, pep_score,
           pep_start, pep_end, pep_miss,
           pep_exp_mr, pep_var_mod_conf) %>%
    dplyr::rename(Accession = prot_acc,
           Sequence = pep_seq,
           Score = pep_score,
           Neutral_mass = pep_exp_mr) %>%
    distinct()

  # Round Neutral_mass to get match between df and mascot
  mascot <- mutate(mascot, Neutral_mass = round(Neutral_mass, 4))
  df$peakData <- mutate(df$peakData, Neutral_mass = round(Neutral_mass, 4))

  # Merge original annotation data table with mascot
  peakData <-
    df$peakData %>%
    left_join(mascot)


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
