# Creata ann_ID in format "AT3G47070.1_T65_T66"
# Accession, phosphorilated amino acid, number

make_annID <- function(df, mascot_file, skip = 71) {

  # Remove everything but Sequence, Accession, start and end position
  # and number of missed cleaveges from mascot output file
  mascot <-
    read.csv(mascot_file,
             skip = skip,
             header = T, sep = ",",
             stringsAsFactors = F,
             colClasses = "character") %>%
    tbl_df() %>%
    select(prot_acc, pep_seq, pep_start, pep_end, pep_miss) %>%
    mutate_each(funs(as.integer), pep_start, pep_end, pep_miss) %>%
    rename(Accession = prot_acc, Sequence = pep_seq) %>%
    distinct()

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
