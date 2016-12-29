library(plyr)
library(dplyr)
library(tidyr)

make_annID <- function(df, mascot_file) {
  
  # Remove everything but Sequence, Accession, start and end position
  # from mascot output file
  mascot <- 
    read.csv(mascot_file,
             skip = 71,
             header = T, sep = ",",
             stringsAsFactors = F,
             colClasses = "character") %>%
    tbl_df() %>%
    select(prot_acc, pep_seq, pep_start, pep_end) %>%
    mutate_each(funs(as.integer), pep_start, pep_end) %>%
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
  
  # Function looking up aminoacid and merging it with position number in protein   
  aa_merge_pos <- function(sequence, pep_start, phosph_pos) {
    paste0(unlist(strsplit(sequence, ""))[phosph_pos],
           pep_start + phosph_pos,
           collapse = "_")
  }

  # Merge accession and aminoacid+position
  ann_ID <- paste(peakData$Accession,
                  mapply(aa_merge_pos, peakData$Sequence, peakData$pep_start, phosph_pos),
                  sep = "_")
  
  # Add new ann_ID to annData table 
  df$peakData <- 
    peakData %>%
    mutate(ann_ID = ann_ID) %>%
    select(peak_ID:Score, ann_ID, everything())
  
  return(df)
}