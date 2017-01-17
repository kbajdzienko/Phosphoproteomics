# Add aAll_accessions column from ions PQI file
# Apply after make_annID to fix I/L issue first!
make_all_acc <- function(df, ions_file) {
  ions <- read.ions(ions_file) %>%
    select(ion_ID, Sequence, Accession, All_accessions) %>%
    distinct()
  df$peakData <-
    df$peakData %>%
    mutate(ion_ID = gsub(",#.+$", "", Decon_ions)) %>%
    mutate(ion_ID = as.integer(gsub("#", "", ion_ID))) %>%
    left_join(ions) %>%
    select(-Decon_ions, -ion_ID)
  return(df)
}
