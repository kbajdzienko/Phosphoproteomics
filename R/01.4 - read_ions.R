# Function reading ions file

read.ions <- function(ions_file) {
  ions <- read.csv(ions_file,
                   skip = 2,
                   header = T, sep = ",",
                   stringsAsFactors = F,
                   colClasses = "character")
  ions <-
    tbl_df(ions) %>%
    select(1:12) %>%
    setNames(c("ion_ID", "RT", "Charge", "mz",
               "Neutral_mass", "Mass_error",
               "Mass_error_ppm",
               "Score", "Sequence", "Modifications",
               "Accession", "All_accessions")) %>%
    mutate_at(vars(RT, mz:Score), as.numeric) %>%
    mutate_at(vars(ion_ID, Charge), as.integer)

  ions$All_accessions <-
    ions$All_accessions %>%
    strsplit(";") %>%
    sapply(sort) %>%
    sapply(paste, collapse = ";")

  return(ions)
}
