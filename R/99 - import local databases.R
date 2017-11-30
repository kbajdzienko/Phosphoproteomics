#Import databases

import_localdb <- function(up_path = "C:/Users/Thinkpad/ownCloud/EXP/UniProt_Ara_2017Oct.xlsx", 
                               suba_path = "C:/Users/Thinkpad/ownCloud/EXP/Suba4-2017-10-24_18-0.csv"){
  
  UniProt_Ara_2017Oct <- readxl::read_excel(up_path,sheet = 1) %>%
    mutate(Accession = stringr::str_extract(`Gene names  (ordered locus )`,"At.g\\d{5}"),
           Accession = paste(Accession, "1", sep = "."),
           Accession = gsub("t","T", Accession),
           Accession = gsub("g","G", Accession)) %>%
    filter(Accession != "NA.1")
  
  Suba4_2017 <- read_csv(suba_path) %>% select(Accession, location_consensus)
  
  db <- list(uniprot = UniProt_Ara_2017Oct,
            suba = Suba4_2017)
  return(db)
}
