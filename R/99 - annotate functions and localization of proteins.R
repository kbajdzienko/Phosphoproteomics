#Intensity distribution of plastidic proteins

annotate_function_location <- function(df, phospho = T) {
  
if(phospho){
  data <-  df$annData %>%
  mutate(Accession = gsub("_.{1,}",".1", ann_ID)) }
else{
  data <-  df$peakData %>%
  mutate(Accession = gsub(";.{1,}","", Accession))}  

          annotate_data <- function(data, phos){  #this function works on a table that has Accession(AGI) column
          data_annotated <- data %>%
            left_join(select(db$suba, Accession, `location_consensus`))%>%
            left_join(distinct(select(db$uniprot, Accession, `Function [CC]`, `Gene ontology (GO)`, `Protein names`,`Entry name`),Accession,.keep_all=T)) %>%
            mutate(location_simple = location_consensus,
                   location_simple = ifelse(grepl(",",location_simple),"other",gsub("cytosol", "Cytosol", location_simple)),
                   location_simple = gsub("cytosol", "Cytosol", location_simple),
                   location_simple = gsub("Cytosol,nucleus", "Cytosol", location_simple),
                   location_simple = gsub("endoplasmic reticulum", "ER, Golgi", location_simple),
                   location_simple = gsub("extracellular", "other", location_simple),
                   location_simple = gsub("golgi", "ER, Golgi", location_simple),
                   location_simple = gsub("mitochondrion", "Mitochondrion", location_simple, fixed=T),
                   location_simple = gsub("nucleus", "Nucleus", location_simple),
                   location_simple = gsub("peroxisome", "Cytosol", location_simple),
                   location_simple = gsub("plasma membrane", "other", location_simple, fixed=T),
                   location_simple = gsub("plasma membrane,Cytosol", "Cytosol", location_simple),
                   location_simple = gsub("plastid", "Plastid", location_simple, fixed=T),
                   location_simple = gsub("Plastid,Mitochondrion", "Plastid", location_simple),
                   location_simple = gsub("vacuole", "other", location_simple),
                   location_simple = gsub("Cytosol,ER, Golgi", "Cytosol", location_simple),
                   location_simple = gsub("ER, Golgi,ER, Golgi", "ER, Golgi", location_simple),
                   location_simple = gsub("Plastid,Cytosol", "other", location_simple, fixed=T),
                   encoded = ifelse(phos,if_else(grepl("ATC",ann_ID),"Plastid","Nucleus"),
                                             if_else(grepl("ATC",peak_ID),"Plastid","Nucleus")),
                   function_simple = ifelse(grepl("(translation)|(ribosome)",`Gene ontology (GO)`), 
                                            
                                            ifelse(grepl("(ribosomal protein)|(initiation factor)",`Protein names`),
                                                   stringr::str_extract(`Protein names`,"(ribosomal protein)|(initiation factor)"), 
                                                   "translation other"),
                                            "other"), `Entry name` = gsub("_ARATH","",`Entry name`))}

  if(phospho){df$annData <- annotate_data(data, phos = T)}else{
            df$peakData <- annotate_data(data, phos = F)}        
          
  return(df)
}



