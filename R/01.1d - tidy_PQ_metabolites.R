#Import metabolite data

secmet <- "C:/Users/Thinkpad/ownCloud/EXP/EXP_16_160509_AZD_gluc/Omics_Integration_EXP16/Data/20170824_Krzy_Sec_AZD_Simca.xlsx"
primmet <- "C:/Users/Thinkpad/ownCloud/EXP/EXP_16_160509_AZD_gluc/Omics_Integration_EXP16/Data/Merged_TOXID.csv"
primsamples <- "C:/Users/Thinkpad/ownCloud/EXP/EXP_16_160509_AZD_gluc/Omics_Integration_EXP16/Data/GC_samplelist.xlsx"


tidy_PQ_metabolites <- function(secmet, primmet, primsamples){
  #Import GC sample file from Gudrun
  GC_sampleData <- readxl::read_excel(primsamples) %>% 
    rename(sample_ID = `Sample Name`, GC_ID = `GC-MS_ID(10/90)`) %>%
    select(sample_ID, GC_ID) %>%
    mutate(sample_ID = sub("[[:alnum:]]{1,3}_","",sample_ID),
           sample_ID = sub("/","_",sample_ID),
           group = stringr::str_extract(sample_ID, "\\D{1,}"), #"(Control)|(Glucose)|(AZD_Glu)|(Full_nutrition)"),
           time = stringr::str_extract(sample_ID, "0|(15)|(60)|(240)")) %>%
    group_by(group, time) %>%
    mutate(replicate = seq_along(along.with=sample_ID),
           sample_ID = sub("_\\d{1,2}$","",sample_ID),
           sample_ID = paste(sample_ID, replicate, sep="_"),
           GC_ID = gsub(":","_", GC_ID))%>%
    arrange(sample_ID)
  #Import GC intensities from ToxID
  GC_intData <- readr::read_csv(primmet) %>%
    gather(key=GC_ID, value=intensity,contains("kg")) %>%
    left_join(GC_sampleData) %>%
    mutate(peak_ID = `Compound Name`,
           intensity = as.numeric(intensity))%>%
    select(peak_ID, sample_ID, intensity)
   #Calculate normalization factor from ribitol peak (factor = distance of ribitol in each sample from mean intensity of measurement) 
  GC_intData_IS <- filter(GC_intData,peak_ID=="RIBITOL_217") %>%
    mutate(factor = intensity/mean(intensity)) %>%
    select(sample_ID, factor)
  # Re-calculate intensities based on normalization factor
  GC_intData <- left_join(GC_intData, GC_intData_IS) %>%
    mutate(intensity = intensity/factor) %>%
    select(peak_ID, sample_ID, intensity)
  # Import LC metabolite data from QI (this was already normalized by Mohamed)
  # Append GC intentsity data to the table
  intData <- readxl::read_excel(secmet) %>% 
    slice(-1) %>%
    gather(key=sample_ID, value=intensity,contains("_")) %>%
    mutate(peak_ID = `Sample name`,
           intensity = as.numeric(intensity)) %>%
    select(peak_ID, sample_ID, intensity) %>%
    bind_rows(GC_intData) %>%
    filter(!grepl("(Blank)",sample_ID)) %>%
    mutate(intensity = zero.to.na(intensity))
  # Make sample data table from intData
  # Append GC sample data
  sampleData <- intData %>% 
    select(sample_ID) %>%
    distinct() %>%
    mutate(group = stringr::str_extract(sample_ID, "\\D{1,}"), #"(Control)|(Glucose)|(AZD_Glu)|(Full_nutrition)"),
           time = stringr::str_extract(sample_ID, "0|(15)|(60)|(240)")) %>%
    bind_rows(select(GC_sampleData, sample_ID, group, time)) %>%
    mutate(treatment = sub("_$","",group),
           group = paste(group, time, sep="")) %>%
    filter(!grepl("Blank",sample_ID)) %>%
    distinct()

  df <- list(
    sampleData = sampleData,
    intData = intData,
    peakData = full_join(sampleData,intData))
  
  return(df)
}

