
tidy_MQ_phosphosites <- function(file){
  #Import MQ data that was preprocessed in Perseus (groups added and annotation added)
  df <- read.delim(file,
                   header = F,
                   stringsAsFactors = F,
                   colClasses = "character")
  
  
  # Data about samples
  sampleData <-
    df %>%
    t() %>%
    tbl_df() %>%
    filter(grepl("E",V2))%>%
    select(V1,V3) %>%
    setNames(c("sample_ID", "group")) %>%
    mutate(group = stringr::str_extract(group, "[[:alnum:][-_]]{1,}$"),
           sample_ID2 = sample_ID,
           sample_ID = paste(group, sample_ID, sep="_"))
  
  #Protein Data
  peak.cols <- grep("T$|N$|C$", df[2,])
  
  
  peakData <-
    df %>%
    tbl_df() %>%
    select(peak.cols) %>%
    setNames(.[1, ]) %>%
    slice(-(1:3)) %>%
    group_by(`Unique identifier`)%>%
    mutate(peak_ID = paste0(Protein,"_",`Amino acid`,Position,"_",`Unique identifier`)) %>%
    ungroup() %>%
    select(peak_ID, everything()) 
  
  
  #Intensity Data
  
  int.cols <- grep("Intensity [[:print:]]", df[1,],value = F)
  int.cols2 <- grep("Intensity [[:print:]]", df[1,],value = T)
  
  intData2 <- 
    df %>%
    tbl_df() %>%
    setNames(.[1, ]) %>%
    slice(-(1:3)) %>%
    group_by(`Unique identifier`)%>%
    mutate(peak_ID = paste0(Protein,"_",`Amino acid`,Position,"_",`Unique identifier`)) %>%
    ungroup()
  
  
  
  
  
  intData <- intData2 %>%
    gather(int.cols, key = "sample_ID2", value = "intensity") %>%
    mutate(intensity = as.numeric(intensity)) %>%
    group_by(sample_ID2, peak_ID) %>%
    summarize(intensity = sum(intensity)) %>%
    ungroup()%>%
    left_join(sampleData) %>%
    select(peak_ID, sample_ID, intensity)
  
  df <- list(sampleData = select(sampleData, sample_ID, group),
             peakData = peakData,
             intData = intData)
  
  df$intData$intensity[is.nan(df$intData$intensity)] <- NA
  
  return(df)
}
