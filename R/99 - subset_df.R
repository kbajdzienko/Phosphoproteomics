# Subset by pattern of sample names

subset_df_samples <- function(df, pattern, data="psite") {
  match.arg(data, c("psite","protein"))
  
  if(data=="psite"){
  shitsamples <- filter(df$sampleData, !grepl(pattern, group)) %>% select(sample_ID)
  df$sampleData <- anti_join(df$sampleData, shitsamples)
  df$intData <- anti_join(df$intData, shitsamples)
  df$annIntData <- anti_join(df$annIntData, shitsamples)}
  else if(data=="protein"){
  shitsamples <- filter(df$sampleData, !grepl(pattern, group)) %>% select(sample_ID)
  df$sampleData <- anti_join(df$sampleData, shitsamples)
  df$intData <- anti_join(df$intData, shitsamples)}
  return(df)
}

# Subset by the vector of accessions

subset_df_acc <- function(df, acc_list,data="psite") {
  
  match.arg(data, c("psite","protein", "peak_ID"))
  
  if(data=="psite"){
  df$annData <- filter(df$annData, Accession %in% acc_list)
  df$peakData <- filter(df$peakData, Accession %in% acc_list)
  df$intData <- left_join(df$peakData,df$intData, by = "peak_ID") %>% select(peak_ID, sample_ID, intensity)
  df$annIntData <- left_join(df$annData,df$annIntData, by = "ann_ID") %>% select(ann_ID, sample_ID, intensity)}
  else if(data=="protein"){
    df$peakData <- filter(df$peakData, Accession %in% acc_list)
  df$intData <- left_join(df$peakData,df$intData, by = "peak_ID") %>% select(peak_ID, sample_ID, intensity)}
  else if(data=="peak_ID"){
    df$peakData <- filter(df$peakData, peak_ID %in% acc_list)
    df$intData <- left_join(df$peakData,df$intData, by = "peak_ID") %>% select(peak_ID, sample_ID, intensity)}
  return(df)
}
