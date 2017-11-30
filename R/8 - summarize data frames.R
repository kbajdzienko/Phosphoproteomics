make_summarized_df <- function(df, data="phospho", summary_function=median){
  ddf <- df 
  
  if(data=="protein"){
    ddf$intData <- right_join(ddf$intData, ddf$sampleData) %>%
      group_by(peak_ID, group) %>%
      mutate(intensity = summary_function(intensity),
             sample_ID = group) %>%
    ungroup() %>% select(peak_ID, sample_ID, intensity)  %>%
      distinct()
  }
  else if(data == "phospho"){
    
    ddf$annIntData <- right_join(ddf$annIntData, ddf$sampleData) %>%
      group_by(ann_ID, group) %>%
      mutate(intensity = summary_function(intensity),
             sample_ID = group) %>%
    ungroup() %>% select(ann_ID, sample_ID, intensity) %>%
      distinct()
  }
  
  ddf$sampleData <- mutate(ddf$sampleData,sample_ID=group) %>% distinct()

  return(ddf)  
}


make_fc_df <- function(df, data="phospho", control="WT"){
  ddf <- df 
  
  if(data=="protein"){
    ddf$intData <- right_join(ddf$intData, ddf$sampleData) %>%
      group_by(peak_ID, time) %>%
      mutate(intensity = log(intensity/mean(intensity[treatment==control]),2)) %>%
      ungroup() %>%
      select(sample_ID, peak_ID, intensity)
      
  }
  else if(data == "phospho"){
    
    ddf$annIntData <- right_join(ddf$annIntData, ddf$sampleData) %>%
      group_by(ann_ID, time) %>%
      mutate(intensity = log(intensity/mean(intensity[treatment==control]),2)) %>%
      ungroup() %>%
      select(sample_ID, ann_ID, intensity)
  }
  
  return(ddf)  
}
