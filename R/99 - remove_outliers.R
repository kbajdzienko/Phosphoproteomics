remove_outliers <- function(df){
  df2<-df
test <- right_join(df2$sampleData,df2$intData) %>%
  group_by(peak_ID, group) %>%
  mutate(conf1 = boxplot.stats(intensity)$conf[1],
         conf2 = boxplot.stats(intensity)$conf[2]) %>%
  ungroup() %>% group_by(peak_ID, sample_ID) %>% 
  mutate(intensity = ifelse((intensity<conf1-(conf1*0.05)|intensity>conf2+(conf2*0.05)),NA,intensity)) %>%
  select(sample_ID, peak_ID, intensity) %>% 
  ungroup()
df2$intData <- test  
  return(df2)
}

remove_outliers_phospho <- function(df){
  df2<-df
  test <- right_join(df2$sampleData,df2$annIntData) %>%
    group_by(ann_ID,group) %>%
    mutate(conf1 = boxplot.stats(intensity)$conf[1],
           conf2 = boxplot.stats(intensity)$conf[2],
           intensity = ifelse((intensity<conf1|intensity>conf2),NA,intensity)) %>%
    select(sample_ID, ann_ID, intensity) %>% 
    ungroup()
  df2$intData <- test  
  return(df2)
}
