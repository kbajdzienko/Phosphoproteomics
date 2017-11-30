filter_df_cv <- function(df, data="psite", min_groups=5, min_CV=0.7){
  match.arg(data,c("psite","protein"))
  if(data=="psite"){
    df_cv <-
      df$annIntData %>%
      mutate(intensity = zero.to.na(intensity)) %>%
      left_join(df$sampleData) %>%
      group_by(ann_ID, group) %>%
      summarize(CV = sd(intensity, na.rm = T)/mean(intensity, na.rm = T)) %>%
      ungroup()%>%
      filter(CV<min_CV) %>%
      group_by(ann_ID)%>%
      count() %>%
      filter(n>min_groups) %>%
      .$ann_ID
    df_cv <- gsub("_.{1,}","", df_cv)
    df <- subset_df_acc(df, df_cv)
    df$intData <- distinct(df$intData)
    df$annIntData <- distinct(df$annIntData)
  }
  else{
    df_cv <-
      df$intData %>%
      mutate(intensity = zero.to.na(intensity)) %>%
      left_join(df$sampleData) %>%
      group_by(peak_ID, group) %>%
      summarize(CV = sd(intensity, na.rm = T)/mean(intensity, na.rm = T)) %>%
      ungroup()%>%
      filter(CV<min_CV) %>%
      group_by(peak_ID)%>%
      count() %>%
      filter(n>min_groups) %>%
      .$peak_ID
    df <- subset_df_acc(df, df_cv, data="peak_ID")
    df$intData <- distinct(df$intData)
  }
  return(df)
}
