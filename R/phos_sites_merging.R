splitSites<- function (ddf) {
  ann <- unlist(strsplit(ddf$ann_ID, "_"))
  ann_ID_new <- paste(ann[1], ann[-1], sep = "_")
  ddf <-
    ddf[rep(1, length(ann_ID_new)), ] %>%
    mutate(ann_ID = ann_ID_new)
  #paste0(ddf[1,1], seq(length(ann_ID)), sep = "_")
  return(ddf)
}


df <- tidy_PQI("exp19-3 peptides.csv")
df <- make_annID(df, "exp19-3-mascot.csv", skip = 0)

mult_sites <-
  df$peakData %>%
  filter(grepl("_[[:alnum:]]+_", ann_ID)) %>%
  ddply("peak_ID", splitSites) %>%
  tbl_df()

df$peakData <-
  df$peakData %>%
  filter(!grepl("_[[:alnum:]]+_", ann_ID)) %>%
  bind_rows(mult_sites)

sumANN(df)

sumANN <- function(df) {
  df$intData <-
    df$intData %>%
    left_join(select(df$peakData, peak_ID, ann_ID)) %>%
    group_by(ann_ID, sample_ID) %>%
    summarize_at(vars(intensity), sum, na.rm = T) %>%
    mutate(intensity = zero.to.na(intensity)) %>%
    ungroup()
  return(df)
}
