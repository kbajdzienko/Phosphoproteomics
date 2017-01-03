#Here are scripts that we used to determine the metrics to assess good quality spectra
#and filter out peak entries that have high likeliyhood of being false positives

#We have calculated CV values within each replicate group and looked at the distribution of
#those values. Very high CV among detected features indicated that many of those should be discarded
df <- tidy_PQI("exp16_peptides_phos.csv")
df <- make_annID(df, "exp16_mascot.csv")
df_cv <- peakCV(df)
df_cv_sum <- annCV(df)
df_NA <- filter_NA(df)
df_NA_cv <- peakCV(df_NA)
df_NA_cv_sum <- annCV(df_NA)




zero.to.na <- function(x) {
  x[x==0] <- NA
  return(x)
}

#Calculate coeficient of variation within replicates of each group
peakCV <- function(df) {
  df_cv <-
    df$intData %>%
    mutate(intensity = zero.to.na(intensity)) %>%
    left_join(df$sampleData) %>%
    group_by(peak_ID, group) %>%
    summarize(CV = sd(intensity, na.rm = T)/mean(intensity, na.rm = T)) %>%
    spread(group, CV)
  return(df_cv)
}

#Sum intensities of duplicate phos site entries in each sample
annCV <- function(df) {
  df_cv_ann <-
    df$intData %>%
    left_join(select(df$peakData, peak_ID, ann_ID)) %>%
    group_by(ann_ID, sample_ID) %>%
    summarize_at(vars(intensity), sum, na.rm = T) %>%
    mutate(intensity = zero.to.na(intensity)) %>%
    left_join(df$sampleData) %>%
    group_by(ann_ID, group) %>%
    summarize(CV = sd(intensity, na.rm = T)/mean(intensity, na.rm = T)) %>%
    spread(group, CV) %>%
    ungroup()
  return(df_cv_ann)
}

#Plot frequency distribution of CV for each group
for (i in 2:7) {
  hist(df_cv_sum[[i]],
       xlim = c(0, 2.5),
       breaks = seq(0, 2.5, length.out = 20),
       main = paste0(gsub("\\W", "_", names(df_cv)[i])))
  par(mfrow=c(2,2))
}


#Is there a correlation between score and CV?
plot(df_cv$'GLU-015', log10(df_cv$Score))
cor(log10(df_cv$Score), df_cv$`AZD/GLU-240`, use = "pairwise.complete.obs")

#How many zero values are in the data frame?
df_sum <- df$intData %>%
  group_by(sample_ID) %>%
  summarize(sum = sum(intensity==0))

#Sum intensities for phos sites from samples
df_sum <-
  df$intData %>%
  left_join(df$peakData) %>%
  filter(grepl("Phos", Modifications)) %>%
  filter(Score > 24) %>%
  group_by(ann_ID, sample_ID) %>%
  summarize_at(vars(intensity), sum) %>%
  spread(sample_ID, intensity)
