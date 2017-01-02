#Here are scripts that we used to determine the metrics to assess good quality spectra
#and filter out peak entries that have high likeliyhood of being false positives


#We have calculated CV values within each replicate group and looked at the distribution of
#those values. Very high CV among detected features indicated that many of those should be discarded




#Plot frequency distribution of CV for each group  
#and save as png file in the working directory

plot_cv <- function(df) {
  df_cv <- cv(df)
  df_cv_sum <- cv_sum(df)
  for (group in names(df_cv)[-1]) {
    png(file = paste0(gsub("\\W", "_", group), "_cv.png"))
    hist(df_cv[[group]], main = paste0(gsub("\\W", "_", group)))
    dev.off()
    png(file = paste0(gsub("\\W", "_", names(df_cv_sum)[i]), "_cv_sum.png"))
    hist(df_cv_sum[[i]], main = paste0(gsub("\\W", "_", group)))
    dev.off()
  }
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

df_NA_score24 <- filter_NA(df, threshold = 0.2)
df_NA_score24 <- filter_score(df_NA_score24, score_threshold = 24)
df_NA_score5 <- filter_NA(df, threshold = 0.2)
df_NA_score5 <- filter_score(df_NA_score5, score_threshold = 5)
df_NA <- filter_NA(df, threshold = 0.2)
df_score24 <- filter_score(df, score_threshold = 24)
