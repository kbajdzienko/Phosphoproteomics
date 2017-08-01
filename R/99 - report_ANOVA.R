# One-way ANOVA, for treatment
# Between subjects.
#
# p.cor - method of p-value correction, can be "holm", "hochberg", "hommel",
# "bonferroni", "BH", "BY", "fdr", "none"

anova <- function (df, p.cor = "fdr") {
  
    # Function, conducting two-way ANOVA/ANCOVA on a subset data frame
  # and getting three p-values
  pv.anova <- function(ddf) {
    aov(intensity ~ group, data = ddf) %>%
      summary() %>%
      unlist() %>%
      .[9]
  }
  
    pv_df <-
      right_join(df$intData, df$sampleData) %>%
      plyr::ddply("peak_ID", pv.anova) %>%
      tbl_df() %>%
      setNames(c("peak_ID", "pv_treatment")) %>%
      mutate_at(vars(-peak_ID), p.adjust, method = p.cor)
  

  return(pv_df)
}

