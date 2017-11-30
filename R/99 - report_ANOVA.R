# One-way ANOVA, for treatment
# Between subjects.
#
# p.cor - method of p-value correction, can be "holm", "hochberg", "hommel",
# "bonferroni", "BH", "BY", "fdr", "none"

report_df_anova <- function (df, p.cor = "fdr",psite=TRUE) {
  
  match.arg(p.cor, c("bonferroni", "BH", "BY", "fdr", "none"))
    # Function, conducting one-way ANOVA/ANCOVA on a subset data frame
  # and getting three p-values
  pv.anova <- function(ddf) {
    aov(intensity ~ group, data = ddf) %>%
      summary() %>%
      unlist() %>%
      .[9]
  }
  
  if(psite){
    pv_df <-
      right_join(df$annIntData, df$sampleData) %>%
      plyr::ddply("ann_ID", pv.anova) %>%
      tbl_df() %>%
      setNames(c("ann_ID", "pv_treatment")) %>%
      mutate_at(vars(-ann_ID), p.adjust, method = p.cor)
    }
  else{
    pv_df <-
      right_join(df$intData, df$sampleData) %>%
      plyr::ddply("peak_ID", pv.anova) %>%
      tbl_df() %>%
      setNames(c("peak_ID", "pv_treatment")) %>%
      mutate_at(vars(-peak_ID), p.adjust, method = p.cor)
    }
  
  
  return(pv_df)
}

