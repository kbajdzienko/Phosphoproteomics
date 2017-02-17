# Two-way ANOVA, for treatment and time variables
# Between subjects.
#
# p.cor - method of p-value correction, can be "holm", "hochberg", "hommel",
# "bonferroni", "BH", "BY", "fdr", "none"

anova2 <- function (df, p.cor = "fdr") {

  # Function, conducting two-way ANOVA/ANCOVA on a subset data frame
  # and getting three p-values
  pv.anova2 <- function(ddf) {
    aov(intensity ~ treatment*time, data = ddf) %>%
      summary() %>%
      unlist() %>%
      .[17:19]
  }

  pv_df <-
    right_join(df$intData, df$sampleData) %>%
    ddply("peak_ID", pv.anova2) %>%
    tbl_df() %>%
    setNames(c("peak_ID", "pv_treatment", "pv_time", "pv_inter")) %>%
    mutate_at(vars(-peak_ID), p.adjust, method = p.cor)

  return(pv_df)
}

# Function plotting two-way ANOVA venn diagram of significantly changing features

plot_venn_anova2 <- function(df, threshold = 0.05, p.cor = "fdr") {
  pv_df <- anova2(df, p.cor)
  venn_list <- lapply(pv_df[-1], function(x) pv_df$peak_ID[x <= 0.05])
  dev.off()
  VennDiagram::venn.diagram(venn_list,
                            category.names = c("Treatment", "Time", "Interaction"),
                            fill = c('yellow', 'purple', 'green'),
                            rotation = 1,
                            main = "Two-way ANOVA: significant features",
                            main.fontface = "bold",
                            main.fontfamily = "sans",
                            main.cex = 1.4,
                            # height = 480,
                            # width = 480,
                            # resolution = 300,
                            # compression = "lzw",
                            lwd = 2,
                            lty = 'blank',
                            cex = 2,
                            # fontface = "bold",
                            fontfamily = "sans",
                            # cat.cex = 0.6,
                            cat.fontface = "bold",
                            cat.fontfamily = "sans",
                            cat.default.pos = "outer",
                            # cat.pos = c(-27, 27, 135),
                            # cat.dist = c(0.055, 0.055, 0.085),
                            filename = NULL) %>%
    grid::grid.draw()
}
