# Two-way ANOVA, for treatment and time variables
# Between subjects.
#
# p.cor - method of p-value correction, can be "holm", "hochberg", "hommel",
# "bonferroni", "BH", "BY", "fdr", "none"

anova2 <- function (df, p.cor = "fdr", data="annInt") {

  match.arg(data, c("int", "annInt"))
  
  # Function, conducting two-way ANOVA/ANCOVA on a subset data frame
  # and getting three p-values
  pv.anova2 <- function(ddf) {
    aov(intensity ~ treatment*time, data = ddf) %>%
      summary() %>%
      unlist() %>%
      .[17:19]
  }
if(data == "int"){
  pv_df <-
    right_join(df$intData, df$sampleData) %>%
    plyr::ddply("peak_ID", pv.anova2) %>%
    tbl_df() %>%
    setNames(c("peak_ID", "pv_treatment", "pv_time", "pv_inter")) %>%
    mutate_at(vars(-peak_ID), p.adjust, method = p.cor)
}
else if (data == "annInt"){
  pv_df <-
    right_join(df$annIntData, df$sampleData) %>%
    plyr::ddply("ann_ID", pv.anova2) %>%
    tbl_df() %>%
    setNames(c("ann_ID", "pv_treatment", "pv_time", "pv_inter")) %>%
    mutate_at(vars(-ann_ID), p.adjust, method = p.cor)
  
}
  return(pv_df)
}

# Function plotting two-way ANOVA venn diagram of significantly changing features

plot_venn_anova2 <- function(df, threshold = 0.05, p.cor = "fdr", data="annInt") {
  match.arg(data, c("int", "annInt"))
  pv_df <- anova2(df, p.cor)
  if (data=="int"){
  venn_list <- lapply(pv_df[-1], function(x) pv_df$peak_ID[x <= 0.05])
  }
  else if (data=="annInt"){
    venn_list <- lapply(pv_df[-1], function(x) pv_df$ann_ID[x <= 0.05])}
  
  venn <- VennDiagram::venn.diagram(venn_list,
                            category.names = c("Treatment", "Time", "Interaction"),
                            fill = c("#DBC094", "#8BD3B0", "#A2C8EC"),
                            alpha=0.6,
                            lty = 'blank',
                            cex = 1.1,
                            fontfamily = "serif",
                            cat.cex = 1.1,
                            cat.default.pos = "outer",
                            # cat.pos = c(-27, 27, 135),
                            # cat.dist = c(0.055, 0.055, 0.085),
                            # fontface = "bold",
                            #cat.fontface = "bold",
                            #cat.fontfamily = "serif",
                            #main = "Two-way ANOVA: significant features",
                            #main.fontface = "bold",
                            #main.fontfamily = "serif",
                            #main.cex = 1.1,
                            # height = 480,
                            # width = 480,
                            # resolution = 300,
                            # compression = "lzw",
                            #lwd = 2,
                            filename = NULL,
                            rotation = 1)
  
  gridExtra::grid.arrange(grid::gTree(children=venn) )
               #top=textGrob("Two-way ANOVA: significant features", gp=gpar(fontsize=20,fontfamily="serif"))
                 
}
