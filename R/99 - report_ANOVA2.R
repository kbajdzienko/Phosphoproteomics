# Two-way ANOVA, for treatment and time variables
# Between subjects.
# Default is type III ANOVA for unbalanced data sets (unequal numbers of observations for each level of a factor).
# It is okay to use type III for balanced data sets as well since the factors are orthogonal, 
# and types I, II and III all give the same results.

# Alternative is type II ANOVA, which should be used when type III will give no interaction between factors. 

# p.cor - method of p-value correction, can be "holm", "hochberg", "hommel",
# "bonferroni", "BH", "BY", "fdr", "none"

anova2 <- function (df, p.cor = "fdr", data="phospho", type = "I") {

  match.arg(data, c("protein", "phospho"))
  match.arg(type, c("I", "II", "III"))
  
  # Function, conducting two-way ANOVA/ANCOVA on a subset data frame
  # and getting three p-values
  pv.anovaIII <- function(ddf) {
    car::Anova(lm(intensity ~ treatment*time, data = ddf, contrasts = list(treatment=contr.sum, time=contr.sum)),
               type=3)%>%
      unlist() %>%
      .[16:18]
  }
  pv.anovaII <- function(ddf) {
    car::Anova(lm(intensity ~ treatment*time, data = ddf),type=2) %>%
      unlist() %>%
      .[13:15]
  }
  
  pv.anovaI <- function(ddf) {
    aov(intensity~treatment*time, data=ddf) %>%
      summary() %>%
      unlist()%>%
      .[17:19]
  }
  
if(data == "protein"){
  pv_df <-
    right_join(df$intData, df$sampleData) %>%
    mutate(time = as.factor(time)) %>%
    plyr::ddply("peak_ID", ifelse(type=="III", pv.anovaIII, ifelse(type=="II", pv.anovaII, ifelse(type=="I", pv.anovaI, pv.anovaIII)))) %>%
    tbl_df() %>%
    setNames(c("peak_ID", "pv_treatment", "pv_time", "pv_interaction")) %>%
    mutate_at(vars(-peak_ID), p.adjust, method = p.cor)
}
  
else if (data == "phospho"){
  pv_df <-
    right_join(df$annIntData, df$sampleData) %>%
    mutate(time = as.factor(time)) %>%
    plyr::ddply("ann_ID", ifelse(type=="III", pv.anovaIII, ifelse(type=="II", pv.anovaII, ifelse(type=="I", pv.anovaI, pv.anovaIII)))) %>%
    tbl_df() %>%
    setNames(c("ann_ID", "pv_treatment", "pv_time", "pv_interaction")) %>%
    mutate_at(vars(-ann_ID), p.adjust, method = p.cor)
}
  return(pv_df)
}

# Function plotting two-way ANOVA venn diagram of significantly changing features

plot_venn_anova2 <- function(df, threshold = 0.05, p.cor = "fdr", data="phospho", type="I") {
  match.arg(data, c("protein", "phospho"))
  
  if (data=="protein"){
    pv_df <- anova2(df, p.cor, data="int", type=type)
  venn_list <- lapply(pv_df[-1], function(x) pv_df$peak_ID[x <= threshold])
  }
  else if (data=="phospho"){
    pv_df <- anova2(df, p.cor, type=type)
    venn_list <- lapply(pv_df[-1], function(x) pv_df$ann_ID[x <= threshold])}
  
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
