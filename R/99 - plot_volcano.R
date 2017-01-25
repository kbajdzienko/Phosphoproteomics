# Volcano plot

plot_peak_volcano <- function(df, group1, group2,
                         method = "t.test",
                         FC.threshold = 2,
                         p.threshold = 0.05) {

  match.arg(method, c("t.test", "wilcox.test"))

  if (!all(c(group1, group2) %in% df$sampleData$group))
    stop("Selected group names are not found")

  volcano <-
    df$sampleData %>%
    filter(group %in% c(group1, group2)) %>%
    left_join(df$intData) %>%
    group_by(peak_ID) %>%
    summarize(p.value = do.call(method, list(intensity[group == group1], intensity[group == group2]))$p.value,
              FC = mean(intensity[group == group1], na.rm = T)/mean(intensity[group == group2], na.rm = T)) %>%
    mutate(log10p = -log10(p.value),
           logFC = log2(FC)) %>%
    mutate(signif = abs(p.value) < p.threshold & abs(logFC) > log2(FC.threshold))

  ggplot(volcano, aes(x = logFC, y = log10p, color = signif)) +
    geom_point(size = 1.5) +
    scale_color_manual(values = c("black", "red")) +
    #geom_vline(xintercept = c(-1, 1)*log2(FC.threshold), colour = "grey", linetype = "longdash") +
    scale_y_continuous('-log10(p)') +
    scale_x_continuous('log2(FC)', limits = c(-1, 1)*max(abs(volcano$logFC))) +
    guides(colour = FALSE) +
    labs(title = paste(group1, "/", group2)) +
    theme(
      panel.background = element_rect(fill = 'white', colour = "black"),
      panel.grid.minor = element_line(colour = "grey"),
      strip.background = element_rect(fill = 'gray95'),
      strip.text.x = element_text(colour = c("#00B0F6"), size = 14),
      axis.text.x = element_text(size = 10, colour="black"),
      axis.text.y = element_text(size = 10, colour="black"),
      axis.ticks = element_line(colour = "black"),
      axis.title.x = element_text(size = 15, vjust = -0.4),
      axis.title.y = element_text(size = 15, vjust = 0.3),
      title = element_text(size = 18, vjust = 1.5)
    )
}

plot_ann_volcano <- function(df, group1, group2,
                             method = "t.test",
                             FC.threshold = 2,
                             p.threshold = 0.05) {
  
  match.arg(method, c("t.test", "wilcox.test"))
  
  if (!all(c(group1, group2) %in% df$sampleData$group))
    stop("Selected group names are not found")
  
  volcano <-
    df$sampleData %>%
    filter(group %in% c(group1, group2)) %>%
    left_join(df$annIntData) %>%
    group_by(ann_ID) %>%
    summarize(p.value = do.call(method, list(intensity[group == group1], intensity[group == group2]))$p.value,
              FC = mean(intensity[group == group1], na.rm = T)/mean(intensity[group == group2], na.rm = T)) %>%
    mutate(log10p = -log10(p.value),
           logFC = log2(FC)) %>%
    mutate(signif = abs(p.value) < p.threshold & abs(logFC) > log2(FC.threshold))
  
  ggplot(volcano, aes(x = logFC, y = log10p, color = signif)) +
    geom_point(size = 1.5) +
    scale_color_manual(values = c("black", "red")) +
    #geom_vline(xintercept = c(-1, 1)*log2(FC.threshold), colour = "grey", linetype = "longdash") +
    scale_y_continuous('-log10(p)') +
    scale_x_continuous('log2(FC)', limits = c(-1, 1)*max(abs(volcano$logFC))) +
    guides(colour = FALSE) +
    labs(title = paste(group1, "/", group2)) +
    theme(
      panel.background = element_rect(fill = 'white', colour = "black"),
      panel.grid.minor = element_line(colour = "grey"),
      strip.background = element_rect(fill = 'gray95'),
      strip.text.x = element_text(colour = c("#00B0F6"), size = 14),
      axis.text.x = element_text(size = 10, colour="black"),
      axis.text.y = element_text(size = 10, colour="black"),
      axis.ticks = element_line(colour = "black"),
      axis.title.x = element_text(size = 15, vjust = -0.4),
      axis.title.y = element_text(size = 15, vjust = 0.3),
      title = element_text(size = 18, vjust = 1.5)
    )
}
