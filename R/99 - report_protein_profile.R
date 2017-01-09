# Plot MSstat protein profile
# Example:
# plot_acc(df, "AT3G60600.1")

plot_acc <- function(df, protein) {
  require(ggplot2)
  y.limup <- ceiling(max(log2(df$intData$intensity), na.rm=TRUE) + 3)
  y.limdown <- -1
  tempGroupName <- df$sampleData # unique(datafeature[, c("GROUP_ORIGINAL", "RUN")])
  groupAxis <- as.numeric(xtabs(~group, tempGroupName))
  cumGroupAxis <- cumsum(groupAxis)
  lineNameAxis <- cumGroupAxis[-nlevels(as.factor(df$sampleData$group))]
  groupName <- data.frame(sample_ID = c(0, lineNameAxis) + groupAxis / 2 + 0.5,
                          intensity = rep(y.limup - 1, length(groupAxis)),
                          Name = levels(as.factor(df$sampleData$group)))
  data <-
    logTransform(df)$intData %>%
    left_join(df$sampleData) %>%
    left_join(df$peakData) %>%
    filter(Accession == protein)
  data.sum <-
    data %>%
    group_by(sample_ID) %>%
    summarize(intensity = median(intensity))

  ptempall <-
    ggplot() +
    geom_point(size = 1.5, color = "grey")+
    geom_line(data = filter(data, Accession == protein),
              aes(x = sample_ID, y = intensity, group = peak_ID),
              size = 0.5, color = "lightgrey") +
    geom_line(data = data.sum,
              aes(x = sample_ID, y = intensity, group = 1),
              size = 0.5, color = "darkred", linetype = "longdash") +
    scale_y_continuous('Log2-intensities', limits = c(y.limdown, y.limup)) +
    scale_x_discrete('MS runs', breaks = cumGroupAxis) +
    geom_vline(xintercept = lineNameAxis + 0.5, colour = "grey", linetype = "longdash") +
    geom_text(data = groupName,
              aes(x = sample_ID, y = intensity,  label = Name),
              size = 4, angle = 0, color = "black") +
    labs(title = protein)+
    theme(
      panel.background = element_rect(fill = 'white', colour = "black"),
      legend.key = element_rect(fill = 'white', colour = 'white'),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = 'gray95'),
      strip.text.x = element_text(colour = c("#00B0F6"), size = 14),
      axis.text.x = element_text(size = 10, colour="black"),
      axis.text.y = element_text(size = 10, colour="black"),
      axis.ticks = element_line(colour = "black"),
      axis.title.x = element_text(size = 15, vjust = -0.4),
      axis.title.y = element_text(size = 15, vjust = 0.3),
      title=element_text(size = 18, vjust = 1.5),
      legend.position = "top",
      legend.text = element_text(size = 7),
      legend.title = element_blank()
    )

  print(ptempall)
}
