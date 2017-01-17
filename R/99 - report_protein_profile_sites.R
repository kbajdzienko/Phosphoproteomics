# Plot MSstat general intensity profile accross the samples for quality check

plot_acc_ann <- function(df, protein) {
  require(ggplot2)
  y.limup <- ceiling(max(log2(df$annIntData$intensity), na.rm=TRUE) + 3)
  y.limdown <- -1
  tempGroupName <- df$sampleData # unique(datafeature[, c("GROUP_ORIGINAL", "RUN")])
  groupAxis <- as.numeric(xtabs(~group, tempGroupName))
  cumGroupAxis <- cumsum(groupAxis)
  lineNameAxis <- cumGroupAxis[-nlevels(as.factor(df$sampleData$group))]
  groupName <- data.frame(sample_ID = c(0, lineNameAxis) + groupAxis / 2 + 0.5,
                          intensity = rep(y.limup - 1, length(groupAxis)),
                          Name = levels(as.factor(df$sampleData$group)))
  data <-
    df$annIntData %>%
    mutate(intensity = log2(intensity)) %>%
    left_join(df$sampleData) %>%
    left_join(df$annData) %>%
    filter(Accession == protein) %>%
    mutate(ann_ID = gsub("^[[:alnum:]]+.*[[:alnum:]]*_", "", ann_ID)) %>%
    mutate(Position = as.integer(gsub("[[:alpha:]]+", "", ann_ID))) %>%
    arrange(Position)
  # data.sum <-
  #   data %>%
  #   group_by(sample_ID) %>%
  #   summarize(intensity = median(intensity, na.rm = T))

  ptempall <-
    ggplot() +
    geom_point(size = 1.5, color = "grey")+
    geom_line(data = filter(data, Accession == protein),
              aes(x = sample_ID, y = intensity, group = ann_ID, color = ann_ID),
              size = 0.5) +
    # geom_line(data = data.sum,
    #           aes(x = sample_ID, y = intensity, group = 1),
    #           size = 0.5, color = "darkred", linetype = "longdash") +
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
