# Plot MSstat protein profile
# Example:
# plot_acc(df, "AT3G60600.1")

plot_protein_profile_ratio <- function(df, protein ="AT3G60600.1", 
                                       treatment_group="AZD", control_group="DMSO") {
  require(ggplot2)

  # Log transform if not
  #if (!is.log(df$intData$intensity)) df <- logTransform(df)

  data <-
    df$annIntData %>%
    left_join(df$sampleData) %>%
    group_by(ann_ID, treatment, time) %>%
    summarize(int_mean = median(intensity, na.rm = T)) %>%
    group_by(ann_ID, time) %>%
    summarize(ratio = int_mean[treatment == treatment_group]/int_mean[treatment == control_group]) %>%
    ungroup() %>%
    #left_join(df$annData) %>%
    filter(grepl(protein, ann_ID)) %>%
    #mutate(ann_ID = gsub("^[[:alnum:]]+.*[[:alnum:]]*_", "", ann_ID)) %>%
    mutate(Position = as.integer(gsub("[[:alpha:]]+", "", ann_ID)),
           time = as.integer(time)) %>%
    arrange(Position)

  # Reorder levels of factor ann_ID by position number for plot legend
  data$ann_ID <- reorder(data$ann_ID, data$Position)


  # Plot
  ggplot(data, aes(x = time, y = ratio, group = ann_ID, color = ann_ID)) +
    geom_smooth()+
    scale_y_continuous('Log2-ratio') +
    scale_x_continuous('Time, min') +
    labs(title = protein) +
      theme_minimal()+
      scale_color_brewer(palette = "Set1")
    # theme(
    #   panel.background = element_rect(fill = 'white', colour = "black"),
    #   legend.key = element_rect(fill = 'white', colour = 'white'),
    #   panel.grid.minor = element_blank(),
    #   strip.background = element_rect(fill = 'gray95'),
    #   strip.text.x = element_text(colour = c("#00B0F6"), size = 14),
    #   axis.text.x = element_text(size = 10, colour="black"),
    #   axis.text.y = element_text(size = 10, colour="black"),
    #   axis.ticks = element_line(colour = "black"),
    #   axis.title.x = element_text(size = 15, vjust = -0.4),
    #   axis.title.y = element_text(size = 15, vjust = 0.3),
    #   title=element_text(size = 18, vjust = 1.5),
    #   #legend.position = "top",
    #   legend.text = element_text(size = 7),
    #   legend.title = element_blank()
    # )
}
