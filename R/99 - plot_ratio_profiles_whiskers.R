# Plot MSstat protein profile
# Example:
# plot_acc(df, "AT3G60600.1")
#
# Plot difference between means of log2 intensity betwee
# +- standard deviation of log2 intensities

plot_protein_profile_ratio <- function(df,
                                       protein="AT3G60600.1",
                                       smooth = TRUE,
                                       whiskers = TRUE,
                                       treatment_group="AZD",
                                       control_group="DMSO",
                                       data="p-site") {
  require(ggplot2)
  match.arg(data, c("protein", "p-site"))
  # Log transform if not
  if (!is.log(df$intData$intensity)) df <- logTransform(df)
  
  if (data == "p-site") {
  data <-
    df$annIntData %>%
    left_join(df$sampleData) %>%
    group_by(ann_ID, treatment, time) %>%
    summarise(int_mean = mean(intensity, na.rm = T),
              int_sd = sd(intensity, na.rm = T)) %>%
    group_by(ann_ID, time) %>%
    summarize(ratio = int_mean[treatment == treatment_group]-int_mean[treatment == control_group],
              int_sd = int_sd[treatment == treatment_group] + int_sd[treatment == control_group]) %>%
    ungroup() %>%
    left_join(df$annData) %>%
    filter(Accession == protein) %>%
    mutate(ann_ID = gsub("^[[:alnum:]]+.*[[:alnum:]]*_", "", ann_ID)) %>%
    mutate(Position = as.integer(gsub("[[:alpha:]]+", "", ann_ID))) %>%
    arrange(Position)
  
  # Reorder levels of factor ann_ID by position number for plot legend
  data$ann_ID <- reorder(data$ann_ID, data$Position)
  
  
  # Plot
  ggplot(data, aes(x = jitter(time, 2), y = ratio, group = ann_ID, color = ann_ID,
                   ymin = ratio - int_sd, ymax = ratio + int_sd)) +
    do.call(paste0("geom_", if_else(smooth, "smooth", "line")), list(size = 0.5)) +
    geom_errorbar(alpha = 0.5*whiskers, width = 5) +
    scale_y_continuous('Log2-ratio') +
    scale_x_continuous('Time, min') +
    labs(title = protein) +
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
      #legend.position = "top",
      legend.text = element_text(size = 7),
      legend.title = element_blank()
    )
 
  }
  else if (data == "protein") {
    data <-
      df$intData %>%
      left_join(df$sampleData) %>%
      group_by(peak_ID, treatment, time) %>%
      summarise(int_mean = mean(intensity, na.rm = T),
                int_sd = sd(intensity, na.rm = T)) %>%
      group_by(peak_ID, time) %>%
      summarize(ratio = int_mean[treatment == treatment_group]-int_mean[treatment == control_group],
                int_sd = int_sd[treatment == treatment_group] + int_sd[treatment == control_group]) %>%
      ungroup() %>%
      left_join(df$peakData) %>%
      filter(Accession == protein) %>%
      mutate(peak_ID = gsub("^[[:alnum:]]+.*[[:alnum:]]*_", "", peak_ID)) %>%
      mutate(Position = as.integer(gsub("[[:alpha:]]+", "", peak_ID))) %>%
      arrange(Position)
 
  # Reorder levels of factor peak_ID by position number for plot legend
  #data$peak_ID <- reorder(data$peak_ID, data$Position)


  # Plot
  ggplot(data, aes(x = jitter(time, 2), y = ratio, group = peak_ID, color = peak_ID,
                   ymin = ratio - int_sd, ymax = ratio + int_sd)) +
    do.call(paste0("geom_", if_else(smooth, "smooth", "line")), list(size = 0.5)) +
    geom_errorbar(alpha = 0.5*whiskers, width = 5) +
    scale_y_continuous('Log2-ratio') +
    scale_x_continuous('Time, min') +
    labs(title = protein) +
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
      #legend.position = "top",
      legend.text = element_text(size = 7),
      legend.title = element_blank()
    )
}
}