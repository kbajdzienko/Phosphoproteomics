# Plot two samples one against another at scatterplot

plot_2samples <- function(df, sample1, sample2) {
  if (!all(c(sample1, sample2) %in% df$sampleData$sample_ID))
    stop("Samples are not found")

  ggplot(intTable(df), aes_string(x = sample1, y = sample2)) +
    geom_point(shape = 1)
}
