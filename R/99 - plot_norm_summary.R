# Plot two summary plot, one before normalization, one after
# for each plot top is box plot, bottom is a density plot

plot_norm_summary <- function(df.pre, df.post, mode) {
  match.arg(mode, c("peak", "sample"))

  layout(matrix(c(1,1,1,2,3,3,3,4), 4, 2, byrow = FALSE))

  # fig 1
  op <- par(mar = c(5,7,4,0), xaxt = "s");
  plot_boxint(df.pre, mode)
  mtext("Before Normalization",3, 1);

  # fig 2
  op <- par(mar = c(7,7,0,0), xaxt = "s");
  plot_density(df.pre, mode)
  mtext("Density", 2, 5);
  mtext("Intensity", 1, 5);

  # fig 3
  op <- par(mar = c(5,7,4,2), xaxt = "s");
  plot_boxint(df.post, mode)
  mtext("After Normalization", 3, 1);

  # fig 4
  op <- par(mar = c(7,7,0,2), xaxt = "s");
  plot_density(df.post, mode)
  mtext("Normalized Intensity", 1, 5);

}


# Plot density plot
plot_density <- function(df, mode) {
  match.arg(mode, c("peak", "sample"))
  df$intData %>%
    group_by_(paste0(mode, "_ID")) %>%
    summarize(intensity = mean(intensity, na.rm = TRUE)) %>%
    .$intensity %>%
    density(na.rm = TRUE) %>%
    plot(col = 'darkblue',
         las = 2, lwd = 2,
         main = "", xlab = "", ylab = "")
}

# Plot boxplot of random subset 40 peaks/samples
# since there may be too many compounds/samples

plot_boxint <- function(df, mode) {
  match.arg(mode, c("peak", "sample"))

  if (mode == "peak") {
    intMatrix <- t(intMatrix(df))
    # Use for boxplots only peaks without NAs
    names.vec <-
      df$intData %>%
      group_by(peak_ID) %>%
      filter(!any(is.na(intensity))) %>%
      .$peak_ID
  } else if (mode == "sample") {
    intMatrix <- intMatrix(df)
    names.vec <- colnames(intMatrix)
  }

  size <- min(40, length(names.vec))
  names.vec <- sample(names.vec, size = size, replace = FALSE)

  intMatrix %>%
    .[, names.vec] %>%
    boxplot(names = substr(names.vec, 1, 12),
            ylim = range(intMatrix, na.rm = TRUE),
            las = 2,
            col = "lightgreen",
            horizontal = T)
}


