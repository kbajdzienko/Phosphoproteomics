# Example of usage
# Create list of all your tables you want to plot
df_list <-
  list(df_cv,
       df_cv_sum,
       df_NA_cv,
       df_NA_cv_sum) %>%
  setNames(c("df_cv",
             "df_cv_sum",
             "df_NA_cv",
             "df_NA_cv_sum"))

# Plot the fucking histograms
plot_list_hist(df_list, "CV.pdf")

# Function to plot CV histograms to pdf file
plot_list_hist <- function(df_list, file = "CV_report.pdf") {

  # Get the number of columns in dataframes
  n <- unique(sapply(df_list, ncol))
  if (length(n) != 1) stop("Number of columns in data frames is different")

  # Check if all the names are the same between data frames
  name_test <-
    sapply(df_list, function(df) names(df)[-1]) %>%
    apply(1, function(group_ID) length(unique(group_ID)) == 1) %>%
    all()
  if (!name_test) stop("Group names are different in data frames")


  # Open pdf to plot
  pdf(file)
  par(mfrow=c(3,2),
      oma = c(0, 0, 2, 0))

  for (i in 2:n) {
    for (j in 1:length(df_list)) {
      CV <- df_list[[j]][[i]]
      hist(CV,
           xlim = c(0, 2.5),
           breaks = seq(0, 2.5, length.out = 19),
           xlab = "CV",
           main = names(df_list)[j])
      abline(v = median(CV, na.rm = T), col = "red")
    }
    mtext(names(df_list[[1]][i]), outer = TRUE)
    par(mfrow=c(3,2))
  }
  dev.off()
}
