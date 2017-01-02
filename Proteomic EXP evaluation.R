#Function to return number of unique prot IDs per LC run
pep_per_run <- function(mascot_file) {
  df <- 
    read.csv(mascot_file, skip = 71, header = T, sep = ",",stringsAsFactors = F) %>%
    select(prot_acc, pep_scan_title) %>%
    mutate(pep_scan_title = substr(pep_scan_title, 1, 10)) %>%  #would be good to be able to 
    distinct() %>%
    group_by(pep_scan_title) %>%
    tally() %>%
    return()
}


#Plot histograms of coeficients of variation in each replicate group in EXP
#Requires peptide tables to be processed by tidy_PQI function
#plot_list_hist_cv(file = "CV.pdf", df, df_NA)


# Function to plot histograms of several data frames' columns to pdf file
plot_list_hist <- function(df_list, file) {
  
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


# Function to plot CV histograms to pdf file
plot_list_hist_cv <- function(..., file = "CV_report.pdf") {
  df_names <- as.character(substitute(...()))
  df_list <- c(sapply(list(...), function(df) list(peakCV(df), annCV(df))))
  names(df_list) <- c(sapply(df_names,
                             paste0, c("_cv", "_cv_sum")))
  plot_list_hist(df_list, file)
}
