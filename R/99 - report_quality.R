# As for the EXP summary:
# Getting 2 tables with numbers:

# And plots:
# Abundance plots IN msstats: box plot on page 14 and profile plots on page 15
# I think those profile plots are quite informative (and colorful).



mascot <- read.csv(mascot_file,
                   skip = skip,
                   header = T,
                   sep = ",",
                   stringsAsFactors = F)
mascot <- tbl_df(mascot)

# All identified proteins
prot.num <- summarize(mascot, Proteins = n_distinct(prot_acc))

# Phosphorylated proteins
prot.phos.num <-
  filter(mascot, grepl("Phosp", pep_var_mod)) %>%
  summarize(PhosphoProteins = n_distinct(prot_acc))

# Identified unique proteins within each run -- CV
prot.per.run.cv <-
  mascot %>%
  mutate(sample_ID = sapply(strsplit(pep_scan_title, "\\."), '[', 1)) %>%
  distinct(sample_ID, prot_acc) %>%
  count(sample_ID) %>%
  summarize(UniqueProteinsCV = sd(n)/mean(n))

# Function to count number of phosphorilation sites
countSites <- function(pep_start, pep_pos) {
  sites <-
    strsplit(pep_pos, "") %>%
    sapply(function(x) which(x == 3)) %>%
    mapply(function(start, pos) start + pos - 3, pep_start, .) %>%
    Reduce(union, .) %>%
    length()
  return(sites)
}

# Identified phosphorylated sites
phos.sites <-
  mascot %>%
  group_by(prot_acc) %>%
  summarize(n_sites = countSites(pep_start, pep_var_mod_pos)) %>%
  summarize(PhosphoSites = sum(n_sites))

# Proportion: non-phosphorylated / phosphorylated peptides
pep.phos.ratio <-
  mascot %>%
  distinct(pep_seq, pep_var_mod) %>%
  count(Phosph = grepl("Phosp", pep_var_mod)) %>%
  summarize('Nonphosph/Phosph' = n[!Phosph]/n[Phosph])


table1 <- bind_cols(prot.num,
                    prot.phos.num,
                    prot.per.run.cv,
                    phos.sites,
                    pep.phos.ratio)

grid.newpage()
grid.draw(tableGrob(table1, rows = NULL))

# Table with median_CV for each sample_group in the experiment
median.peak.cv <-
  peakCV(df) %>%
  summarize_at(vars(-peak_ID), median, na.rm = T)

grid.newpage()
grid.draw(tableGrob(median.peak.cv, rows = NULL))


# Histograms: Missed cleavages, Phos peptide score distribution, Peptides per protein
barplot(title, table(mascot$pep_miss), space = 0)
hist(mascot$pep_score)

pep.per.prot <-
  mascot %>%
  distinct(prot_acc, pep_seq) %>%
  count(prot_acc)
hist(pep.per.prot$n, breaks = 20)

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
