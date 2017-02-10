
if (!require("plyr", character.only = T)) install.packages("plyr")
if (!require("tidyverse", character.only = T)) install.packages("tidyverse")

library(plyr)
library(tidyverse)
library(RColorBrewer)
library(treemap)
library(ggplot2)
library(e1071)


source("R/00 - data_processing_pipeline.R")
source("R/01.1 - tidy_PQI.R")
source("R/01.2 - read_mascot.R")
source("R/01.3 - intData_formats.R")
source("R/01.4 - read_ions.R")
source("R/02.1 - make_annID.R")
source("R/02.2 - make_conf.R")
source("R/02.3 - make_all_acc.R")
source("R/02.4 - make_sampleData.R")
source("R/03 - filter_unique_peptides.R")
source("R/04 - norm_transform_fillNA.R")
source("R/05 - filter_NA_filter_score.R")
source("R/06 - phos_sites_merging.R")
source("R/99 - report_PQI.R")
source("R/99 - report_plot_QC_mascot.R")
source("R/99 - report_plot_QC_df.R")
source("R/99 - plot_CV_hist.R")
source("R/99 - plot_abundance_distribution.R")
source("R/99 - plot_protein_profile.R")
source("R/99 - plot_mass_accuracy.R")
source("R/99 - plot_norm_summary.R")
source("R/99 - report_plot_2samples.R")
source("R/99 - report_protein_profile_sites.R")
source("R/99 - plot_pca.R")
source("R/99 - plot_ratio_profiles.R")
source("R/99 - plot_pls.R")
source("R/99 - plot_pls_cv.R")
source("R/99 - plot_volcano.R")
source("R/99 - subset_df.R")
source("R/99 - plot_profile_clust.R")
source("R/99 - anova and tukey post-hoc.R")
source("R/99 - report_FC_p.R")



# Convert zeros to NAs
zero.to.na <- function(x) {
  x[x==0] <- NA
  return(x)
}
# Convert NAs to zeros
na.to.zero <- function(x) {
  x[is.na(x)] <- 0
  return(x)
}

# Function returning TRUE for elements repeating more than once
dupl <- function(x) duplicated(x) | duplicated(x, fromLast = T)

# Check if data are log-transformed
# (should work for MS data, for bases >2 for sure)
is.log <- function(intensity) all(intensity < 1000, na.rm = T)
