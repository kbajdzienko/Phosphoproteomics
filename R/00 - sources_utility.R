setwd("C:/Users/Thinkpad/ownCloud/EXP/EXP_20_Phosphoproteomics_Data_Analysis/R_pipeline2/R_pipeline")
if (!require("plyr", character.only = T)) install.packages("plyr")
if (!require("tidyverse", character.only = T)) install.packages("tidyverse")

library(plyr)
library(tidyverse)
library(ggplot2)


source("R/00 - data_processing_pipeline.R")
source("R/01.1a - tidy_PQI_peptides.R")
source("R/01.1b_tidy_PQI_proteins.R")
source("R/01.2 - read_mascot.R")
source("R/01.3 - intData_formats.R")
source("R/01.4 - read_ions.R")
source("R/02.1 - make_annID.R")
source("R/02.2 - make_conf.R")
source("R/02.3 - make_all_acc.R")
source("R/02.4 - make_sampleData.R")
source("R/02.5 - remove_samples.R")
source("R/03 - filter_unique_peptides.R")
source("R/04 - norm_transform_fillNA.R")
source("R/05 - filter_NA_filter_score.R")
source("R/06 - phos_sites_merging.R")
source("R/07 - join_df.R")
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
source("R/99 - plot_ratio_profiles_whiskers.R")
source("R/99 - plot_pls.R")
source("R/99 - plot_pls_cv.R")
source("R/99 - plot_volcano.R")
source("R/99 - subset_df.R")
source("R/99 - plot_profile_clust.R")
source("R/99 - anova and tukey post-hoc.R")
source("R/99 - report_FC_p.R")
source("R/99 - report_ANOVA2.R")
source("R/99 - plot_correlation.R")


# Convert zeros to NAs
zero.to.na <- function(x) {
  x[x==0] <- NA
  return(x)
}
# Convert empty cells to NAs
empty.to.na <- function(x){
  if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
  ifelse(as.character(x)!="", x, NA)
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

#Function to format numbers on y scale
fmt_scale <- function(x){format(x,nsmall = 1,scientific = T, digits=2)}
#Function to force equal numbers of breaks on the plot scale
equal_breaks <- function(n = 3, s = 0.05, ...){
  function(x){
    # rescaling
    d <- s * diff(range(x)) / (1+2*s)
    seq(min(x)+d, max(x)-d, length=n)
  }
}
