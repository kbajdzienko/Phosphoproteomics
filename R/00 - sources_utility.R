for (pckg in c("plyr", "dplyr", "tidyr"))  {
  if (!require(pckg, character.only = T)) install.packages(pckg)
}

library(plyr)
library(dplyr)
library(tidyr)

source("R/01 - tidy_PQI.R")
source("R/02 - make_annID.R")
source("R/03 - filter_unique_peptides.R")
source("R/04 - norm_transform_fillNA.R")
source("R/report_PQI.R")
source("R/report_quality.R")
source("R/report_abundance_distribution.R")
source("R/report_plot_CV_hist.R")

zero.to.na <- function(x) {
  x[x==0] <- NA
  return(x)
}
