if (!require("plyr", character.only = T)) install.packages("plyr")
if (!require("tidyverse", character.only = T)) install.packages("tidyverse")

library(plyr)
library(tidyverse)


source("R/01.1 - tidy_PQI.R")
source("R/01.2 - mascot_pep_query.R")
source("R/01.3 - intData_formats.R")
source("R/02 - make_annID.R")
source("R/03 - filter_unique_peptides.R")
source("R/04 - norm_transform_fillNA.R")
source("R/05 - filter_NA_filter_score.R")
source("R/06 - phos_sites_merging.R")
source("R/99 - report_PQI.R")
#source("R/99 - report_quality.R")
source("R/99 - report_quality_DF.R")
source("R/99 - report_abundance_distribution.R")
source("R/99 - report_plot_CV_hist.R")

zero.to.na <- function(x) {
  x[x==0] <- NA
  return(x)
}
