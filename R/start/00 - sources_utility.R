
#Setup the environment
if (!require("tidyverse", character.only = T)) install.packages("tidyverse")
source("https://bioconductor.org/biocLite.R")
biocLite("pcaMethods")

#Load all scripts in R folder
#HOME
file.sources = list.files(pattern="*.R$", ignore.case=TRUE,
      path = "C:/Users/Thinkpad/ownCloud/EXP/EXP_20_Phosphoproteomics_Data_Analysis/R_pipeline2/R_pipeline/R/",full.names = T) #Insert path of the scripts
#MPIMP
file.sources = list.files(pattern="*.R$", ignore.case=TRUE,
                          path = "W:/Krzysztof/EXP/EXP_20_Phosphoproteomics_Data_Analysis/R_pipeline2/R_pipeline/R/",full.names = T) #Insert path of the scripts




sapply(file.sources,source,.GlobalEnv)

