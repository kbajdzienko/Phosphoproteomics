ions <- read.csv("exp19-3 peptide ions.csv",
               skip = 2,
               header = T, sep = ",",
               stringsAsFactors = F,
               colClasses = "character")

ions <- tbl_df(ions)

ions <-
  ions %>%
  select(1:2,5,8:12) %>%
  setNames(c("ion_ID", "RT", "Neutral_mass",
             "Score", "Sequence", "Modifications",
             "Accession", "All_accessions")) %>%
  mutate_at(vars(RT:Score), as.numeric)

dupl_ions <- filter(ions, grepl(";", ions$All_accessions))

df_redup <- df
df_redup$peakData <- anti_join(df$peakData, dupl_ions, by = "Sequence")

# Get the list of proteins not having a unique peptide
prot.nonuniq <- df$peakData %>% distinct(Accession)
prot.uniq <- df_redup$peakData %>% distinct(Accession)
setdiff(prot.nonuniq, prot.uniq)

###
# Test if all sets of accessions are the same for a certain sequence
diffAcc <- function(all_acc) {
  is.diff <-
    all_acc %>%
    strsplit(";") %>%
    outer(., ., Vectorize(setequal)) %>%
    all() %>%
    !.
  return(is.diff)
}

dupl_ions %>%
  group_by(Sequence) %>%
  summarize(diff_acc = diffAcc(All_accessions)) %>%
  filter(diff_acc)

###

