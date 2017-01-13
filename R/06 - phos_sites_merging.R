# Split not cleaved peptides into single phosphorilation sites and sum intensitites
# for the same ann_ID (protein/site)
# Save intensity data in annIntData table
# Save annotation data in annData table

sitesMerge <- function(df) {
  # Function splitting ann_ID into by phosphorilation sites and replicating row as
  # many times as there are sites
  splitSites<- function (ddf) {
    # Split ann_ID
    ann <- unlist(strsplit(ddf$ann_ID, "_"))
    # ann[1] is an accession, ann[-1] is a vector of sites
    # paste them this way to couple accession with each site
    ann_ID_new <- paste(ann[1], ann[-1], sep = "_")
    # Replicate the raw as many times as necessary and change ann_ID
    ddf <-
      ddf[rep(1, length(ann_ID_new)), ] %>%
      mutate(ann_ID = ann_ID_new)

    #paste0(ddf[1,1], seq(length(ann_ID)), sep = "_")
    return(ddf)
  }

  # Get peptides with multiple sites from peakData and split them into single ones
  # Add this table to peakData for single-site-phosphorilated peptides
  mult_sites <-
    df$peakData %>%
    filter(grepl("_[[:alnum:]]+_", ann_ID)) %>%
    ddply(.var = "peak_ID", splitSites) %>%
    tbl_df() %>%
    bind_rows(filter(df$peakData, !grepl("_[[:alnum:]]+_", ann_ID)))

  # Get site data
  df$annData <- distinct(mult_sites, ann_ID, Accession, Description)

  # Sum intensities for the same ann_ID (same protein/phosphorilation site)
  df$annIntData <-
    df$intData %>%
    left_join(distinct(mult_sites, peak_ID, ann_ID)) %>%
    group_by(ann_ID, sample_ID) %>%
    summarize_at(vars(intensity), sum, na.rm = T) %>%
    ungroup() %>%
    mutate(intensity = zero.to.na(intensity))

  return(df)
}


# # Sum intensities for the same protein, the same modification
# sumANN <- function(df) {
#   df$annIntData <-
#     df$intData %>%
#     left_join(select(df$peakData, peak_ID, ann_ID)) %>%
#     group_by(ann_ID, sample_ID) %>%
#     summarize_at(vars(intensity), sum, na.rm = T) %>%
#     mutate(intensity = zero.to.na(intensity)) %>%
#     ungroup()
#   return(df)
# }

