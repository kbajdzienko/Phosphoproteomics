# Function, calculating QC statistical metrics based on mascot output file:
# number of identified phosphorilated proteins, sites

QC_stat <- function(df) {

  # All identified proteins
  prot.num <- summarize(df$peakData, Proteins = n_distinct(Accession))

  # Phosphorylated proteins
  prot.phos.num <-
    filter(df$peakData, grepl("Phosp", Modifications)) %>%
    summarize(PhosphoProteins = n_distinct(Accession))

  # # Identified unique proteins within each run -- CV
  # prot.per.run.cv <-
  #   df %>%
  #   mutate(sample_ID = sapply(strsplit(pep_scan_title, "\\."), '[', 1)) %>%
  #   distinct(sample_ID, prot_acc) %>%
  #   count(sample_ID) %>%
  #   summarize(UniqueProteinsCV = sd(n)/mean(n))

  # Identified phosphorylated sites
  if (is.null(df$annIntData)) {
    phos.sites <-
      sitesMerge(df)$annIntData %>%
      summarize(Proteins = n_distinct(ann_ID))
  } else {
    phos.sites <-
      df$annIntData %>%
      summarize(PhosphoSites = n_distinct(ann_ID))
  }

  # Proportion: non-phosphorylated / phosphorylated peptidestable(df$peakData$pep_miss)
  pep.phos.ratio <-
    df$peakData %>%
    distinct(Sequence, Modifications) %>%
    summarize('Nonphosph/Phosph' = sum(!grepl("Phosph", Modifications))/sum(grepl("Phosph", Modifications)))

  return(bind_cols(prot.num,
                   prot.phos.num,
                   #prot.per.run.cv,
                   phos.sites,
                   pep.phos.ratio))
}

# Plot several quality control histograms
plot_QC_hist <- function(df) {

  par(mfrow=c(2,2))

  # Histogram: Missed cleavages
  pep.miss <-
    df$peakData %>%
    distinct(Sequence, pep_miss) %>%
    .$pep_miss
  barplot(table(pep.miss)/length(pep.miss),
          space = 0,
          main = "Missed cleavages number",
          col = "white",
          ylab = "Frequency")

  # Histogram: Phosphorilation sites per peptide
  phosh.sites.counts <-
    df$peakData %>%
    distinct(Sequence, Modifications) %>%
    .$Modifications %>%
    stringr::str_count("Phospho")
  barplot(table(phosh.sites.counts)/length(phosh.sites.counts),
          space = 0,
          main = "Phosphorilation sites per peptide",
          col = "white",
          ylab = "Frequency")

  # Histogram: Peptides per protein
  pep.per.prot <-
    df$peakData %>%
    distinct(Accession, Sequence) %>%
    count(Accession)
  barplot(table(pep.per.prot$n)/nrow(pep.per.prot),
          space = 0,
          main = "Peptides per protein",
          col = "white",
          ylab = "Frequency")

  # Histogram: Phos peptide score distribution
  hist(as.numeric(df$peakData$Score),
       main = "Score distribution",
       xlab = NULL)

  par(mfrow = c(1, 1))
}
