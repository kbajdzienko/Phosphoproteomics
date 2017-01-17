# Function, calculating QC statistical metrics based on mascot output file:
# number of identified phosphorilated proteins, sites

QC_stat <- function(df) {

  # All identified proteins
  prot.num <- summarize(df$peakData, Proteins = n_distinct(Accession))

  # Phosphorylated proteins
  prot.phos.num <-
    filter(df$peakData, grepl("Phosp", Modifications)) %>%
    summarize('(P)-proteins' = n_distinct(Accession)) %>%
    '/'(prot.num/100) %>%
    mutate_all(funs(paste0(round(., digits = 1), "%")))

  # All identified peptides
  pep.num <- summarize(df$peakData, Peptides = n_distinct(Sequence, Modifications))

  # Phosphorylated peptides
  pep.phos.num <-
    filter(df$peakData, grepl("Phosp", Modifications)) %>%
    summarize('(P)-peptides' = n_distinct(Sequence, Modifications)) %>%
    '/'(pep.num/100) %>%
    mutate_all(funs(paste0(round(., digits = 1), "%")))

  # Identified phosphorylated sites
  phos.sites <-
    sitesMerge(df)$annData %>%
    group_by(Accession) %>%
    summarize(n_sites = n()) %>%
    summarize('(P)-Sites' = sum(n_sites), 'Sites/protein' = sum(n_sites)/n())

  return(bind_cols(prot.num,
                   prot.phos.num,
                   pep.num,
                   pep.phos.num,
                   phos.sites))
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
          ylab = "Frequency",
          ylim = c(0,1)
          )

  # Histogram: Phosphorilation sites per peptide
  phosh.sites.counts <-
    df$peakData %>%
    distinct(Sequence, Modifications) %>%
    .$Modifications %>%
    stringr::str_count("Phospho")
  barplot(table(phosh.sites.counts)/length(phosh.sites.counts),
          space = 0,
          main = "Phosphorylation sites per peptide",
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
       xlab = NULL,
       breaks = seq(0,max(as.numeric(df$peakData$Score)+10), by = 10))

  par(mfrow = c(1, 1))
}
