# Function, calculating QC statistical metrics based on mascot output file:
# number of identified phosphorilated proteins, sites

QC_stat <- function(mascot_file) {
  mascot <- read.mascot(mascot_file, "pep")

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

  return(bind_cols(prot.num,
                   prot.phos.num,
                   prot.per.run.cv,
                   phos.sites,
                   pep.phos.ratio))
}

# Plot several quality control histograms
plot_QC_hist <- function(mascot_file, skip = 71, score = 0) {

  mascot <- read.csv(mascot_file,
                     skip = skip,
                     header = T,
                     sep = ",",
                     stringsAsFactors = F)
  mascot <- tbl_df(mascot) %>% filter(pep_score > score)
  
  par(mfrow=c(2,2),
      oma = c(0, 0, 2, 0))

  # Histogram: Missed cleavages
  hist (mascot$pep_miss,
        freq = FALSE,
        #main =" ", 
        xlab ="Missed cleaveges", 
        ylab ="Density",
        #ylim = c(1,5000),
        xlim = c(0,3),
        breaks = seq(0,3, by=1),
        xaxt='n',
        right = FALSE)
  axis(side=1, at=seq(0.5,3,1), labels=seq(0,2,1))

  # Histogram: Phos peptide score distribution
  hist(as.numeric(mascot$pep_score), freq = FALSE)

  # Histogram: Peptides per protein
  pep.per.prot <-
    mascot %>%
    distinct(prot_acc, pep_seq) %>%
    count(prot_acc)
  hist(pep.per.prot$n, breaks = 20, freq = FALSE)

  frame()
}


