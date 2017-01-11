# Function, calculating QC statistical metrics based on mascot output file:
# number of identified phosphorilated proteins, sites

QC_stat_mascot <- function(mascot_file) {
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
plot_QC_hist_mascot <- function(mascot_file, score = 0) {
  
  mascot <- read.mascot(mascot_file, "pep")
  
  mascot <- tbl_df(mascot) %>% filter(pep_score > score)
  
  par(mfrow=c(2,2),
      oma = c(0, 0, 2, 0))

     # Histogram: Missed cleavages
    pep.miss <-
      mascot %>%
      distinct(pep_seq, pep_miss) %>%
      .$pep_miss
    barplot(table(pep.miss)/length(pep.miss),
            space = 0,
            main = "Missed cleavages number",
            col = "white",
            ylab = "Frequency",
            ylim = c(0, 1))
    
    # Histogram: Phosphorilation sites per peptide
    phosh.sites.counts <-
      mascot %>%
      distinct(pep_seq, pep_var_mod_pos) %>%
      .$pep_var_mod_pos %>%
      stringr::str_count("3")
    barplot(table(phosh.sites.counts)/length(phosh.sites.counts),
            space = 0,
            main = "Phosphorilation sites per peptide",
            col = "white",
            ylab = "Frequency")
    
    # Histogram: Peptides per protein
    pep.per.prot <-
      mascot %>%
      distinct(prot_acc, pep_seq) %>%
      count(prot_acc)
    barplot(table(pep.per.prot$n)/nrow(pep.per.prot),
            space = 0,
            main = "Peptides per protein",
            col = "white",
            ylab = "Frequency")
    
    # Histogram: Phos peptide score distribution
    hist(as.numeric(mascot$pep_score),
         main = "Score distribution",
         xlab = NULL,
         breaks = seq(0,max(as.numeric(mascot$pep_score)+10), by = 10))
    
    
    par(mfrow = c(1, 1))
  }
  

