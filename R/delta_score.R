# Mascot delta-score calculation

mascot <-
  read.csv("EXP19-3 Mascot for MD delta score.csv",
           skip = 72,
           header = T, sep = ",",
           stringsAsFactors = F,
           colClasses = "character")

mascot <- tbl_df(mascot)

delta_score <-
  mascot %>%
  select(pep_scan_title, pep_score, pep_seq, pep_var_mod_pos, pep_query, everything()) %>%
  filter(grepl("Phosph", pep_var_mod)) %>%
  mutate_at(vars(pep_score), as.numeric) %>%
  arrange(pep_scan_title, desc(pep_score)) %>%
  group_by(pep_scan_title) %>%
  filter(n() > 1) %>%
#  group_by(pep_scan_title, pep_seq, pep_var_mod_pos, prot_acc) %>%
#  filter(n() > 1) %>%
  summarize(delta_score = first(pep_score) - nth(pep_score, 2)) %>%
  ungroup()
