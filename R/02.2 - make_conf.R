make_conf <- function(df, mascot_file) {

  mascot <- read.mascot(mascot_file, "pep")
  query <- read.mascot(mascot_file, "query")

  # To queries with only one annotation assign confidence score 100%
  query$pep_var_mod_conf[!dupl(query$query_number)] <- 100

  # Join pep_var_mod_conf from query table to peptide table
  mascot <-
    left_join(mascot, query) %>%
    rename(Accession = prot_acc,
           Sequence = pep_seq,
           Score = pep_score,
           Neutral_mass = pep_exp_mr) %>%
    distinct()

  # Create the future dataframe for conf scores
  df_conf <- NULL
  # Copy objects to change in cycle
  peaks. <- df$peakData
  mascot. <- mascot

  # To get match between df and mascot, round Neutral_mass to different extent
  for (digits in 6:0) {
    mascot.round <- mutate(mascot., Neutral_mass = round(Neutral_mass, digits = digits))
    peaks.round <- mutate(peaks., Neutral_mass = round(Neutral_mass, digits = digits))

    # Tests of matching
    mismatch <- nrow(anti_join(peaks.round, mascot.round))
    extramatch <- nrow(left_join(peaks.round, mascot.round)) - nrow(peaks.)
    message("Digits = ", digits)
    message("Number of non-matches = ", mismatch)
    message("Number of potential extra-matches = ", extramatch)

    # Get all unique matches
    matched <-
      inner_join(peaks.round, mascot.round) %>%
      filter(!dupl(peak_ID))

    # Update table of peak_ID - conf score
    df_conf <-
      select(matched, peak_ID, pep_var_mod_conf) %>%
      bind_rows(df_conf)

    # Filter peaks and mascot by already matched peaks
    peaks. <- anti_join(peaks., matched, by = "peak_ID")
    mascot. <- anti_join(mascot., matched, by = c("query_number", "Sequence",
                                                  "pep_var_mod_pos", "Accession"))

    # Add the cases when it's a multiple match, but conf score is the same
    conf_dupl <-
      inner_join(peaks.round, mascot.round) %>%
      filter(dupl(peak_ID)) %>%
      group_by(peak_ID) %>%
      # Keep peaks with the same confidence score
      filter(isTRUE(all.equal(pep_var_mod_conf, rep(first(pep_var_mod_conf), n())))) %>%
      ungroup() %>%
      distinct(peak_ID, pep_var_mod_conf)

    # Add those peaks to result as well and exclude from peaks table
    df_conf <- bind_rows(df_conf, conf_dupl)
    peaks. <- anti_join(peaks., conf_dupl, by = "peak_ID")
  }

  # TEST: multiple matching peaks left
  extramatched <-
    inner_join(peaks.round, mascot.round) %>%
    filter(dupl(peak_ID)) %>%
    distinct(peak_ID) %>%
    anti_join(conf_dupl)
  if (nrow(extramatched) != 0) {
    extramatched %>%
      unlist() %>%
      paste(collapse = "\n") %>%
      message("Peaks with multiple mascot matches:\n", .)
  } else message("No peaks with multiple mascot matches.")

  # Final output
  df$peakData <- left_join(df$peakData, df_conf)
  return(df)
}




