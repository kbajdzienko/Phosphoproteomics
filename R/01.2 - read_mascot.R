# Functions reading mascot files

read.mascot <- function(mascot_file, data = "pep") {

  match.arg(data, c("pep", "query"))
  file.in <- file(mascot_file, 'rt')
  mascot <- readLines(mascot_file, n = 100)
  line_prot <- grep("prot_hit_num", mascot)

  # Read file in chunks till Queries are found
  i <- -1
  while(!any(grepl("Queries", mascot))) {
    i <- i + 1
    mascot <- readLines(file.in, n = 100)
  }
  line_queries <- (i*100 + grep("Queries", mascot))

  if (data == "pep") {
    # Read peptides table till queries table
    skip <- line_prot - 1
    nrows <- line_queries - line_prot - 2
  } else if (data == "query") {
    # Read query table till the end
    skip <- line_queries + 1
    nrows <- -1L
  }

  close(file.in)


  mascot <- tbl_df(read.csv(mascot_file,
                            skip = skip,
                            nrows = nrows,
                            header = T, sep = ",",
                            stringsAsFactors = F,
                            colClasses = "character"))

  if (data == "pep") {
    mascot <-
      mascot %>%
      select(prot_acc, prot_score,
             pep_query, pep_isbold, pep_isunique,
             pep_start, pep_end, pep_miss,
             pep_exp_mz, pep_exp_mr, pep_delta,
             pep_score, pep_seq, pep_var_mod, pep_var_mod_pos,
             pep_scan_title) %>%
      mutate_at(vars(prot_score:pep_miss), as.integer) %>%
      mutate_at(vars(pep_exp_mz:pep_score), as.numeric) %>%
      mutate_at(vars(pep_isbold, pep_isunique), as.logical) %>%
      mutate(prot_acc = sub("\\.\\d","",prot_acc))%>%
      dplyr::rename(query_number = pep_query) %>%
      select(query_number, pep_score, pep_seq, prot_acc, everything())

  } else if (data == "query") {
    mascot <-
      mascot %>%
      filter(pep_seq != "") %>%
      select(query_number, pep_score, pep_seq,
             pep_var_mod, pep_var_mod_pos, pep_scan_title,
             pep_var_mod_conf) %>%
      mutate_at(vars(query_number), as.integer) %>%
      mutate_at(vars(pep_score), as.numeric) %>%
      mutate(pep_var_mod_conf = as.numeric(gsub("\\%","", pep_var_mod_conf)))
  }

  return(mascot)
}
