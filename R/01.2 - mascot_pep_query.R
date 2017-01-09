# Functions reading mascot files

read.mascot <- function(mascot_file, data = "pep") {
  match.arg(data, c("pep", "query"))
  mascot <- readLines(mascot_file)

  if (data == "pep") {
    # Read peptides table till queries table or till the end if there is no query data
    skip <- grep("prot_hit_num", mascot) - 1
    if (any(grepl("Queries", mascot))) {
      nrows <- grep("Queries", mascot) - grep("prot_hit_num", mascot) - 1
    } else {
      nrows <- -1L
    }
  } else if (data == "query") {
    # Read query table till the end
    skip <- grep("query_number", mascot) - 1
    if (length(skip) == 0) {
      stop("No query data detected in mascot file")
    }
    nrows <- -1L
  }

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
      rename(query_number = pep_query)

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
