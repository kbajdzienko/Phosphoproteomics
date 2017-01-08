# Functions reading mascot files

read.mascot <- function(mascot_file, data = c("pep", "query")) {
  mascot <- readLines(mascot_file)
  if (data == "pep") {
    skip <- grep("prot_hit_num", mascot) - 1
    nrows <- grep("Queries", mascot) - grep("prot_hit_num", mascot) - 1
  } else if (data == "query") {
    skip <- grep("query_number", mascot) - 1
    nrows = -1L
  }

  mascot <- tbl_df(read.csv(mascot_file,
                            skip = skip,
                            nrows = nrows,
                            header = T, sep = ",",
                            stringsAsFactors = F,
                            colClasses = "character"))
  return(mascot)
}


