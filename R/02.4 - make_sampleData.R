# Add grouping variables (time, treatment) and colors to sampleData
# One function for each dataset

#General sample data
make_sampleData <- function(df) {
  sampleData <-
    df$sampleData %>%
    mutate(group = gsub("(?=\\b[[:digit:]]{2}$)", "0", group, perl = TRUE)) %>%
    mutate(group = gsub("(?=[[:digit:]]{3})", "", group, perl = TRUE)) %>%
    mutate(sample_ID = stringr::str_extract(sample_ID, "\\d{1,3}$")) %>%
    mutate(treatment = stringr::str_extract(group, "[[:upper:]]{3,4}")) %>%
    mutate(time = na.to.zero(as.integer(gsub("^.*-", "", group))))
    
  df$intData <- 
    df$intData %>%
    mutate(sample_ID = stringr::str_extract(sample_ID, "\\d{1,3}$"))
  
  # Add colors
  timeline <- sort(unique(sampleData$time))
  colorVector <- colorRampPalette(c("white", "green"))(length(timeline))
  names(colorVector) <- timeline
  
  df$sampleData <-
    sampleData %>%
    mutate(color = colorVector[as.character(time)])
  
  return(df)
}


#EXP16
make_sampleData_exp16 <- function(df) {
  sampleData <-
    df$sampleData %>%
    mutate(group = gsub("Starvation", "-000", group)) %>%
    # mutate(sample_ID = paste0(gsub("-", "_", group),
    #                          gsub("^.*_", "_", sample_ID))) %>%
    mutate(treatment = gsub("-[[:alnum:]]*$", "", group)) %>%
    mutate(time = na.to.zero(as.integer(gsub("\\D", "", group))))

  azd_control <-
    filter(sampleData, grepl("Control", treatment)) %>%
    mutate(treatment = "AZD/GLU")

  glu_control <-
    filter(sampleData, grepl("Control", treatment)) %>%
    mutate(treatment = "GLU")

  sampleData <-
    sampleData %>%
    filter(!grepl("Control", treatment)) %>%
    bind_rows(azd_control, glu_control)

  # Add colors
  timeline <- sort(unique(sampleData$time))
  colorMatrix <-
    laply(c("red", "blue"),
          function(color) colorRampPalette(c("white", color))(length(timeline)))
  rownames(colorMatrix) <- unique(sampleData$treatment)
  colnames(colorMatrix) <- timeline

  df$sampleData <-
    sampleData %>%
    mutate(color = colorMatrix[cbind(treatment, as.character(time))])

  return(df)
}



#EXP18

make_sampleData_exp18 <- function(df) {
  sampleData <-
    df$sampleData %>%
    mutate(group = gsub("(?=\\b[[:digit:]]{2}$)", "0", group, perl = TRUE)) %>%
    mutate(group = gsub("(?=[[:digit:]]{3})", "PF4708671-", group, perl = TRUE)) %>%
    mutate(group = gsub("C", "Control-000", group)) %>%
    # mutate(sample_ID = paste0(gsub("-", "_", group),
    #                           gsub("^.*_", "_", sample_ID))) %>%
    mutate(time = na.to.zero(as.integer(gsub("^.*-", "", group))))

  # Add colors
  timeline <- sort(unique(sampleData$time))
  colorVector <- colorRampPalette(c("white", "green"))(length(timeline))
  names(colorVector) <- timeline

  df$sampleData <-
    sampleData %>%
    mutate(color = colorVector[as.character(time)])

  return(df)
}

