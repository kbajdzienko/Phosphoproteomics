# Add grouping variables (time, treatment) and colors to sampleData
# One function for each dataset

#EXP16
make_sampleData_exp16 <- function(df) {
  sampleData <-
    df$sampleData %>%
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

  # Add colors and shapes
  timeline <- sort(unique(sampleData$time))
  colorMatrix <-
    laply(c("red", "blue"),
          function(color) colorRampPalette(c("white", color))(length(timeline)))
  rownames(colorMatrix) <- unique(sampleData$treatment)
  colnames(colorMatrix) <- sort(unique(sampleData$time))

  df$sampleData <-
    sampleData %>%
    mutate(color = colorMatrix[cbind(treatment, as.character(time))])

  return(df)
}



