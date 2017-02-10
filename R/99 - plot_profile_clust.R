# 1. Run Dmin()
# 2. Check for number of centeres with fuzz_c_selection()
# 2. in fuzz_c_selection All values in subsequent column shoud be equal to value in c:value
# 3. Run plot_fuzz()

# Function for plotting the clusters
plot_fuzz <- function(df, centers, m = fuzz_m_estimate(df), min.memberhip = 0,
                      nrow = 4, ncol = 4, color, time.labels,
                      file = "fuzz_clust.pdf") {

  cl <- fuzz_clust(df, centers, m)

  # Add cluster number to intData
  data <- mutate(df$intData, cluster = cl$cluster[peak_ID])

  # Replace all membership below min by negative value
  cl$membership[cl$membership < min.memberhip] <- -1

  # Create a color rainbow scale and replace matrix of membership by matrix of colors
  color <- rainbow(60, start = 0.1)
  cl$color <-
    cl$membership %>%
    aaply(1, findInterval, seq(0, 1, length.out = 60)) %>%
    aaply(1, function(x) color[x])

  # List for storage of ggplots
  plots <- list()

  for (cl.idx in sort(unique(cl$cluster))) {
    # Add colour to data and leave only peaks from certain cluster
    tmp <-
      data %>%
      mutate(colour = cl$color[peak_ID, cl.idx]) %>%
      filter(cluster == cl.idx)

    plots[[cl.idx]] <-
      ggplot(data = tmp,
             aes(x = sample_ID, y = intensity, group = peak_ID, colour = colour)) +
      geom_line() +
      guides(colour = FALSE) +
      labs(title = paste("Cluster", cl.idx)) +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.background=element_rect(fill='white', colour="black"),
            panel.grid.minor = element_blank())
  }

  pdf(file)
  npages <- 1 + (length(plots) - 1)%/%(nrow*ncol)
  for (j in 1:npages) {
    do.call(gridExtra::grid.arrange,
            c(plots[((j - 1)*ncol*nrow + 1):min(j*ncol*nrow, length(plots))],
              nrow = nrow, ncol = ncol))
  }
  dev.off()
}


# Fuzzy c-means clustering
fuzz_clust <- function(df, centers, m = fuzz_m_estimate(df), ...){
  cl <- e1071::cmeans(intMatrix(df),
                      centers = centers,
                      method = "cmeans",
                      m = m,
                      ...)
  return(cl)
}

# Define parameters

# Estimate m
fuzz_m_estimate <- function(df) {
  N <- nrow(df$intData)
  D <- ncol(df$intData)
  m.sj <- 1 + (1418/N + 22.05)*D^(-2) + (12.33/N +0.243)*D^(-0.0406*log(N) - 0.1134)
  return(m.sj)
}

# Check number of empty clusters for disfferent cluster numbers
fuzz_c_selection <- function(df, m = fuzz_m_estimate(df), crange = seq(4, 32, 4),
                       repeats = 5, visu = TRUE, ...){

  Nonempty <- matrix(0, ncol = length(crange), nrow = repeats)
  i <- 0

  for (c in crange){
    i <- i + 1
    for (ii in  1:repeats){
      Utmp <- fuzz_clust(df, centers = c, m = m)$membership #,...)[[4]]
      for (jj in 1:ncol(Utmp)){
        if((sum(Utmp[ , jj] > 0.5)) > 0){
          Nonempty[ii, i] <-  Nonempty[ii, i] + 1
        }
      }
    }
  }
  dimnames(Nonempty) <- list(paste0("repeats:", c(1:repeats)),
                             paste0("c:", crange))
  if (visu) {
    plot(crange, Nonempty[1, ], pch="x")
    for (i in 1:nrow(Nonempty)) points(crange, Nonempty[i, ], pch = "x")
    lines(c(0, max(crange)), c(0, max(crange)), col = "red")
  }
  Nonempty
}


#Plot minimal distance between clusters for different cluster numbers
Dmin <- function(df, m = fuzz_m_estimate(df),
                 crange = seq(4, 40, 4), repeats = 3, visu=TRUE){
  DminM <- matrix(0, nrow=length(crange), ncol=repeats)
  for (ii in 1:repeats) {
    j <- 0
    for (i in crange){
      cl <- fuzz_clust(df, c = i, m = m)
      DminM[j <- j + 1, ii] <- min(dist(cl$centers))
    }
  }
  DminMav <- apply(DminM, 1, mean)

  if (visu) plot(crange, DminMav, xlab="Cluster number", ylab="Min. centroid distance")

  return(DminMav)
}



