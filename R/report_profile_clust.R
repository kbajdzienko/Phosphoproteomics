# Fuzzy c-means clustering

# Define parameters

# Don't load Mfuzz library: masking objects from dplyr!


df <-
  tidy_PQI("DATA/exp16_peptides_phos.csv") %>%
  #make_annID("DATA/exp16_mascot.csv") %>%
  filter_NA() %>%
  fillNA("ppca")

mfuzz <- function(eset,centers,m,...){
  cl <- e1071::cmeans(intMatrix(df),
                      centers = centers,
                      method = "cmeans",
                      m = m,
                      ...)
  return(cl)
}

cl <- mfuzz(df, 16, 1.15)

mfuzz.plot <- function(eset, cl, mfrow=c(1,1),color,min.mem = 0,time.labels,new.window=TRUE){
  # function for plotting the clusters
  data <-
    df$intData %>%
    mutate(cluster = cl$cluster[peak_ID])

  color <- rainbow(60, start = 0.1)

  colMatrix <-
    cl$membership %>%
    aaply(1, findInterval, seq(0, 1, length.out = 60)) %>%
    aaply(1, function(x) color[x])

  plots <- list()
  for (cl.idx in sort(unique(cl$cluster))) {
    tmp <- filter(data, cluster == cl.idx)
    plots[cl.idx] <- ggplot(data = tmp,
                            aes(x = sample_ID,
                                y = intensity,
                                group = peak_ID)) +
    geom_line(colour = colMatrix[rownames(colMatrix) %in% tmp$peak_ID, cl.idx])
    }
  plots[i] <- ggplot(data = data, aes(x = sample_ID,
                          y = intensity,
                          group = peak_ID)) +
           geom_line() +
    scales::


  clusterindex <- cl$cluster[peak_ID]
  memship <- cl[[4]]
  memship[memship < min.mem] <- -1
  colorindex <- integer(dim(exprs(eset))[[1]])
  if (missing(colo)){
    colo <- c("#FF8F00",
              "#FFA700", "#FFBF00", "#FFD700", "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00",
              "#AFFF00", "#97FF00", "#80FF00", "#68FF00", "#50FF00", "#38FF00", "#20FF00",
              "#08FF00", "#00FF10", "#00FF28", "#00FF40", "#00FF58", "#00FF70", "#00FF87",
              "#00FF9F", "#00FFB7", "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF", "#00CFFF",
              "#00B7FF", "#009FFF", "#0087FF", "#0070FF", "#0058FF", "#0040FF", "#0028FF",
              "#0010FF", "#0800FF", "#2000FF", "#3800FF", "#5000FF", "#6800FF", "#8000FF",
              "#9700FF", "#AF00FF", "#C700FF", "#DF00FF", "#F700FF", "#FF00EF", "#FF00D7",
              "#FF00BF", "#FF00A7", "#FF008F", "#FF0078", "#FF0060", "#FF0048", "#FF0030",
              "#FF0018")

  }

  colorseq <- seq(0,1,length=length(colo))


  for (j in 1:max(clusterindex)){
    tmp <- exprs(eset)[clusterindex==j, , drop=FALSE]# thanks Ian for the fix
    tmpmem <- memship[clusterindex==j,j]

    if (((j-1)%% (mfrow[1] * mfrow[2]))==0){

      if (new.window) X11()
      par(mfrow=mfrow)

      if (sum(clusterindex==j)==0) {
        ymin <- -1; ymax <- +1;
      } else {
        ymin <- min(tmp);ymax <- max(tmp);
      }

      plot.default(x=NA,xlim=c(1,dim(exprs(eset))[[2]]), ylim= c(ymin,ymax),
                   xlab="Time",ylab="Expression changes",main=paste("Cluster",j),axes=FALSE)
      if (missing(time.labels)){
        axis(1, 1:dim(exprs(eset))[[2]],c(1:dim(exprs(eset))[[2]]))
        axis(2)
      } else {
        axis(1, 1:dim(exprs(eset))[[2]],time.labels)
        axis(2)
      }


    } else {

      if (sum(clusterindex==j)==0) {
        ymin <- -1; ymax <- +1;
      } else {
        ymin <- min(tmp);ymax <- max(tmp);
      }


      plot.default(x=NA,xlim=c(1,dim(exprs(eset))[[2]]), ylim= c(ymin,ymax),
                   xlab="Time",ylab="Expression changes",main=paste("Cluster",j),axes=FALSE)

      if (missing(time.labels)){
        axis(1, 1:dim(exprs(eset))[[2]],c(1:dim(exprs(eset))[[2]]))
        axis(2)
      } else {
        axis(1, 1:dim(exprs(eset))[[2]],time.labels)
        axis(2)
      }


    }


    if (!(sum(clusterindex==j)==0)){
      for (jj in 1:(length(colorseq)-1)){
        tmpcol <- (tmpmem >= colorseq[jj] & tmpmem <= colorseq[jj+1])
        if (sum(tmpcol)> 0) {
          tmpind <- which(tmpcol)
          for (k in 1:length(tmpind)){
            lines(tmp[tmpind[k],],col=colo[jj])
          }
        }
      }}
  }
}
