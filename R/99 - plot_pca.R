# Plot PCA loadings plot.
# inx1, inx2 - the order numbers of PCs.
# showNames - if TRUE then points at the plot are labeled.

plot_PCA_loadings <- function(df, inx1 = 1, inx2 = 2, showNames = F){
  pca <- PCA_anal(df)
  loadings <- pca$rotation[, c(inx1, inx2)]
  colnames(loadings) <- paste("Loadings", c(inx1, inx2))
  par(mar = c(6, 5, 2, 6))
  plot(loadings[, 1], loadings[, 2], las = 2,
       xlab = colnames(loadings)[1], ylab = colnames(loadings)[2]);
  grid(col = "lightgray", lty = "dotted", lwd = 1)
  points(loadings[, 1], loadings[, 2], pch = 19, col = "red");
  if(showNames){
    text(loadings[, 1], loadings[, 2],
         labels = substr(rownames(loadings), 1, 12),
         pos = 4, col = "blue", xpd = T);
  }
}

# Plot 2D score plot.
# inx1,inx2 - the order numbers of PCs.
# reg - set the confidence level for plotting confidence region ellipse.
# showNames - if TRUE then points at the plot are labeled.

plot_PCA_scores <- function(df, inx1 = 1, inx2 = 2,
                            reg = 0.95,
                            showNames = TRUE,
                            color = NULL) {
  pca <- PCA_anal(df)
  xlabel = paste("PC", inx1, "(", round(100*pca$variance[inx1], 1), "%)");
  ylabel = paste("PC", inx2, "(", round(100*pca$variance[inx2], 1), "%)");
  pc1 = pca$x[, inx1];
  pc2 = pca$x[, inx2];
  text.lbls <- substr(names(pc1), 1, 14) # some names may be too long

  op <- par(mar = c(5,5,3,3));

  color <-  df$sampleData$color
  plot(pc1, pc2, xlab = xlabel, ylab=ylabel, type='n', main="Scores Plot")
  points(pc1, pc2, pch = 23, col = "black", bg = color, cex = 2)
  text(pc1, pc2, label = text.lbls, pos=4, col ="blue", xpd=T, cex=0.8)

  par(op)
}



PCA_anal <-function(df){
  pca <-
    df %>%
    filter_NA() %>%
    fillNA() %>%
    intMatrix() %>%
    t() %>%
    prcomp(center = T, scale = T)

  # obtain variance explained
  sum.pca <- summary(pca);
  imp.pca <- sum.pca$importance;
  std.pca <- imp.pca[1,]; # standard deviation
  var.pca <- imp.pca[2,]; # variance explained by each PC
  cum.pca <- imp.pca[3,]; # cummulated variance explained

  # store the item to the pca object
  pca <- append(pca, list(std = std.pca,
                          variance = var.pca,
                          cum.var = cum.pca))
  return(pca);
}



#
#
#   if (dataSet$cls.type == "disc"){
#     # obtain ellipse points to the scatter plot for each category
#     lvs <- unique(df$sampleData$group)
#     pts.array <- array(0, dim = c(100, 2, length(lvs)),
#                        dimnames = list(NULL, NULL, lvs))
#
#     for(lvl in lvs) {
#       smpls <-
#         df$sampleData %>%
#         filter(group == lvl) %>%
#         .$sample_ID
#       groupVar <- var(cbind(pc1[smpls], pc2[smpls]), na.rm = T)
#       groupMean <- cbind(mean(pc1[smpls], na.rm = T), mean(pc2[smpls], na.rm = T))
#       pts.array[ , , lvl] <- ellipse::ellipse(groupVar,
#                                          centre = groupMean,
#                                          level = reg)
#     }
#
#     # Define limits and extend them a bit for plotting
#     xrg <- range (pc1, pts.array[, 1, ]);
#     yrg <- range (pc2, pts.array[, 2, ]);
#     xlims <- xrg + c(-0.08, 0.08)*(xrg[2] - xrg[1])
#     ylims <- yrg + c(-0.08, 0.08)*(yrg[2] - yrg[1])
#
#     # cols <- GetColorSchema(dataSet, gray.scale);
#     # uniq.cols <- unique(cols);
#
#     plot(pc1, pc2,
#          xlab = xlabel, xlim = xlims, ylim = ylims, ylab = ylabel,
#          type  ='n', main = "Scores Plot",
#          #col = cols,
#          pch = as.numeric(as.factor(df$sampleData$group)) + 1); ## added
#     grid(col = "lightgray", lty = "dotted", lwd = 1);
#
#     # make sure name and number of the same order DO NOT USE levels, which may be different
#     legend.nm <- unique(df$sampleData$group);
#     ## uniq.cols <- unique(cols);
#
#     # ## BHAN: when same color is choosen; it makes an error
#     # if ( length(uniq.cols) > 1 ) {
#     #   names(uniq.cols) <- legend.nm;
#     # }
#
#     # draw ellipse
#     for(lvl in lvs){
#       polygon(pts.array[,,lvl], col = adjustcolor("black", alpha.f = 0.25), border = NA);
#     }
#
#       if (length(uniq.cols) > 1) {
#         polygon(pts.array[,,i], col=adjustcolor(uniq.cols[lvs[i]], alpha.f=0.25), border=NA);
#       } else {
#
#       }
#       if(gray.scale) {
#         lines(pts.array[,,i], col=adjustcolor("black", alpha.f=0.5), lty=2);
#       }
#     }
#
#     pchs <- GetShapeSchema(dataSet, show, gray.scale);
#     if(gray.scale) {
#       cols <- rep("black", length(cols));
#     }
#     if(show){
#       text(pc1, pc2, label=text.lbls, pos=4, xpd=T, cex=0.75);
#       points(pc1, pc2, pch=pchs, col=cols);
#     }else{
#       if(length(uniq.cols) == 1){
#         points(pc1, pc2, pch=pchs, col=cols, cex=1.0);
#       }else{
#         if(gray.scale | (!is.null(dataSet$shapeVec) && all(dataSet$shapeVec>0))){
#           points(pc1, pc2, pch=pchs, col=cols, cex=1.8);
#         }else{
#           points(pc1, pc2, pch=21, bg=cols, cex=2);
#         }
#       }
#     }
#     uniq.pchs <- unique(pchs);
#
#     legend("topright", legend = legend.nm, pch=uniq.pchs, col=uniq.cols);
#
#   par(op)
# }

