#Functions setup
  anova.setup <- function(x,y) anova(lm(x ~ y))
  tukey.setup <- function(x,y) TukeyHSD(aov(x ~ y))
  ##Pairwise comparison
  getPairs <- function(g) {
  z <- levels(g)
  out <- character(length(z)*(length(z)-1)/2)
  k <- 1
  for(i in 1:(length(z)-1))
    for(j in (i+1):length(z)) {
      out[k] <- paste(z[j],z[i],sep='-')
      k <- k+1
    }
  out
  }
  ## extract the TukeyHSD p-values
  tuk.p.values <- sapply(tuk, function(x) {
  z <- x$y[,"p adj"]
  z <- z[gp]
  names(z) <- gp
  z
  })


df_anova <- function(df){
groups <- factor(arrange(df$sampleData, group)$group, 
                 levels = unique(arrange(df$sampleData, group)$group)) 
p.values <- apply(intMatrix(df),1, anova.2, groups)
df_p.values <- p.adjust(p.values, "BH")
return(as.data.frame(df_p.values))
}

df_tukey <- function(df){
groups <- factor(arrange(df$sampleData, group)$group, 
                 levels = unique(arrange(df$sampleData, group)$group)) 
tuk <- apply(intMatrix(df),1, tukey.2, groups)
gp <- getPairs(groups)
tuk.p.values <- t(tuk.p.values)
return(tuk.p.values)
}