#Functions setup ($Pr[1])
  anova.setup <- function(x,y) anova(lm(x ~ y))$Pr[1]
  tukey.setup <- function(x,y) TukeyHSD(aov(x ~ y))
  ##Pairwise comparison (no idea what it does)
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
 
  

#Calculate anova for groups in DF
df_anova <- function(df){
groups <- factor(arrange(df$sampleData, group)$group, 
                 levels = unique(arrange(df$sampleData, group)$group)) 
p.values <- apply(intMatrix(df),1, anova.setup, groups)
df_p.values <- p.adjust(p.values, "BH")
df_p.values <- as.data.frame(df_p.values) %>% filter(df_p.values <0.05)
df_p.values <- cbind(peakID = rownames(df_p.values), df_p.values)
rownames(df_p.values) <- NULL
return(df_p.values)
}
#Post hoc test for DF
df_tukey <- function(df){
groups <- factor(arrange(df$sampleData, group)$group, 
                 levels = unique(arrange(df$sampleData, group)$group)) 
tuk <- apply(intMatrix(df),1, tukey.setup, groups)
gp <- getPairs(groups)
    ## extract the TukeyHSD p-values
    tuk.p.values <- sapply(tuk, function(x) {
      z <- x$y[,"p adj"]
      z <- z[gp]
      names(z) <- gp
      z
      })
df.tuk <- t(tuk.p.values)
df.tuk<- as.data.frame(df.tuk)
df.tuk <- cbind(peakID = rownames(df.tuk), df.tuk)
rownames(df.tuk) <- NULL
df.tuk <- gather(df.tuk, contains("-"), key="Comparison", value="p.val") %>% 
  filter(p.val<0.05) %>% 
  arrange(p.val)

return(df.tuk)
}