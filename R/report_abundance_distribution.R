library(ggplot2)

text.size <- 4
text.angle <- 0
x.axis.size=10
y.axis.size=10
yaxis.name <- 'Log2-intensities'
y.limup <- ceiling(max(log2(df$intData$intensity), na.rm=TRUE) + 3)
y.limdown <- -1

tempGroupName <- df$sampleData # unique(datafeature[, c("GROUP_ORIGINAL", "RUN")])
groupAxis <- as.numeric(xtabs(~group, tempGroupName))
cumGroupAxis <- cumsum(groupAxis)
lineNameAxis <- cumGroupAxis[-nlevels(as.factor(df$sampleData$group))]

groupName <- data.frame(sample_ID = c(0, lineNameAxis)+groupAxis / 2 + 0.5,
                        intensity = rep(y.limup - 1, length(groupAxis)),
                        Name = levels(as.factor(df$sampleData$group)))

label.color <- c("lightblue")
df$intData <-
  df$intData %>%
  mutate(LABEL = factor("Endogenous")) %>%
  mutate(intensity = log2(intensity))

ptemp <-
  ggplot(aes_string(x = 'sample_ID', y='intensity'), data = df$intData) +
  facet_grid(~LABEL) +
  geom_boxplot(aes_string(fill = 'LABEL'), outlier.shape = 1, outlier.size = 1.5) +
  scale_fill_manual(values = label.color, guide = "none")+
  scale_x_discrete('MS runs', breaks = cumGroupAxis)+
  scale_y_continuous(yaxis.name, limits = c(y.limdown, y.limup))+
  geom_vline(xintercept = lineNameAxis + 0.5, colour = "grey", linetype = "longdash") +
  #labs(title = "All proteins")+
  geom_text(data = groupName, aes(x = sample_ID, y = intensity, label = Name),
            size = text.size, angle = text.angle, color = "black")+
  theme(
    panel.background=element_rect(fill='white', colour="black"),
    legend.key=element_rect(fill='white', colour='white'),
    panel.grid.minor = element_blank(),
    strip.background=element_rect(fill='gray95'),
    strip.text.x=element_text(colour=c("#00B0F6"), size=14),
    axis.text.x=element_text(size=x.axis.size,colour="black"),
    axis.text.y=element_text(size=y.axis.size,colour="black"),
    axis.ticks=element_line(colour="black"),
    axis.title.x=element_text(size=x.axis.size+5, vjust=-0.4),
    axis.title.y=element_text(size=y.axis.size+5, vjust=0.3),
    title=element_text(size=x.axis.size+8, vjust=1.5)
  )

print(ptemp)
