plot_df_prot_boxplot <- function(df) {
ggplot(left_join(df$intData,df$sampleData), aes(x=sample_ID, y=intensity))+
  geom_boxplot(na.rm = F)+
  labs(y="log10 intensity", x="group")+
  coord_flip()+
  scale_y_log10()+
  facet_grid(group~., switch = "y", scales = "free_y", space = "free_y")+
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        axis.text.y= element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.y = element_text(angle=180)) 
}
