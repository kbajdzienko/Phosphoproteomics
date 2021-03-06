plot_correlation <- function(df){
  
  get_upper_tri <- function(cor_mat){
    cor_mat[upper.tri(cor_mat)]<- NA
    return(cor_mat)
  }
  
  df1_mat <- logTransform(df)%>% intMatrix()
  
  sID <- df$sampleData %>% 
    select(sample_ID2 = sample_ID, group2 = group) 
  
  #Calculate correlation matrix
  cor_mat <- cor(df1_mat, use="pairwise.complete.obs", method="pearson") %>%
    get_upper_tri() %>%
    reshape2::melt(na.rm=TRUE)%>% 
    mutate(sample_ID = as.character(Var1), sample_ID2 = as.character(Var2)) %>%
    left_join(df$sampleData) %>%
    select(sample_ID, sample_ID2, group, value)%>%
    left_join(sID)%>%
    #mutate(sample_ID = paste(group,sample_ID, sep="_"),sample_ID2 = paste(group2,sample_ID2, sep="_"))%>%
      arrange(group, group2)
  #cor_mat$sample_ID <- reorder(cor_mat$sample_ID, order(cor_mat$group))
  #cor_mat$sample_ID2 <- reorder(cor_mat$sample_ID2, order(cor_mat$group2))
  
  ggplot(data = cor_mat, aes(x=sample_ID, y=sample_ID2,fill = value))+
    geom_tile(color = "white")+
    #geom_text(aes(sample_ID,sample_ID2, label = round(value, digits = 2)), color = "black", size = 3) +
    scale_fill_gradient2(low = "skyblue3",mid = "white" ,high = "red", 
                         midpoint = 0,
                        limit = c(-1,1), space = "Lab", 
                        name="Pearson\nCorrelation") +
    scale_x_discrete(labels = distinct(cor_mat,sample_ID, .keep_all=T)$group) +
    scale_y_discrete(labels = distinct(cor_mat,sample_ID2, .keep_all=T)$group2, position = "right") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 10, hjust = 1))+
    coord_fixed() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(0, 1),
      #legend.position = c(0.6, 0.7),
      legend.direction = "vertical")+
    guides(fill = guide_colorbar(barwidth = 1, barheight = 7,
                                 title.position = "bottom", title.hjust = 0.5))
}

plot_correlation_ann <- function(df){
  
  get_upper_tri <- function(cor_mat){
    cor_mat[upper.tri(cor_mat)]<- NA
    return(cor_mat)
  }
  
  df1_mat <- logTransform_ann(df)%>% intMatrix_ann()
  
  sID <- df$sampleData %>% 
    select(sample_ID2 = sample_ID, group2 = group) 
  
  #Calculate correlation matrix
  cor_mat <- cor(df1_mat, use="pairwise.complete.obs", method="pearson") %>%
    get_upper_tri() %>%
    reshape2::melt(na.rm=TRUE)%>% 
    mutate(sample_ID = as.character(Var1), sample_ID2 = as.character(Var2)) %>%
    left_join(df$sampleData) %>%
    select(sample_ID, sample_ID2, group, value)%>%
    left_join(sID)%>%
    #mutate(sample_ID = paste(group,sample_ID, sep="_"),sample_ID2 = paste(group2,sample_ID2, sep="_"))%>%
    arrange(group, group2)
  #cor_mat$sample_ID <- reorder(cor_mat$sample_ID, order(cor_mat$group))
  #cor_mat$sample_ID2 <- reorder(cor_mat$sample_ID2, order(cor_mat$group2))
  
  ggplot(data = cor_mat, aes(x=sample_ID, y=sample_ID2,fill = value))+
    geom_tile(color = "white")+
    #geom_text(aes(sample_ID,sample_ID2, label = round(value, digits = 2)), color = "black", size = 3) +
    scale_fill_gradient2(low = "skyblue3",mid = "white" ,high = "red", 
                         midpoint = 0,
                         limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    scale_x_discrete(labels = distinct(cor_mat,sample_ID, .keep_all=T)$group) +
    scale_y_discrete(labels = distinct(cor_mat,sample_ID2, .keep_all=T)$group2, position = "right") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 10, hjust = 1))+
    coord_fixed() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(0, 1),
      #legend.position = c(0.6, 0.7),
      legend.direction = "vertical")+
    guides(fill = guide_colorbar(barwidth = 1, barheight = 7,
                                 title.position = "bottom", title.hjust = 0.5))
}
