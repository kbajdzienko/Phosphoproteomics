report_FC_p <- function(df, group1, group2,
                             method = "t.test",
                            data.type = "peakID",
                             FC.threshold = 1.5,
                             p.threshold = 0.05) {
  
  match.arg(method, c("t.test(var.equal=TRUE)", "wilcox.test"))
  match.arg(data.type, c("peakID","psite"))
  
  if (!all(c(group1, group2) %in% df$sampleData$group))
    stop("Selected group names are not found")
  
  if (data.type=="psite"){
  stat_mat <-
    df$sampleData %>%
    filter(group %in% c(group1, group2)) %>%
    left_join(df$annIntData) %>%
    group_by(ann_ID) %>%
    summarize(p.value = do.call(method, list(intensity[group == group1], intensity[group == group2]))$p.value,
              FC = median(intensity[group == group1], na.rm = T)/median(intensity[group == group2], na.rm = T)) %>%
    mutate(signif = abs(p.value) < p.threshold & abs(FC) > FC.threshold) %>%
    mutate(q.value = p.adjust(p.value, "BH"))
  }  
  else if(data.type=="peakID"){
  stat_mat <-  filter(df$sampleData, group %in% c(group1, group2)) %>%
  left_join(df$intData) %>%
    group_by(peak_ID) %>%
    summarize(p.value = do.call(method, list(intensity[group == group1], intensity[group == group2]))$p.value,
              FC = mean(intensity[group == group1], na.rm = T)/mean(intensity[group == group2], na.rm = T)) %>%
    mutate(signif = abs(p.value) < p.threshold & abs(FC) > FC.threshold) %>%
    mutate(q.value = p.adjust(p.value, "BH"))}
  return(stat_mat)
}
  