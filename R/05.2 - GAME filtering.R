#General Additive Model Evaluation (GAME)

#This function will report a table with information about which part of permuted data fits the GAM model
#better than the original data. Three values are reported: residual sum of squares, GCV, deviance explained

#First the function will create fold-change data frame from generic data(df list by Dmitri)
#sample_ID has to contain treatment information
#Here the fc_mode means whether the FC should be calculated based on the treatment (treatment/control for each time point)
# or based on time (example are two treatments, and the control_group is time point 0)
#In second case(time fc_mode, it is important to set the group_to_model for the first treatment, and then repeat GAME 
#for the second treatment

#After FC table is computed, specified number of permutations of the original data is calculated
#Here its important to set the time_grouping correctly (time point 0 has to be included if its present in the data)
#It would be advised to first test the function with small (100) permutations and if its working with 1000-4000 permutations

#Then original and permuted data are combined and to each group of replicates for each protein/p-site GAM model is fitted
#and fittness parameters calculated

#Lastly fitness parameters are compared between permuted and original groups of replicates and p-value is calculated

report_game <- function(df, 
                        data_type="phospho", 
                        fc_mode="treatment", 
                        control_group = "DMSO", 
                        group_to_model = "AZD",
                        permutations= 1000, 
                        time_grouping = c(5,15,60,240)) {

  match.arg(fc_mode,c("time","treatment"))
  match.arg(data_type, c("phospho", "protein"))
            
  
#Import data
  if(data_type=="phospho") {
df2 <- fillNA_phos(df)    
data <- right_join(df2$sampleData, df2$annIntData) %>%
  rename(ID = ann_ID)
} else if(data_type == "protein"){
df2 <- fillNA(df)  
  data <- right_join(df2$sampleData, df2$intData) %>%
    rename(ID = peak_ID)
}
  
#Calculate FC 
if(fc_mode == "treatment"){
 data2 <-  mutate(data, group="original") %>%
    group_by(ID, time) %>%
    mutate(FC = log(intensity/mean(intensity[treatment==control_group]),2)) %>%
    ungroup()%>%
   filter(treatment %in% group_to_model) %>%
   select(sample_ID, ID, time, FC,group)
 
 } else if (fc_mode == "time") {
  data2 <-  mutate(data, group="original") %>%
    group_by(ID, treatment) %>%
    mutate(FC = log(intensity/mean(intensity[time==control_group]),2)) %>%
    ungroup()%>%
    filter(treatment %in% group_to_model) %>%
    select(sample_ID, ID, time, FC,group)
}
original.data<-data2
data <- data2

#Add permuted values
library(boot)
permutation_mean <- function(data, indices) {         
  mean(data[unique(indices), "FC"], na.rm = TRUE)  
}                                                     

## bootstrapping with 1000 replications - i.e. we are creating at each timepoint 1000 'replicates' of the permuted data
results <- plyr::ddply(data, plyr::.(ID), function(each){
  data.frame(sample_ID = NA, FC = boot(data = each, statistic = permutation_mean, R = permutations)$t)})


results <- group_by(arrange(results,ID),ID) %>% 
  mutate(time=rep(time_grouping,(permutations/length(time_grouping))))%>%
  ungroup()%>%
  group_by(ID,time)%>%
  mutate(group = paste("bootstrap",seq(1,(permutations/length(time_grouping))), sep="_")) %>%
  ungroup()



## Now we append the newly created data to the bottom of your data
data <- rbind(data.frame(original.data),  # I'm creating a 'group' factor to differentiate both data
              data.frame(results[,colnames(original.data)])) # Here I'm calling the 'results' data frame in the same column order as the "data"


#Create gam model for each trait (ID)

RSS_gam <- function(ddf) {
  
  mgcv::gam(FC ~ s(time, k=3),data= ddf) %>%
    summary()%>%
    unlist()%>%
    .[c(56,26,12)] %>%
    unlist()
}


pv_df <-  plyr::ddply(data, plyr::.(ID, group), RSS_gam, .progress = "text") %>%
  tbl_df() %>%
  setNames(c("ID", "group", "GCV","dev.expl","r.sq"))


pv_df2 <- pv_df %>%
  group_by(ID) %>%
  summarize(r.sq_ratio = sum(r.sq>r.sq[group=="original"])/permutations,
            GCV_ratio = sum(GCV<GCV[group=="original"])/permutations,
            dev_expl_ratio = sum(dev.expl>dev.expl[group=="original"])/permutations)

return(pv_df2)}

#mutate(data,sample_ID = stringr::str_extract(sample_ID, "\\d{2,3}[[:punct:]]\\d{1,2}$"))