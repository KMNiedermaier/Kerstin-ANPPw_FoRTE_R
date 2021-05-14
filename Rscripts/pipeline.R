
# change date to proper Data format
dendro_data$date <- as.Date(dendro_data$date)

# arrange by tag and week for organization
dendro_data <- arrange(dendro_data, tag, week)

#selected only for columns that I will use 
dendro_data <- select(dendro_data, week, subplot, tag, 
                      fate, species, DBH_cm, bands_in, DOY, date, notes) %>%
  mutate(bands_cm= bands_in*2.54) %>% 
  filter(!is.na(bands_cm))

#create uniqueID for dendrodata 
dendro_data$uniqueID <- paste(dendro_data$subplot, dendro_data$week,
                              dendro_data$species, dendro_data$fate, sep = "_")


#bring in all trees tagged experiment wide (all trees within experimental subplots)
all_trees <- read.csv("data/FoRTE_all_trees.csv") %>% arrange(Tag)

#select columns
all_trees <- all_trees %>% 
  select(SubplotID, Species, Tag, fate, dbh)

# rename column names to match dendro_data
names(all_trees) <- c("subplot", "species", "tag", "fate", "DBH_cm")

# creating  "weeks" (# of data points for each tag) to prep for join with
# dendro_data; mutate in the allometric parameters for each species, 
# then join and fill in DOY.  
# THIS COULD BE CLEANER: Needs to be changed everytime new data is collected 
# need to add a and b to this df because POTR, TSCA, AMEL are not in my sampled population
# and therefore will not have a or b's when I join this table with sample pop. 

weeks <- nrow(filter(dendro_data, tag == 3))

all_trees_weeks <- all_trees %>% 
  #dplyr::setdiff(all_reps)
  slice(rep(1:n(), each = weeks)) %>% 
  group_by(tag) %>% 
  mutate(week = row_number()) %>% 
  mutate(a = case_when(
    species == "ACPE" ~ 0.03117,
    species == "ACRU" ~ 0.03177,
    species == "ACSA" ~ 0.1693,
    species == "AMEL" ~ 0.1630,
    species == "BEPA" ~ 0.0301,
    species == "FAGR" ~ 0.1892,
    species == "PIRE" ~ 0.0526,
    species == "PIST" ~ 0.0408,
    species == "POGR" ~ 0.1387,
    species == "POTR" ~ 0.0589,
    species == "QURU" ~ 0.0398,
    species == "ABBA" ~ 0.0705,
    species == "TSCA" ~ 0.1617
  )) %>% 
  mutate(b = case_when(
    species == "ACPE" ~ 2.7780,
    species == "ACRU" ~ 2.7780,
    species == "ACSA" ~ 2.3436,
    species == "AMEL" ~ 2.4940,
    species == "BEPA" ~ 2.8387,
    species == "FAGR" ~ 2.3097,
    species == "PIRE" ~ 2.5258,
    species == "PIST" ~ 2.5735,
    species == "POGR" ~ 2.3498,
    species == "POTR" ~ 2.6235,
    species == "QURU" ~ 2.7734,
    species == "ABBA" ~ 2.4970,
    species == "TSCA" ~ 2.1536
  ))

#create subplot_week_species, fate unique_ID for application of DBH - RGR adjustment
#further down the pipe
all_trees_weeks$uniqueID <- paste(all_trees_weeks$subplot, all_trees_weeks$week,
                                  all_trees_weeks$species, all_trees_weeks$fate, sep = "_")

#Rename dendro_data to df1 for simplicity. Calculate increment in circumference cm's 
# and divide by pi for radial inc.
#create new column DBH_cm_new for the first band read of all sampled trees 
# (the first record of each grou of "tag"); this new column will be populated 
# with the new DBH for each sampled tree (based on measured growth).
df1 <- dendro_data %>%   
  mutate(inc_cm = (bands_cm - lag(bands_cm, default = first(bands_cm)))/pi) %>% 
  group_by(tag) %>% 
  mutate(DBH_cm_new = case_when(
    row_number()==1 ~ DBH_cm)) 

#zero out negative growth (assuming trees can't shrink)
df1$inc_cm[df1$inc_cm < 0] <- 0

##################################################################################
# time to radially grow my sampled trees based on measured radial growth 
# further down the pipeline I will use this measured growth to calculate RGR's 
# for each species, week, subplot, and fate 

# split df1 into tags to apply the increment growth functions and loops
df1 <- data.frame(df1)
df1.list <- split(df1, df1$tag)
# Create a function that grows measured trees based on previous DBH and DBH increment
increment <- function(x) {
  for (i in 1:nrow(x)) {
    for(j in 2:nrow(x)){
      x$DBH_cm_new[j] = x$DBH_cm_new[j - 1] + (x$inc_cm[j])
    }
    return(x)
  }
}
# apply this function to the split list of tags (trees)
df1.list_grow <- lapply(df1.list, increment)
# paste lists back together as a df 
df1 <- plyr::ldply(df1.list_grow, data.frame)

#remove some columns to make df more manageable 
df1 <- select(df1, -.id, -bands_in, -bands_cm)

# remove clutter from global environment
remove(df1.list, df1.list_grow)

###############################################################################
#calculate RGR and number of days between dates (inc_days)
RGR_plot <- df1 %>% 
  group_by(tag) %>% 
  mutate(BETWEEN0=as.numeric(difftime(date,lag(date,1))), 
         inc_days=ifelse(is.na(BETWEEN0),0,BETWEEN0),
         total_inc=cumsum(as.numeric(inc_days)))%>%
  select(-BETWEEN0)  %>%
  mutate(RGR = inc_cm / lag(DBH_cm_new, default = first(DBH_cm_new))) 

# generate mean RGR for each subplot, wk, species, fate to use as RGR for unsampled trees
RGR_wk_sp_fate <- RGR_plot %>%
  group_by(tag) %>% 
  filter(row_number()!=1) %>% ungroup() %>% # this filters out the first record for each tag (RGR is NaN)
  group_by(subplot, week, date, species, fate) %>% 
  summarise(mean_RGR_sp= mean(RGR))

# generate mean RGR for each subplot, week, fate to use as RGR for unsampled trees for which
# there is no sampled species in the subplot (i.e. i did not sample oaks in A03E, but 
# there ARE oaks in A03E)
RGR_wk_fate <- RGR_plot %>%
  group_by(tag) %>% 
  filter(row_number()!=1) %>% ungroup() %>% 
  group_by(subplot, week, date, fate) %>% 
  summarise(mean_RGR_plot= mean(RGR))

# generate mean RGR for each subplot, week (excluding fate for the 7-8 trees/weeks 
# that no species, week, fate was measured)
RGR_wk <- RGR_plot %>%
  group_by(tag) %>% 
  filter(row_number()!=1) %>% ungroup() %>% 
  group_by(subplot, week, date) %>% 
  summarise(mean_RGR_wk= mean(RGR))

# join above df's with RGR_plot so that each tree has either mean species RGR, mean_plot RGR,
# or measured RGR
RGR_mean_wk_sp <- right_join(RGR_wk_sp_fate, RGR_plot)
RGR_mean <- right_join(RGR_wk_fate, RGR_mean_wk_sp)
RGR_plot_mean <- right_join(RGR_wk, RGR_mean) %>% arrange(tag, week)

# remove some clutter from global environment
remove(all_trees, dendro_2019, dendro_2020, dendro_data, df1)

##################################################################################
#This chunk takes sampled trees  and uses regression equations to measure in which
#subplot, week, fate, species DBH had a significant relationship with growth. THis 
#information is then used to adjust RGR's that were significantly (p<0.05) related to
#DBH

#select for RGR, tag, date, subplot, species, DBH
RGR_plot_regression <- RGR_plot_mean %>%
  filter(!is.na(RGR)) %>% 
  select(subplot, week, tag, species, 
         fate, DBH_cm_new, DOY, RGR) %>% 
  group_by(tag) 

#scale to NPP by subplot, week, species, AND DBH (where necessary)
#runs regression for each subplot, week, and species AND returns object with 
#coefficients; SHOULD INCLUDE FATE HERE~~!!
regressions <- RGR_plot_regression %>% 
  group_by(subplot, week, species, fate) %>% 
  do(mod = summary(lm(RGR ~ DBH_cm_new, data = .))) #%>% #-1 removes the intercept

regressions$uniqueID <- paste(regressions$subplot, regressions$week, 
                              regressions$species, regressions$fate, sep = "_")

final_output <- matrix(nrow= nrow(regressions), ncol=4, dimnames = list(c(), 
                                                                        c("uniqueID", "pvalue", 
                                                                          "slope", "intercept")))

for(i in 1:nrow(regressions)) {
  
  final_output[i, 4] <- regressions$mod[[i]][4]$coefficients[1]
  final_output[i, 3] <- regressions$mod[[i]][4]$coefficients[2]
  final_output[i, 2] <- regressions$mod[[i]][4]$coefficients[8]
  final_output[i, 1] <- regressions$uniqueID[[i]] #may need [,i] or [i,]
  
} 

#convert matrix to data frame 
final_output_df <- as.data.frame(final_output) 

#remove NA's: NEED TO FIGURE OUT WHY NA'S VS NaN'S 
final_output_df_NAomit <-   filter(final_output_df, pvalue != "NaN", !is.na(pvalue)) 

#coerce pvalue, slope, and intercept into character and then numeric
final_output_df_NAomit$pvalue <-  as.numeric(as.character(final_output_df_NAomit$pvalue))
final_output_df_NAomit$slope <-  as.numeric(as.character(final_output_df_NAomit$slope))
final_output_df_NAomit$intercept <-  as.numeric(as.character(final_output_df_NAomit$intercept))

#filter for pvalue's < 0.05; subplot, weeks, and species where DBH is having a sig effect 
sig_DBH_effect <- filter(final_output_df_NAomit, pvalue < 0.05)

#####################################################################################
#create new RGR_obs_est_df for calculated RGR (sampled trees), avg_RGR (unsampled trees), 
#adj_RGR (for unique id's where DBH had significant effect on RGR)
#first bring in ALL unsampled trees and/or missing weeks and assign a RGR using
#group_by and fill functions
# DISCLAIMER: total_inc is off on all "new bands" 
RGR_join_fill <- RGR_plot_mean %>% 
  right_join(all_trees_weeks) %>% arrange(tag, week) %>% 
  group_by(subplot, week) %>% 
  fill(c(DOY, inc_days, total_inc, date, mean_RGR_wk), .direction = "updown") %>% 
  group_by(subplot, week, fate) %>% 
  fill(c(mean_RGR_plot), .direction = "updown") %>% 
  group_by(subplot, week, species, fate) %>% 
  fill(mean_RGR_sp, .direction = "updown") %>% 
  #this is a quick fix for total_inc and DBH_cm_new; take them out and redo
  #so that newly joined records from all_trees_weeks have them
  select(-inc_days, -total_inc, -DBH_cm_new) %>% 
  group_by(tag) %>% 
  fill(notes, .direction = "down") %>% #ensures "dead" and "felled" trees are not grown
  mutate(BETWEEN0=as.numeric(difftime(date,lag(date,1))), 
         inc_days=ifelse(is.na(BETWEEN0),0,BETWEEN0),
         total_inc=cumsum(as.numeric(inc_days))) %>%
  select(-BETWEEN0) %>% 
  left_join(sig_DBH_effect) %>% 
  mutate(DBH_cm_new = case_when(
    row_number()==1 ~ DBH_cm)) 

# Create df with an RGR_obs_est column that has the appropriate RGR for every record
# Each record will now have recorded RGR, species mean RGR, plot mean RGR, OR 
# adjusted RGR (for DBH effect on RGR)
RGR_obs_est_df <-  RGR_join_fill %>% 
  mutate(RGR_obs_est = case_when(
    !is.na(RGR) ~ RGR, 
    notes == "felled" ~ 0, #manually checked all "dead"/"felled" trees and weeks
    notes == "dead" ~ 0, 
    notes == "typo" ~ 0, #sloppy, but tag 2093, wk 16 had "typo" instead of "dead"
    #species == "ACRU" & week != 1 & fate == "kill" ~ 0, #SWELLING ADJUSTMENT!!; removes ACRU production
    is.na(RGR) & !is.na(mean_RGR_sp) & is.na(slope) & 
      !is.na(mean_RGR_plot) ~ mean_RGR_sp,
    is.na(RGR) & is.na(mean_RGR_sp) & is.na(slope) & 
      !is.na(mean_RGR_plot) ~ mean_RGR_plot,
    is.na(RGR) & !is.na(slope) &!is.na(mean_RGR_sp) & 
      !is.na(mean_RGR_plot) ~ DBH_cm*slope+intercept,
    is.na(RGR) & is.na(mean_RGR_sp) & is.na(slope) & 
      is.na(mean_RGR_plot) ~ mean_RGR_wk)) 

#################################################################################
#grow the D of all trees based on appropriate RGR's 
df <- data.frame(RGR_obs_est_df)
# this should set the RGR to 0 if negative
df$RGR_obs_est[df$RGR_obs_est < 0] <- 0
#split into list of data frames
df.list <- split(df, df$tag)
df.list <- lapply(df.list, function(x){
  for (i in 1:nrow(x)) {
    for(j in 2:nrow(x)){
      x$DBH_cm_new[j] = (x$DBH_cm_new[j - 1] * x$RGR_obs_est[j]) + x$DBH_cm_new[j-1]
    }
    return(x)
  }
})
df <- plyr::ldply(df.list, data.frame)

#create week and month columns to look at different temporal scales 
df$month <- as.Date(cut(df$date, breaks = "month"))
df$week <- as.Date(cut(df$date, breaks = "week", start.on.monday = FALSE))

#remove clutter from global environment
remove(df.list, final_output,final_output_df,
       final_output_df_NAomit, regressions, RGR_obs_est_df, RGR_plot, 
       RGR_plot_regression, RGR_wk, RGR_mean_wk_sp, RGR_wk_fate, RGR_wk_sp_fate,
       sig_DBH_effect, RGR_mean, RGR_plot_mean, all_trees_weeks, RGR_join_fill)