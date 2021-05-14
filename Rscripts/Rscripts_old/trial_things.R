#set library to dyplr
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotrix)

# use source to run subcanopy scripts 
source("Rscripts/1-8cm_increment&graph.R")
source("Rscripts/0-1cm_seedlings.R")

#importing csv
dendro_data <- read.csv("data/FoRTE_dendro_data_clean.csv")

all_reps_all_trees <- read.csv("data/allreps_dist_treatment_output_update.csv")

#selected only for columns that I will use 
all_reps <- select(dendro_data, week, subplot, tag, 
                   fate, species, DBH_cm, bands_in, DOY, date, notes) %>%
  mutate(bands_cm= bands_in*2.54) %>% 
  filter(!is.na(bands_cm))

#select my good good columns
all_trees <- all_reps_all_trees %>% 
  select(SubplotID, Species, Tag, fate, dbh)

#rename column names to match all_reps
names(all_trees) <- c("subplot", "species", "tag", "fate", "DBH_cm")

#creating 16 of each tag to build a week column and prep for join; then join and fill in DOY 
all_trees_weeks <- all_trees %>% 
  slice(rep(1:n(), each = 17)) %>% 
  group_by(tag) %>% 
  mutate(week = row_number()) %>% 
  left_join(all_reps) %>% 
  group_by(subplot, week) %>% 
  fill(DOY) %>%
  fill(DOY, .direction = c("up")) %>% 
  ungroup()

all_trees_weeks$uniqueID <- paste(all_trees_weeks$subplot, all_trees_weeks$week,
                                  all_trees_weeks$species, sep = "_")
#building df with all trees sampled and unsampled
#moving from raw increment to both delta biomass AND RGR 
RGR_plot <- all_trees_weeks %>%
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
  )) %>% 
  #group_by(species) %>% 
  #fill(a) %>% fill(a, .direction = "up") %>% 
  #fill(b) %>% fill(b, .direction = "up") %>% ungroup() %>%
  filter(week == 1 | notes == "new band" | notes == "new band; 182.2=after adjustment") %>%
  mutate(biomass_a= a*DBH_cm^b) %>% #initial biomass
  mutate(bandread_i= bands_cm) %>% #november 2018 band read 
  right_join(all_trees_weeks, by = c("subplot", "species", "tag", "DBH_cm", "week", 
                                     "bands_in", "DOY", "notes", "fate", "bands_cm", "uniqueID"),
             suffix = c("", ".y")) %>% #joining filtered table w/new biomass_a col
  group_by(tag) %>%
  fill(bandread_i) %>% fill(a) %>% fill(b) %>% #fill in missing values
  mutate(inc_cm = case_when(
    species == "ACRU" ~ (bands_cm - lag(bands_cm, default = first(bands_cm)))*0.7647059,
    species == "ACSA" ~ (bands_cm - lag(bands_cm, default = first(bands_cm)))*0.5857143,
    TRUE              ~ bands_cm - lag(bands_cm, default = first(bands_cm))
  )) %>% 
  mutate(DBH_cm_new = case_when(
    !is.na(biomass_a) ~ DBH_cm)) 

#mutate(band_diff= bands_cm-bandread_i) # comment this out for swelling adjutments

#create function for band_diff (how much it has grown since initial bandread)
# banddiff <- function(bands_cm, bandread_i, species, fate){
#   if (species == "ACRU" & fate == "kill"){
#     band_diff <- (bands_cm - bandread_i) * 0.7647059
#   } else if (species == "ACSA" & fate == "kill"){
#     band_diff <- (bands_cm - bandread_i) * 0.5857143
#   } else {
#     band_diff <- bands_cm - bandread_i
#   }
#   return(band_diff)
# }

##########################This is all an attempt to calc band diff is a more effcient way 
# create an inc_cm column with 0's for all week 1's because there is not growth
# anywhere there is a biomass_a means it was the first measurement; therefore initial
# inc is 0 
# initial_inc <- function(biomass_a){
#   if (!is.na(biomass_a)){
#     inc_cm <- 0
#   } else {
#     inc_cm <- NA
#   }
#   return(inc_cm)
# }
# 
# # loop initial inc through RGR_plot to create the inc_cm column 
# for (i in 1:nrow(RGR_plot)){
#   RGR_plot$inc_cm[i] <- initial_inc(RGR_plot$biomass_a[i])
# }
# 
# split RGR_plot into tags to apply the increment growth functions and loops
RGR_plot <- data.frame(RGR_plot)
RGR_plot.list <- split(RGR_plot, RGR_plot$tag)
# 
increment <- function(x) {
  for (i in 1:nrow(x)) {
    for(j in 2:nrow(x)){
      x$DBH_cm_new[j] = x$DBH_cm_new[j - 1] + x$inc_cm[j]
    }
    return(x)
  }
}
# 
RGR_plot.list_grow <- lapply(RGR_plot.list, increment)
# 
RGR_plot <- plyr::ldply(RGR_plot.list_grow, data.frame)

# lag(biomass_b_kg, default = first(biomass_b_kg)))

# loop the function through RGR_plot to create new band_diff column adjusted for
# swelling in the "kill" ACRU
# for (i in 1:nrow(RGR_plot)){
#   RGR_plot$band_diff[i] <- banddiff(RGR_plot$bands_cm[i], RGR_plot$bandsread_i[i],
#                                     RGR_plot$species[i], RGR_plot$fate[i])
# }

# continue pipe to get RGR 
RGR_plot <- RGR_plot %>% 
  mutate(biomass_b_kg= a*(DBH_cm_new)^b) %>%
  mutate(increment_kg= biomass_b_kg-lag(biomass_b_kg, default = first(biomass_b_kg))) %>%
  mutate(RGR= increment_kg/biomass_b_kg) %>% 
  mutate(inc_days= abs(DOY-lag(DOY, default = first(DOY)))) 

###commented out previous pipeline, not accounting for swell ######
# RGR_plot <- all_trees_weeks %>%
#   group_by(species) %>% 
#   fill(a) %>% fill(a, .direction = "up") %>% 
#   fill(b) %>% fill(b, .direction = "up") %>% ungroup() %>%
#   filter(week == 1 | notes == "new band" | notes == "new band; 182.2=after adjustment") %>%
#   mutate(biomass_a= a*DBH_cm^b) %>% #initial biomass
#   mutate(bandread_i= bands_cm) %>% #november band read 
#   right_join(all_trees_weeks, by = c("subplot", "species", "tag", "DBH_cm", "week", 
#                                               "bands_in", "DOY", "notes", "bands_cm", "uniqueID"),
#                       suffix = c("", ".y")) %>% #joining filtered table w/new biomass_a col
#   group_by(tag) %>%
#   fill(bandread_i) %>% fill(a) %>% fill(b) %>% #fill in missing values 
#   mutate(band_diff= bands_cm-bandread_i) %>% # maybe add in swell adj here 
#   mutate(biomass_b_kg= a*(DBH_cm+band_diff)^b) %>%
#   mutate(increment_kg= biomass_b_kg-lag(biomass_b_kg, default = first(biomass_b_kg))) %>%
#   mutate(RGR= increment_kg/biomass_b_kg) %>% 
#   mutate(inc_days= abs(DOY-lag(DOY, default = first(DOY))))

# remove negatives before generating plot-week mean RGR for each species
#RGR_plot$RGR[RGR_plot$RGR < 0] <- 0

# generate mean RGR for each subplot, wk, species to use as RGR for unsampled trees
RGR_mean_subplot_wk_sp <- RGR_plot %>%
  filter(!is.na(RGR)) %>% #filter for a measured trees that have a RGR
  group_by(subplot, week, species) %>% 
  #summarise(mean_delta_biomass= mean(increment_kg)) %>% 
  summarise(mean_RGR_sp= mean(RGR))

# generate mean RGR for each subplot, week to use as RGR for unsampled trees for which
# there is no sampled species in the subplot (i.e. i did not sample oaks in A03E, but 
# there ARE oaks in A03E)
RGR_mean_subplot_wk <- RGR_plot %>%
  filter(!is.na(RGR)) %>% #filter for a measured trees that have a RGR
  group_by(subplot, week) %>% 
  summarise(mean_RGR_plot= mean(RGR))

# join above df's with RGR_plot so that each tree has either mean species RGR, mean_plot RGR,
# or measured RGR
RGR_plot_mean_sp <- right_join(RGR_mean_subplot_wk_sp, RGR_plot)
RGR_plot_mean <- right_join(RGR_mean_subplot_wk, RGR_plot_mean_sp)

##################################################################################


#create new RGR_obs_est column for calculated RGR (sampled trees), avg_RGR (unsampled trees), 
#adj_RGR (for unique id's where DBH had significant effect on RGR)
RGR_obs_est_df <- RGR_plot_mean %>% 
  mutate(RGR_obs_est = case_when(
    !is.na(RGR) ~ RGR,
    is.na(RGR) & !is.na(mean_RGR_sp) ~ mean_RGR_sp,
    is.na(RGR) & is.na(mean_RGR_sp) ~ mean_RGR_plot))

RGR_obs_est_df <- RGR_obs_est_df %>% # biomass_b_kg should be used for sampled trees
  mutate(biomass_b_kg_new = biomass_a + biomass_a*RGR_obs_est) %>%
  select(uniqueID, week, subplot, tag, species, DBH_cm, DBH_cm_new, biomass_a,
         biomass_b_kg_new, fate, RGR_obs_est, DOY, inc_days) %>%
  arrange(tag)

df <- data.frame(RGR_obs_est_df)
# this should set the RGR to 0 if negative
# df$RGR_obs_est[df$RGR_obs_est < 0] <- 0
#split into list of data frames
df.list <- split(df, df$tag)
df.list <- lapply(df.list, function(x){
  for (i in 1:nrow(x)) {
    for(j in 2:nrow(x)){
      x$biomass_a[j] = x$biomass_b_kg_new[j - 1]
      x$biomass_b_kg_new[j] = (x$biomass_a[j] * x$RGR_obs_est[j]) + x$biomass_a[j]
    }
    return(x)
  }
})
df <- plyr::ldply(df.list, data.frame)

# let's get some NPP bebe
inc_days <- RGR_plot %>% 
  group_by(subplot, week) %>% # MUST ADD SPECIES HERE TOO!!! 
  summarise(inc_days= mean(inc_days))

NPP_final <- df %>% 
  filter(!is.na(biomass_a)) %>% #NA's are from species that are not present in a subplot; need to remove them 
  mutate(biomass_diff = biomass_b_kg_new - biomass_a) %>% 
  group_by(subplot, week) %>% # can ADD SPECIES in here if i want to look at effect of species 
  summarise(subplot_biomass_kg = sum(biomass_diff), DOY = mean(DOY)) %>% 
  mutate(biomass_per_ha = subplot_biomass_kg*10) %>% 
  mutate(kgC_per_ha = biomass_per_ha*0.48) %>% 
  right_join(inc_days) %>% 
  mutate(kgC_per_ha_day = kgC_per_ha / inc_days) %>% # this is kgC/ha/day; call NPP_can for the gather function
  mutate(replicate = substr(subplot, 1, 1)) %>% 
  filter(!is.na(kgC_per_ha_day)) %>% 
  mutate(severity = case_when(
    subplot == "A01E" ~ 0.85, subplot == "A01W" ~ 0.85, subplot == "A02E" ~ 0.45,
    subplot == "A02W" ~ 0.45, subplot == "A03E" ~ 0.65, subplot == "A03W" ~ 0.65,
    subplot == "A04E" ~ 0.00, subplot == "A04W" ~ 0.00, subplot == "B01E" ~ 0.00,
    subplot == "B01W" ~ 0.00, subplot == "B02E" ~ 0.45, subplot == "B02W" ~ 0.45,
    subplot == "B03E" ~ 0.85, subplot == "B03W" ~ 0.85, subplot == "B04E" ~ 0.65,
    subplot == "B04W" ~ 0.65, subplot == "C01E" ~ 0.00, subplot == "C01W" ~ 0.00,
    subplot == "C02E" ~ 0.65, subplot == "C02W" ~ 0.65, subplot == "C03E" ~ 0.85,
    subplot == "C03W" ~ 0.85, subplot == "C04E" ~ 0.45, subplot == "C04W" ~ 0.45, 
    subplot == "D01E" ~ 0.00, subplot == "D01W" ~ 0.00, subplot == "D02E" ~ 0.85,
    subplot == "D02W" ~ 0.85, subplot == "D03E" ~ 0.45, subplot == "D03W" ~ 0.45,
    subplot == "D04E" ~ 0.65, subplot == "D04W" ~ 0.65
  )) %>% 
  mutate(treatment = case_when(
    subplot == "A01E" ~ "bottom", subplot == "A01W" ~ "top", subplot == "A02E" ~ "top",
    subplot == "A02W" ~ "bottom", subplot == "A03E" ~ "bottom", subplot == "A03W" ~ "top",
    subplot == "A04E" ~ "bottom", subplot == "A04W" ~ "top", subplot == "B01E" ~ "bottom",
    subplot == "B01W" ~ "top", subplot == "B02E" ~ "top", subplot == "B02W" ~ "bottom",
    subplot == "B03E" ~ "bottom", subplot == "B03W" ~ "top", subplot == "B04E" ~ "top",
    subplot == "B04W" ~ "bottom", subplot == "C01E" ~ "top", subplot == "C01W" ~ "bottom",
    subplot == "C02E" ~ "bottom", subplot == "C02W" ~ "top", subplot == "C03E" ~ "bottom",
    subplot == "C03W" ~ "top", subplot == "C04E" ~ "top", subplot == "C04W" ~ "bottom", 
    subplot == "D01E" ~ "bottom", subplot == "D01W" ~ "top", subplot == "D02E" ~ "top",
    subplot == "D02W" ~ "bottom", subplot == "D03E" ~ "bottom", subplot == "D03W" ~ "top",
    subplot == "D04E" ~ "bottom", subplot == "D04W" ~ "top"
  )) 

# getting DOY into dates and correcting the november 2018 date 
dates <- data.frame(as.factor(as.Date(NPP_final$DOY, "2019-01-01")))
names(dates) <- c("date")
dates <- data.frame(recode_factor(dates$date, '2019-11-03' = "2018-11-15"))
names(dates) <- c("date")
dates$date <- as.Date(dates$date,"%Y-%m-%d")

NPP_final <- NPP_final %>% 
  bind_cols(dates) 

#create week and month columns to look at different temporal scales 
NPP_final$month <- as.Date(cut(NPP_final$date, breaks = "month"))
NPP_final$week <- as.Date(cut(NPP_final$date, breaks = "week", start.on.monday = FALSE))

#rename kgC_per_ha_day to NPP_canopy
# NPP_final <- rename(NPP_final, NPP_canopy = kgC_per_ha_day)

################3
#scale whole season increment (nov. 2018 to Nov 2019) for cumulative NPP
# NPP_cumulative <- df %>% 
#   filter(!is.na(biomass_a)) %>% #NA's are from species that are not present in a subplot; need to remove them 
#   mutate(biomass_diff = biomass_b_kg_new - biomass_a) %>% 
#   group_by(subplot) %>% 
#   summarise(subplot_biomass_kg = sum(biomass_diff))%>% 
#   mutate(biomass_per_ha = subplot_biomass_kg*10) %>% 
#   mutate(annualNPP = biomass_per_ha*0.48) %>% 
#   mutate(replicate = substr(subplot, 1, 1)) 

# calculate cumulative NPP starting with the yearly (nov 2018-nov 2019) increment
annual_inc <- df %>% 
  filter(week == 1) %>% 
  mutate(biomass_i = biomass_a) %>% #biomass_i is nov 2018 biomass
  right_join(df) %>% 
  group_by(tag) %>% 
  fill(biomass_i) %>% 
  select(subplot, week, tag, species, DOY, biomass_i) %>% #, biomass_b_kg_new
  filter(week == 17) %>% 
  mutate(annual_biomass = biomass_b_kg_new - biomass_i) %>% ungroup() %>% 
  group_by(subplot) %>% 
  summarise(annual_subplot = sum(annual_biomass)) %>% 
  mutate(biomass_per_ha = annual_subplot*10) %>% 
  mutate(kgC_ha_yr = biomass_per_ha* 0.48) %>% 
  mutate(severity = case_when(
    subplot == "A01E" ~ 0.85, subplot == "A01W" ~ 0.85, subplot == "A02E" ~ 0.45,
    subplot == "A02W" ~ 0.45, subplot == "A03E" ~ 0.65, subplot == "A03W" ~ 0.65,
    subplot == "A04E" ~ 0.00, subplot == "A04W" ~ 0.00, subplot == "B01E" ~ 0.00,
    subplot == "B01W" ~ 0.00, subplot == "B02E" ~ 0.45, subplot == "B02W" ~ 0.45,
    subplot == "B03E" ~ 0.85, subplot == "B03W" ~ 0.85, subplot == "B04E" ~ 0.65,
    subplot == "B04W" ~ 0.65, subplot == "C01E" ~ 0.00, subplot == "C01W" ~ 0.00,
    subplot == "C02E" ~ 0.65, subplot == "C02W" ~ 0.65, subplot == "C03E" ~ 0.85,
    subplot == "C03W" ~ 0.85, subplot == "C04E" ~ 0.45, subplot == "C04W" ~ 0.45, 
    subplot == "D01E" ~ 0.00, subplot == "D01W" ~ 0.00, subplot == "D02E" ~ 0.85,
    subplot == "D02W" ~ 0.85, subplot == "D03E" ~ 0.45, subplot == "D03W" ~ 0.45,
    subplot == "D04E" ~ 0.65, subplot == "D04W" ~ 0.65
  )) %>% 
  mutate(treatment = case_when(
    subplot == "A01E" ~ "bottom", subplot == "A01W" ~ "top", subplot == "A02E" ~ "top",
    subplot == "A02W" ~ "bottom", subplot == "A03E" ~ "bottom", subplot == "A03W" ~ "top",
    subplot == "A04E" ~ "bottom", subplot == "A04W" ~ "top", subplot == "B01E" ~ "bottom",
    subplot == "B01W" ~ "top", subplot == "B02E" ~ "top", subplot == "B02W" ~ "bottom",
    subplot == "B03E" ~ "bottom", subplot == "B03W" ~ "top", subplot == "B04E" ~ "top",
    subplot == "B04W" ~ "bottom", subplot == "C01E" ~ "top", subplot == "C01W" ~ "bottom",
    subplot == "C02E" ~ "bottom", subplot == "C02W" ~ "top", subplot == "C03E" ~ "bottom",
    subplot == "C03W" ~ "top", subplot == "C04E" ~ "top", subplot == "C04W" ~ "bottom", 
    subplot == "D01E" ~ "bottom", subplot == "D01W" ~ "top", subplot == "D02E" ~ "top",
    subplot == "D02W" ~ "bottom", subplot == "D03E" ~ "bottom", subplot == "D03W" ~ "top",
    subplot == "D04E" ~ "bottom", subplot == "D04W" ~ "top"
  ))%>% 
  group_by(severity) %>% 
  summarise(annual_NPP = mean(kgC_ha_yr))

# make severity %'s for x axis
annual_inc$severity <- recode(annual_inc$severity, "0.00" = "0%",
                              "0.45" = "45%", "0.65" = "65%",
                              "0.85" = "85%")
