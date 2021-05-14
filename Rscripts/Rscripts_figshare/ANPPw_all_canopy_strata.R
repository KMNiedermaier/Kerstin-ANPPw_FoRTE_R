#set library to dyplr
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotrix)

# use source to run subcanopy and seedling_saplings scripts 
# edit path (...) to access source scripts and project data 
source(".../subcanopy.R")
source(".../seedling_saplings.R")

#importing csv
dendro_data <- read.csv(".../canopy_dendrobands.csv")

all_reps_all_trees <- read.csv(".../FoRTE_all_trees.csv")

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

# make a df with same columns as all_trees but with only sampled trees to compare to
# all_trees and filter for only unsampled trees below 
#sample_unsampled <- select(all_reps, subplot, species, tag, fate,)

#creating 16 of each tag to build a week column and prep for join; then join and fill in DOY 
all_trees_weeks <- all_trees %>% 
  #dplyr::setdiff(all_reps)
  slice(rep(1:n(), each = 17)) %>% 
  group_by(tag) %>% 
  mutate(week = row_number()) %>% 
  left_join(all_reps) %>% 
  group_by(subplot, week) %>% 
  fill(DOY) %>%
  fill(DOY, .direction = c("up")) %>% 
  fill(notes) %>% 
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
  filter(week == 1 | notes == "new band" | notes == "new band; 182.2=after adjustment") %>%
  mutate(biomass_a= a*DBH_cm^b) %>% #initial biomass
  mutate(bandread_i= bands_cm) %>% #november 2018 band read 
  right_join(all_trees_weeks, by = c("subplot", "species", "tag", "DBH_cm", "week", 
                                     "bands_in", "DOY", "notes", "fate", "bands_cm", "uniqueID"),
             suffix = c("", ".y")) %>% #joining filtered table w/new biomass_a col
  group_by(tag) %>%
  fill(bandread_i) %>% fill(a) %>% fill(b) %>% #fill in missing values
  mutate(inc_cm = case_when(
    species == "ACRU" ~ ((bands_cm - lag(bands_cm, default = first(bands_cm)))*0.7647059)/pi,
    species == "ACSA" ~ ((bands_cm - lag(bands_cm, default = first(bands_cm)))*0.5857143)/pi,
    species != "ACSA" & species != "ACRU" ~ 
      (bands_cm - lag(bands_cm, default = first(bands_cm)))/pi)
    ) %>% 
  mutate(DBH_cm_new = case_when(
    !is.na(biomass_a) ~ DBH_cm)) 

RGR_plot$inc_cm[RGR_plot$inc_cm < 0] <- 0

# split RGR_plot into tags to apply the increment growth functions and loops
RGR_plot <- data.frame(RGR_plot)
RGR_plot.list <- split(RGR_plot, RGR_plot$tag)
# 
increment <- function(x) {
  for (i in 1:nrow(x)) {
    for(j in 2:nrow(x)){
      x$DBH_cm_new[j] = x$DBH_cm_new[j - 1] + (x$inc_cm[j])
    }
    return(x)
  }
}
# 
RGR_plot.list_grow <- lapply(RGR_plot.list, increment)
# 
RGR_plot <- plyr::ldply(RGR_plot.list_grow, data.frame)

RGR_plot$DOY <- as.numeric(RGR_plot$DOY)

# continue pipe to get RGR 
RGR_plot <- RGR_plot %>% 
  group_by(tag) %>% 
  mutate(RGR = inc_cm / lag(DBH_cm_new, default = first(DBH_cm_new))) %>% 
  mutate(inc_days= case_when(
    DOY == 98 ~ 157,
    DOY == 99 ~ 158,
    DOY == 100 ~ 159,
    DOY != 98 |DOY != 99 |DOY!= 100 ~ abs(DOY-lag(DOY, default = first(DOY)))
  )) 

# generate mean RGR for each subplot, wk, species to use as RGR for unsampled trees
RGR_mean_subplot_wk_sp <- RGR_plot %>%
  filter(!is.na(RGR)) %>% #filter for a measured trees that have a RGR
  group_by(subplot, week, species) %>% 
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
RGR_plot_mean <- RGR_plot_mean %>% 
  group_by(tag) %>% fill(notes)

##################################################################################
##### sampled trees to generate regression equations to estimate growth 
##### of unsampled population. regressions will also be used to identify and adjust 
##### growth for the size class effect on growth (i.e. bigger trees grow more wood)

#selected only for columns that I will use 
all_reps <- select(dendro_data, severity, fate, treatment, week, subplot, tag, replicate, 
                   species, DBH_cm, bands_in, DOY, notes) %>%
  mutate(bands_cm= bands_in*2.54) %>% 
  filter(!is.na(bands_cm))

#starting biomass and band read; 
all_reps_filtered <- all_reps %>%
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
  group_by(tag) %>%
  filter(week == 1 | notes == "new band" | notes == "new band; 182.2=after adjustment") %>%
  mutate(biomass_a= a*DBH_cm^b) %>% #initial biomass
  mutate(bandread_i= bands_cm) %>% #november band read 
  right_join(all_reps) %>% #joining filtered table w/new biomass col
  group_by(tag) %>%
  fill(bandread_i) %>% fill(a) %>% fill(b) %>%
  mutate(inc_cm = case_when(
    species == "ACRU" ~ ((bands_cm - lag(bands_cm, default = first(bands_cm)))*0.7647059)/pi,
    species == "ACSA" ~ ((bands_cm - lag(bands_cm, default = first(bands_cm)))*0.5857143)/pi,
    species != "ACSA" & species != "ACRU" ~ 
      (bands_cm - lag(bands_cm, default = first(bands_cm)))/pi)) %>% 
  mutate(DBH_cm_new = case_when(
    !is.na(biomass_a) ~ DBH_cm)) 

all_reps_filtered$inc_cm[all_reps_filtered$inc_cm < 0] <- 0
  
# split RGR_plot into tags to apply the increment growth functions and loops
all_reps_filtered <- data.frame(all_reps_filtered)
all_reps_filtered.list <- split(all_reps_filtered, all_reps_filtered$tag)
# 
increment <- function(x) {
  for (i in 1:nrow(x)) {
    for(j in 2:nrow(x)){
      x$DBH_cm_new[j] = x$DBH_cm_new[j - 1] + (x$inc_cm[j])
    }
    return(x)
  }
}
# 
all_reps_filtered.list_grow <- lapply(all_reps_filtered.list, increment)
# 
all_reps_filtered <- plyr::ldply(all_reps_filtered.list_grow, data.frame)

all_reps_filtered1 <- all_reps_filtered %>% #calc inc percents 
  group_by(tag) %>%
  mutate(RGR = inc_cm / lag(DBH_cm_new, default = first(DBH_cm_new))) 

#select for RGR, tag, date, subplot, species, DBH
RGR_plot_regression <- all_reps_filtered1 %>%
  select(subplot, week, tag, species, 
         DBH_cm_new, DOY, biomass_a, RGR) %>% #biomass_b_kg, , increment_kg
  group_by(tag) 


RGR_plot$uniqueID <- paste(RGR_plot$subplot, RGR_plot$week,
                           RGR_plot$species, sep = "_")

#scale to NPP by subplot, week, species, AND DBH (where necessary)
#runs regression for each subplot, week, and species AND returns object with 
#coefficients
regressions <- RGR_plot_regression %>% 
  group_by(subplot, week, species) %>% 
  do(mod = summary(lm(RGR ~ DBH_cm_new, data = .))) #%>% #-1 removes the intercept

regressions$uniqueID <- paste(regressions$subplot, regressions$week, 
                              regressions$species, sep = "_")

final_output <- matrix(nrow= 2466, ncol=4, dimnames = list(c(), 
                                                           c("uniqueID", "pvalue", 
                                                             "slope", "intercept")))
regressions$uniqueID <- paste(regressions$subplot, regressions$week, 
                              regressions$species, sep = "_")

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

#create new RGR_obs_est column for calculated RGR (sampled trees), avg_RGR (unsampled trees), 
#adj_RGR (for unique id's where DBH had significant effect on RGR)
RGR_obs_est_df <- RGR_plot_mean %>% 
  left_join(sig_DBH_effect) %>% 
  mutate(RGR_obs_est = case_when(
    !is.na(RGR) ~ RGR,
    notes == "felled" ~ 0,
    notes == "dead" ~ 0,
    is.na(RGR) & !is.na(mean_RGR_sp) & is.na(slope) ~ mean_RGR_sp,
    is.na(RGR) & is.na(mean_RGR_sp) & is.na(slope) ~ mean_RGR_plot,
    is.na(RGR) & !is.na(slope) ~ DBH_cm*slope+intercept))

df <- data.frame(RGR_obs_est_df)
# this should set the RGR to 0 if negative
df$RGR_obs_est[df$RGR_obs_est < 0] <- 0
#split into list of data frames
df.list <- split(df, df$tag)
df.list <- lapply(df.list, function(x){
  for (i in 1:nrow(x)) {
    for(j in 2:nrow(x)){
      x$DBH_cm_new[j] = (x$DBH_cm_new[j - 1] * x$RGR_obs_est[j]) + x$DBH_cm_new[j-1]
      #x$biomass_a[j] = x$biomass_b_kg_new[j - 1]
      #x$biomass_b_kg_new[j] = (x$biomass_a[j] * x$RGR_obs_est[j]) + x$biomass_a[j]
    }
    return(x)
  }
})
df <- plyr::ldply(df.list, data.frame)

# let's get some NPP bebe
inc_days <- RGR_plot %>% 
  group_by(subplot, week) %>% # MUST ADD SPECIES HERE TOO!!! 
  summarise(inc_days= mean(inc_days)) 

# weekly NPP per subplot 
NPP_final <- df %>% 
  group_by(tag) %>% 
  mutate(biomass_new = case_when(
    species == "AMEL" ~ (a*(DBH_cm_new*10)^b)/1000,
    species != "AMEL" ~ a*DBH_cm_new^b
  )) %>%  
  group_by(subplot, week) %>% # can ADD SPECIES in here if you want to look at effect of species 
  summarise(subplot_biomass_kg = sum(biomass_new), DOY = mean(DOY)) %>% 
  mutate(biomass_diff = subplot_biomass_kg - lag(subplot_biomass_kg, 
                                                 default = first(subplot_biomass_kg))) %>% 
  mutate(biomass_per_ha = biomass_diff*10) %>% 
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
    subplot == "D01E" ~ "bottom", subplot == "D01W" ~ "top", subplot == "D02E" ~ "bottom",
    subplot == "D02W" ~ "top", subplot == "D03E" ~ "bottom", subplot == "D03W" ~ "top",
    subplot == "D04E" ~ "top", subplot == "D04W" ~ "bottom"
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

# severity time series df; used for generating figures 
severity_time <- NPP_final %>% 
  group_by(week, severity) %>% 
  summarize(NPP = mean(kgC_per_ha_day), SE = std.error(kgC_per_ha_day)) %>% 
  mutate(DOY = yday(week))

# annual (cumulative) NPP; essentially the sum of all the weeks in each subplot 
annual_inc <- df %>% 
  group_by(tag) %>% # NEW WIH inc_cm
  mutate(biomass_new = case_when(
    species == "AMEL" ~ (a*(DBH_cm_new*10)^b)/1000,
    species != "AMEL" ~ a*DBH_cm_new^b
  )) %>% 
  mutate(biomass_diff = biomass_new - lag(biomass_new, default = first(biomass_new))) %>% 
  group_by(subplot) %>% 
  summarise(annual_subplot = sum(biomass_diff)) %>%
  mutate(biomass_per_ha = annual_subplot*10) %>% 
  mutate(kgC_ha_yr = biomass_per_ha* 0.48) %>% 
  mutate(replicate = substr(subplot, 1, 1)) %>% 
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
    subplot == "D01E" ~ "bottom", subplot == "D01W" ~ "top", subplot == "D02E" ~ "bottom",
    subplot == "D02W" ~ "top", subplot == "D03E" ~ "bottom", subplot == "D03W" ~ "top",
    subplot == "D04E" ~ "top", subplot == "D04W" ~ "bottom"
  )) %>% 
  group_by(subplot, severity, treatment) %>% 
  summarise(NPP_canopy = mean(kgC_ha_yr))

# make severity %'s for x axis
annual_inc$severity <- recode(annual_inc$severity, "0.00" = "0",
                              "0.45" = "45", "0.65" = "65",
                              "0.85" = "85")

###################################################################################
########### This is joining subcanopy and canopy for the cumulative figure
# join NPP subcanopy and seedlings data, create strata column and ID NPP by canopy,
# subcanopy, and seedling
annual_inc1 <- annual_inc %>%
  left_join(NPP_subcanopy) %>% 
  left_join(seedlings_total) %>% 
  mutate(annual_NPP = NPP_canopy + NPP_subcan + NPP_seedlings) %>% 
  ungroup()

annual_inc1$treatment <- recode(annual_inc1$treatment, bottom = "zbottom")

#annual NPP by treatments (for figures)
annual_NPP_tx <- annual_inc1 %>% 
  group_by(severity, treatment) %>% 
  summarise(NPP = mean(annual_NPP), SE = std.error(annual_NPP)) %>% 
  group_by(severity) %>% 
  mutate(errorbars = case_when(
    treatment == "zbottom" ~ NPP,
    treatment == "top" ~ sum(NPP)
  ))
  
# Create df for annual NPP figure 
annual_NPP_fig <- annual_inc1 %>% 
  group_by(severity) %>% 
  summarise(NPP_canopys = mean(NPP_canopy), SE1 = std.error(NPP_canopy),
            NPP_subcans = mean(NPP_subcan), SE2 = std.error(NPP_subcan),
            NPP_zseedlings = mean(NPP_seedlings), SE3 = std.error(NPP_seedlings)) %>% 
  gather(canopy_strata, annual_NPP, NPP_canopys, NPP_subcans, NPP_zseedlings) %>% 
  mutate(SE = case_when(
    canopy_strata == "NPP_zseedlings" ~ SE3,
    canopy_strata == "NPP_subcans" ~ SE2,
    canopy_strata == "NPP_canopys" ~ SE1
  )) %>% 
  group_by(severity) %>% 
  arrange(severity) %>% 
  mutate(errorbars = case_when(
    canopy_strata == "NPP_zseedlings" ~ annual_NPP,
    canopy_strata == "NPP_subcans" ~ sum(annual_NPP) - lag(annual_NPP, default = first(annual_NPP)),
    canopy_strata == "NPP_canopys" ~ sum(annual_NPP)
  ))
