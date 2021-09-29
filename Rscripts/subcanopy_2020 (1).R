# This script brings in subcanopy stem counts (to the species level) that were 
# perfromed in one quadrant (0.025ha) of each subplot, and the subcanopy repeated
# diameter measures from the 2019 growing season. It calculates mean annual
# subplot subcanopy ANPPw from the growth increment between diameter measurements 
# in August 2019 to Nov. 2020

#load necessary packages 
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(plotrix)

# run the subcanopy 219 script (1-8cm_subcanopy.R)
source("Rscripts/1-8cm_subcanopy.R")

#remove clutter from global environmnet from the subcanopy and seedling scripts
remove(subcan_comp, 
       subcan_stem_density, subcanopy_data, subcanopy_select, tree_counts, 
       sub_annual_inc, subcan_comp1, may_2019, 
       aug_2019)

# importing tree count data 
tree_counts <- read.csv("data/subcanopy_stemcounts.csv")

#importing csv of subcanopy diameter measurements (sc = subcanopy)
sc_2019 <- read.csv("data/subcanopy_D.csv")

# filter for august 2019 (latest season measurement for 2019)
sc_aug_2019 <- filter(sc_2019, date == "2019-08-03" | date == "2019-08-04") %>% 
  select(uniqueID, quadrat, species, tag, DBH_mm, date)

# import 2020 subcanopy data 
sc_2020 <- read.csv("data/subcanopy_D_2020.csv")

# filter 2020 data for the Nov. data 
sc_nov_2020 <- filter(sc_2020, date == "2020-11-18" | date == "2020-11-19") %>% 
  select(-X)

# bind aug 2019 and nov 2020
sc_19_20 <- bind_rows(sc_aug_2019, sc_nov_2020)

##################################################################################
# Scale diameter increment to subcanopy ANPPw

# create function for biomass based on allometric parameters (Cooper, 1981) 
# for each species 

biomass_a <- function(species, DBH){
  if (species == "ACRU"){
    biomass <- 0.03117 * (DBH) ^ 2.7780
  } else if (species == "ACPE"){
    biomass <-  0.2040 * DBH ^ 2.2524
  } else if (species == "ACSA"){
    biomass <- 0.1693 * DBH ^ 2.3436
  } else if (species == "AMEL"){
    biomass <- (0.1630 * (DBH * 10) ^ 2.4940)/1000
  } else if (species == "FAGR"){
    biomass <- 0.1892 * DBH ^ 2.3097
  } else if (species == "PIRE"){
    biomass <- 0.0526 * DBH ^ 2.5258
  } else if (species == "PIST"){
    biomass <- 0.0408 * DBH ^ 2.5735
  } else if (species == "POGR"){
    biomass <- 0.1387 * DBH ^ 2.3498
  } else if (species == "QURU"){
    biomass <- 0.0398 * DBH ^ 2.7734
  }
  return(biomass)
}

# vectorize my function becuase it will be used on a vector 
biomass_a <- Vectorize(biomass_a, vectorize.args = c("species", "DBH"))

# remove records where DBH_mm
# was NA (some stems are missing week 1 and 2 of data collection). Next I convert
# DBH in mm to cm, and use my biomass_a function to estimate biomass(kg)
# (DBH allometry) 
sc_19_20_biomass <- sc_19_20 %>% 
  filter(!is.na(DBH_mm)) %>% 
  mutate(DBH_cm = DBH_mm/10) %>%
  mutate(biomass_kg = biomass_a(species, DBH_cm)) 

# create a weeks column based on the date 
sc_19_20_biomass$date <- as.Date(sc_19_20_biomass$date,"%Y-%m-%d")
sc_19_20_biomass$week <- as.Date(cut(sc_19_20_biomass$date, breaks = "week", 
                                     start.on.monday = TRUE))

# now I arrange and group for organization. Next, I do
# some vector math to calculat the biomass increment (biomass_inc). Then, I filter 
# for the end of season measuremnt (no increment for first measurement), and again, 
# select and arrange for organization. Laslty, create a subplot column from uniqueID
sc_20_inc <- sc_19_20_biomass %>% 
  mutate(subplot = substr(uniqueID, 1, 4)) %>% 
  arrange(tag) %>% group_by(tag) %>% 
  mutate(biomass_inc = biomass_kg - lag(biomass_kg, default = first(biomass_kg))) %>% 
  filter(row_number() == n()) %>% 
  select(uniqueID, species, tag, biomass_kg, biomass_inc) %>% 
  arrange(uniqueID) %>% 
  mutate(subplot = substr(uniqueID, 1, 4))

# get rid of the negatives in biomass increments. 
# Negatives likely result from inherent error in measuring diameter with 
# digital calipers 
sc_20_inc$biomass_inc[sc_20_inc$biomass_inc < 0] <- 0

# Now I will move from biomass increment to subcanopy NPP. To scale from my sampled 
# population to the entire population, I first calculate and create a new vector 
# with mean subplot biomass increment. This is then used for subcanopy species that 
# are within the subplot, but that I did not capture in my sample. I then join the 
# tree counts df which creates a new record within subplots for species that were 
# not captures in my sample population. Next, I fill the mean subplot increments, and
# create a new vector that has either the measured increment (sampled species) or the 
# mean increment (unsampled species). Then, I summarize for mean_biomass for each
# subplot, and species, rejoin tree_counts (lost in the summarize), and scale up to 
# the hectare. Lastly, I sum the biomass production of all species within a subplot
# and multiply by a site specific carbon fraction of 0.48 for kgC_ha_yr
NPP_sc_2020 <- sc_20_inc %>% 
  group_by(subplot) %>% 
  summarize(mean_subplot = mean(biomass_inc)) %>% 
  left_join(sc_20_inc) %>% 
  right_join(tree_counts) %>% arrange(subplot, species) %>% 
  group_by(subplot) %>% fill(mean_subplot, .direction = "updown") %>% 
  group_by(subplot, species) %>% 
  mutate(biomass_inc_obs_est = case_when(
    !is.na(biomass_inc) ~ biomass_inc,
    is.na(biomass_inc) ~ mean_subplot
  )) %>% 
  summarize(mean_biomass = mean(biomass_inc_obs_est)) %>% 
  right_join(tree_counts) %>% 
  mutate(sp_biomass_ha = mean_biomass * count * 40) %>% 
  # * species count * 40 b/c counts were in 0.025ha
  group_by(subplot) %>% 
  summarize(kgC_ha_yr = sum(sp_biomass_ha) * 0.48) %>% 
  rename(NPP_2020 = kgC_ha_yr) %>% 
  left_join(NPP_sc_2019) %>% 
  rename(NPP_2019 = NPP_subcan) %>% 
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
  ungroup()

# now that NPP_sc_2020 was joing with 2019, reorganize to form a year column
NPP_sc_alltime <- NPP_sc_2020 %>% 
  gather(year, NPP1, NPP_2019, NPP_2020) 

# recode th eyear column to contain only numbers
NPP_sc_alltime$year <- recode(NPP_sc_alltime$year, NPP_2019 = "2019", NPP_2020 = "2020")

# make severity %'s for x axis
NPP_sc_alltime$severity <- recode(NPP_sc_alltime$severity, "0.00" = "0",
                               "0.45" = "45", "0.65" = "65",
                               "0.85" = "85")

# plot severity by year with NPP as response variable. NPP is KgC/ha/year
ggplot(NPP_sc_alltime, aes(x = factor(severity), y = NPP1, fill = factor(year))) +
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x = element_text(size = 17), axis.text.y = 
          element_text(size = 17), axis.title.x = element_text(size = 17),
          axis.title.y = element_text(size = 17), legend.text = element_text(size = 15),
          legend.title = element_blank()) +
  labs(x = "Disturbance Severity (%)",
    y = expression(paste("ANPP" [w], " ( ",kgC," ",ha^-1," ",day^-1,")"))) 

#ggsave("figures/subcanopyNPP.jpeg", height = 6, width = 8, units = "in", last_plot())
  
# same plot as above but for treatment/disturbance type 
ggplot(NPP_sc_alltime, aes(x = factor(treatment), y = NPP1, fill = factor(year))) +
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x = element_text(size = 17), axis.text.y = 
          element_text(size = 17), axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17), legend.text = element_text(size = 15),
        legend.title = element_blank()) +
  labs(x = "Disturbance Type",
       y = expression(paste("ANPP" [w], " ( ",kgC," ",ha^-1," ",day^-1,")")))

#####################################################################################
# close and clean up
#remove clutter from global environmnet from the subcanopy and seedling scripts
remove(sc_19_20, sc_19_20_biomass, sc_20_inc, sc_2019, sc_2020, sc_aug_2019, 
       sc_nov_2020, tree_counts)


