# Taken from Max's 1-8cm_subcanopy script #blessings`

#This script brings in subcanopy stem counts (to the species level) that were 
# performed in one quadrant (0.025ha) of each subplot, and the subcanopy repeated
# diameter measures from the 2019 growing season. It calculates mean annual
# subplot subcanopy ANPPw from the growth increment between diameter measurements 
# in May 2019 and August 2019
#I (Kerstin) will edit the script to show changes from 2020 to 2021 

#load necessary packages 
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(plotrix)

# importing tree count data 
tree_counts <- read.csv("data_don't_use/subcanopy_stemcounts.csv")
# select columns of interest
tree_counts <-  select(tree_counts, subplot, species, count)

#importing csv of subcanopy diameter measurements 
subcanopy_data_2021 <- read.csv("2021_data/subcanopy_diameter_2021.csv")
subcanopy_data_2020 <- read.csv("data_don't_use/fd_subcanopy_diameter_2020.csv")


# this chunk of code stands alone (not used in ANPPw calcs, but interesting data!)
# It calculates subcanopy stem density per subplotand scaled to the ha. 
# Then summarize stems/ha in each replicate with SE among subplots within a replicate 
subcan_stem_density <- tree_counts %>% 
  mutate(replicate = substr(subplot, 1, 1)) %>%
  group_by(subplot, replicate) %>% 
  summarize(quart_subplot = sum(count)) %>% 
  mutate(stems_per_ha = quart_subplot*40) %>% #counts were done in 0.025 ha
  group_by(replicate) %>% 
  summarize(rep_stem_density = mean(stems_per_ha), SE = std.error(stems_per_ha))

# This chunk also stands alone. It calculate the subcanopy species composition 
# within each replicate.
subcan_comp1 <- tree_counts %>% 
  mutate(replicate = substr(subplot, 1, 1)) %>% 
  select(subplot, replicate, species, count) %>% 
  group_by(replicate, species) %>% 
  summarise(total_stems = sum(count)) #%>% # this is in 0.025 ha plots (1/4 subplots)

subcan_comp <- subcan_comp1 %>%   
  summarise(total_stems_rep = sum(total_stems)) %>% 
  right_join(subcan_comp1) %>% 
  mutate(perc_comp = total_stems/total_stems_rep) 

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

# vectorize my function becuase it will be used on a vetor 
biomass_a <- Vectorize(biomass_a, vectorize.args = c("species", "DBH"))


# select important columns from subcanopy data. Then remove records where DBH_mm
# was NA (some stems are missing week 1 and 2 of data collection). Next I convert
# DBH in mm to cm, and use my biomass_a function to estimate biomass(kg)
# (DBH allometry) 
subcanopy_select <- subcanopy_data_2021 %>% 
  select(uniqueID, species, tag, DBH_mm, date) %>% 
  filter(!is.na(DBH_mm)) %>% 
  mutate(DBH_cm = DBH_mm/10) %>%
  mutate(biomass_kg = biomass_a(species, DBH_cm)) 

subcanopy_2020_select <- subcanopy_data_2020 %>% 
  select(subplot_id, species, tag, dbh_mm, date) %>% 
  filter(!is.na(dbh_mm)) %>% 
  mutate(dbh_cm = dbh_mm/10) %>% 
  mutate(biomass_kg = biomass_a(species, dbh_cm))

# create a weeks column based on the date 
subcanopy_select$date <- as.Date(subcanopy_select$date,"%Y-%m-%d")
subcanopy_select$week <- as.Date(cut(subcanopy_select$date, breaks = "week", 
                                     start.on.monday = TRUE))
subcanopy_2020_select$date <- as.Date(subcanopy_2020_select$date,"%Y-%m-%d")
subcanopy_2020_select$week <- as.Date(cut(subcanopy_2020_select$date, breaks = "week", 
                                     start.on.monday = TRUE))

# create a DOY column
subcanopy_select$DOY <- yday(subcanopy_select$date)

subcanopy_2020_select$DOY <- yday(subcanopy_2020_select$date)

# here is the first measurement of every sampled tree
nov_2020 <- subcanopy_2020_select %>% 
  group_by(tag) %>% 
  filter(row_number() == 1) %>% 
  mutate(biomass_kg = biomass_a(species, dbh_cm)) 

# here is the final measurement of every sampled tree
jun_2021 <- subcanopy_select %>% 
  group_by(tag) %>% 
  filter(row_number() == n()) %>% 
  mutate(biomass_kg = biomass_a(species, DBH_cm)) 

# now I bind these two dfs together, arrange and group for organization. Next, I do
# some vector math to calculate the biomass increment (biomass_inc). Then, I filter 
# for the end of season measurement (no increment for first measurement), and again, 
# select and arrange for organization
sub_annual_inc <- bind_rows(nov_2020, jun_2021) %>% 
  arrange(tag) %>% group_by(tag) %>% 
  mutate(biomass_inc = biomass_kg - lag(biomass_kg, default = first(biomass_kg))) %>% 
  filter(row_number() == n()) %>% 
  select(species, tag, biomass_kg, biomass_inc) %>% 
  arrange(tag)

#get rid of negatives
sub_annual_inc$biomass_inc[sub_annual_inc$biomass_inc < 0] <- 0


