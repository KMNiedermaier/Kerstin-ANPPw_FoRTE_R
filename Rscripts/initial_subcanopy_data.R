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
tree_counts <- read.csv("old_data/subcanopy_stemcounts.csv")
# select columns of interest
tree_counts <-  select(tree_counts, subplot, species, count)

#importing csv of subcanopy diameter measurements 
subcanopy_data_2021 <- read.csv("2021_data/subcanopy_diameter_2021.csv")
subcanopy_data_2020 <- read.csv("old_data/fd_subcanopy_diameter_2020.csv")
subcanopy_data_2019 <- read.csv("old_data/fd_subcanopy_diameter_2019.csv")


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

### Quickly clean up 2021 data

subcanopy_data_2021$nested_subplot <- substr(subcanopy_data_2021$uniqueID, 5, 5)
subcanopy_data_2021$notes <- "NA"

# this deletes the nested subplot ID from the uniqueID column leaving only the subplot_id
subcanopy_data_2021$uniqueID <- gsub('.$', '', subcanopy_data_2021$uniqueID)

# make all column names lowercase 
names(subcanopy_data_2021) <- tolower(names(subcanopy_data_2021))

# now rename columns to match other forte data products
names(subcanopy_data_2021)[names(subcanopy_data_2021) == "uniqueid"] <- "subplot_id"

# now drop unwanted columns and reorder 
subcanopy_data_2021 <- subcanopy_data_2021[c("subplot_id", "nested_subplot", "tag", "species", "dbh_mm", 
                     "date", "notes")]

# select important columns from subcanopy data. Then remove records where DBH_mm
# was NA (some stems are missing week 1 and 2 of data collection). Next I convert
# DBH in mm to cm, and use my biomass_a function to estimate biomass(kg)
# (DBH allometry) 
subcanopy_select <- subcanopy_data_2021 %>% 
  select(subplot_id,  species, tag, dbh_mm, date) %>% 
  filter(!is.na(dbh_mm)) %>% 
  mutate(dbh_cm = dbh_mm/10) %>%
  mutate(biomass_kg = biomass_a(species, dbh_cm)) 

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
  mutate(biomass_kg = biomass_a(species, dbh_cm)) 

# now I bind these two dfs together, arrange and group for organization. Next, I do
# some vector math to calculate the biomass increment (biomass_inc). Then, I filter 
# for the end of season measurement (no increment for first measurement), and again, 
# select and arrange for organization

sub_annual_inc <- bind_rows(nov_2020, jun_2021) %>% 
  arrange(tag) %>% group_by(tag) %>% 
  mutate(biomass_inc = biomass_kg - lag(biomass_kg, default = first(biomass_kg))) %>% 
  filter(row_number() == n()) %>% 
  select(subplot_id, species, tag, biomass_kg, biomass_inc) %>% 
  arrange(subplot_id)

#get rid of negatives
sub_annual_inc$biomass_inc[sub_annual_inc$biomass_inc < 0] <- 0

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
NPP_sc_2021 <- sub_annual_inc %>% 
  group_by(subplot_id) %>% 
  summarize(mean_subplot = mean(biomass_inc)) %>% 
  left_join(sub_annual_inc) %>% 
  right_join(tree_counts) %>% arrange(subplot_id, species) %>% 
  group_by(subplot_id) %>% fill(mean_subplot, .direction = "updown") %>% 
  group_by(subplot_id, species) %>% 
  mutate(biomass_inc_obs_est = case_when(
    !is.na(biomass_inc) ~ biomass_inc,
    is.na(biomass_inc) ~ mean_subplot
  )) %>% 
  summarize(mean_biomass = mean(biomass_inc_obs_est)) %>% 
  right_join(tree_counts) %>% 
  mutate(sp_biomass_ha = mean_biomass * count * 40) %>% 
  # * species count * 40 b/c counts were in 0.025ha
  group_by(subplot_id) %>% 
  summarize(kgC_ha_yr = sum(sp_biomass_ha) * 0.48) %>% 
  rename(NPP_subcan = kgC_ha_yr)




