#### What the heck am I doing??? The goal is to make some initial figures to see
# what the heckity heck is happening out there
# Will I succeed? Only time will tell
# Probably run "initial_subcanopy_data first so you have all the stuff you need 
# in the global environment

#load necessary packages 
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(plotrix)


# Goal #1: Look at subcanopy NPP over the years w/ variance and compare

#### 2019-2020 first:
sc_2019 <- read.csv("old_data/fd_subcanopy_diameter_2019.csv")

## Setting up the increments
sc_aug_2019 <- filter(sc_2019, date == "2019-08-03" | date == "2019-08-04") %>% 
  select(subplot_id, species, tag, dbh_mm, date)

sc_2020 <- read.csv("old_data/fd_subcanopy_diameter_2020.csv")

sc_nov_2020 <- filter(sc_2020, date == "2020-11-18" | date == "2020-11-19")

sc_19_20 <- bind_rows(sc_aug_2019, sc_nov_2020)

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
  filter(!is.na(dbh_mm)) %>% 
  mutate(dbh_cm = dbh_mm/10) %>%
  mutate(biomass_kg = biomass_a(species, dbh_cm)) 

# create a weeks column based on the date 
sc_19_20_biomass$date <- as.Date(sc_19_20_biomass$date,"%Y-%m-%d")
sc_19_20_biomass$week <- as.Date(cut(sc_19_20_biomass$date, breaks = "week", 
                                     start.on.monday = TRUE))

# now I arrange and group for organization. Next, I do
# some vector math to calculate the biomass increment (biomass_inc). Then, I filter 
# for the end of season measurement (no increment for first measurement), and again, 
# select and arrange for organization. 
sc_20_inc <- sc_19_20_biomass %>% 
  arrange(tag) %>% group_by(tag) %>% 
  mutate(biomass_inc = biomass_kg - lag(biomass_kg, default = first(biomass_kg))) %>% 
  filter(row_number() == n()) %>% 
  select(subplot_id, species, tag, biomass_kg, biomass_inc) %>% 
  arrange(subplot_id) 

# get rid of the negatives in biomass increments. 
# Negatives likely result from inherent error in measuring diameter with 
# digital calipers 
sc_20_inc$biomass_inc[sc_20_inc$biomass_inc < 0] <- 0


# Now I will move from biomass increment to subcanopy NPP. 

# To scale from my sampled population to the entire population, I first 
# calculate and create a new vector with mean subplot biomass increment. This is
# then used for subcanopy species that are within the subplot, but that I did 
# not capture in my sample. I then join the tree counts df which creates a new 
# record within subplots for species that were not captured in my sample population. 

tree_counts <- read.csv("old_data//subcanopy_stemcounts.csv")
# select columns of interest
tree_counts <-  select(tree_counts, subplot, species, count)

# Next, I fill the mean subplot increments, and create a new vector that has 
# either the measured increment (sampled species) or the mean increment 
# (unsampled species). Then, I summarize for mean_biomass for each
# subplot, and species, rejoin tree_counts (lost in the summarize), and scale up to 
# the hectare. 

# Lastly, I sum the biomass production of all species within a subplot
# and multiply by a site specific carbon fraction of 0.48 for kgC_ha_yr
NPP_sc_2020 <- sc_20_inc %>% 
  group_by(subplot_id) %>% 
  summarize(mean_subplot = mean(biomass_inc)) %>% 
  left_join(sc_20_inc) %>% 
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
   rename(NPP_2020 = kgC_ha_yr) %>% 
#  left_join(NPP_sc_2019) %>% 
#  rename(NPP_2019 = NPP_subcan) %>% 
  mutate(replicate = substr(subplot_id, 1, 1)) %>% 
  mutate(severity = case_when(
    subplot_id == "A01E" ~ 0.85, subplot_id == "A01W" ~ 0.85, subplot_id == "A02E" ~ 0.45,
    subplot_id == "A02W" ~ 0.45, subplot_id == "A03E" ~ 0.65, subplot_id == "A03W" ~ 0.65,
    subplot_id == "A04E" ~ 0.00, subplot_id == "A04W" ~ 0.00, subplot_id == "B01E" ~ 0.00,
    subplot_id == "B01W" ~ 0.00, subplot_id == "B02E" ~ 0.45, subplot_id == "B02W" ~ 0.45,
    subplot_id == "B03E" ~ 0.85, subplot_id == "B03W" ~ 0.85, subplot_id == "B04E" ~ 0.65,
    subplot_id == "B04W" ~ 0.65, subplot_id == "C01E" ~ 0.00, subplot_id == "C01W" ~ 0.00,
    subplot_id == "C02E" ~ 0.65, subplot_id == "C02W" ~ 0.65, subplot_id == "C03E" ~ 0.85,
    subplot_id == "C03W" ~ 0.85, subplot_id == "C04E" ~ 0.45, subplot_id == "C04W" ~ 0.45, 
    subplot_id == "D01E" ~ 0.00, subplot_id == "D01W" ~ 0.00, subplot_id == "D02E" ~ 0.85,
    subplot_id == "D02W" ~ 0.85, subplot_id == "D03E" ~ 0.45, subplot_id == "D03W" ~ 0.45,
    subplot_id == "D04E" ~ 0.65, subplot_id == "D04W" ~ 0.65
  )) %>% 
  mutate(treatment = case_when(
    subplot_id == "A01E" ~ "bottom", subplot_id == "A01W" ~ "top", subplot_id == "A02E" ~ "top",
    subplot_id == "A02W" ~ "bottom", subplot_id == "A03E" ~ "bottom", subplot_id == "A03W" ~ "top",
    subplot_id == "A04E" ~ "bottom", subplot_id == "A04W" ~ "top", subplot_id == "B01E" ~ "bottom",
    subplot_id == "B01W" ~ "top", subplot_id == "B02E" ~ "top", subplot_id == "B02W" ~ "bottom",
    subplot_id == "B03E" ~ "bottom", subplot_id == "B03W" ~ "top", subplot_id == "B04E" ~ "top",
    subplot_id == "B04W" ~ "bottom", subplot_id == "C01E" ~ "top", subplot_id == "C01W" ~ "bottom",
    subplot_id == "C02E" ~ "bottom", subplot_id == "C02W" ~ "top", subplot_id == "C03E" ~ "bottom",
    subplot_id == "C03W" ~ "top", subplot_id == "C04E" ~ "top", subplot_id == "C04W" ~ "bottom", 
    subplot_id == "D01E" ~ "bottom", subplot_id == "D01W" ~ "top", subplot_id == "D02E" ~ "bottom",
    subplot_id == "D02W" ~ "top", subplot_id == "D03E" ~ "bottom", subplot_id == "D03W" ~ "top",
    subplot_id == "D04E" ~ "top", subplot_id == "D04W" ~ "bottom"
  )) %>% 
  ungroup()

#################################################################################
# Now we're just gonna repeat for 2021 increment
sc_2021 <- read.csv("2021_data/subcanopy_diameter_2021.csv")
#quick clean 2021
sc_2021 <-  sc_2021 %>% 
  mutate(nested_subplot = substr(uniqueID, 5, 5)) %>% 
  mutate(notes = "NA")

# this deletes the nested subplot ID from the uniqueID column
sc_2021$uniqueID <- gsub('.$', '', sc_2021$uniqueID)

# make all column names lowercase 
names(sc_2021) <- tolower(names(sc_2021))

# now rename columns to match other forte data products
names(sc_2021)[names(sc_2021) == "uniqueid"] <- "subplot_id"
names(sc_2021)[names(sc_2021) == "subplotid"] <- "subplot_id"

# now drop unwanted columns and reorder 
sc_2021 <- sc_2021[c("subplot_id", "nested_subplot", "tag", "species", "dbh_mm", 
                     "date", "notes")]

########

## Setting up the increments
sc_jun_2021 <- filter(sc_2021, date == "2021-06-22" | date == "2021-06-23") %>% 
  select(subplot_id, species, tag, dbh_mm, date)

sc_20_21 <- bind_rows(sc_nov_2020, sc_jun_2021)

# remove records where DBH_mm
# was NA (some stems are missing week 1 and 2 of data collection). Next I convert
# DBH in mm to cm, and use my biomass_a function to estimate biomass(kg)
# (DBH allometry) 

sc_20_21_biomass <- sc_20_21 %>% 
  filter(!is.na(dbh_mm)) %>% 
  mutate(dbh_cm = dbh_mm/10) %>%
  mutate(biomass_kg = biomass_a(species, dbh_cm)) 

# create a weeks column based on the date 
sc_20_21_biomass$date <- as.Date(sc_20_21_biomass$date,"%Y-%m-%d")
sc_20_21_biomass$week <- as.Date(cut(sc_20_21_biomass$date, breaks = "week", 
                                     start.on.monday = TRUE))

# now I arrange and group for organization. Next, I do
# some vector math to calculate the biomass increment (biomass_inc). Then, I filter 
# for the end of season measurement (no increment for first measurement), and again, 
# select and arrange for organization. 
sc_21_inc <- sc_20_21_biomass %>% 
  arrange(tag) %>% group_by(tag) %>% 
  mutate(biomass_inc = biomass_kg - lag(biomass_kg, default = first(biomass_kg))) %>% 
  filter(row_number() == n()) %>% 
  select(subplot_id, species, tag, biomass_kg, biomass_inc) %>% 
  arrange(subplot_id) 

# get rid of the negatives in biomass increments. 
# Negatives likely result from inherent error in measuring diameter with 
# digital calipers 
sc_21_inc$biomass_inc[sc_21_inc$biomass_inc < 0] <- 0


# Now I will move from biomass increment to subcanopy NPP for 2021

NPP_sc_2021 <- sc_21_inc %>% 
  group_by(subplot_id) %>% 
  summarize(mean_subplot = mean(biomass_inc)) %>% 
  left_join(sc_21_inc) %>% 
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
  rename(NPP_2021 = kgC_ha_yr) %>% 
  #  left_join(NPP_sc_2019) %>% 
  #  rename(NPP_2019 = NPP_subcan) %>% 
  mutate(replicate = substr(subplot_id, 1, 1)) %>% 
  mutate(severity = case_when(
    subplot_id == "A01E" ~ 0.85, subplot_id == "A01W" ~ 0.85, subplot_id == "A02E" ~ 0.45,
    subplot_id == "A02W" ~ 0.45, subplot_id == "A03E" ~ 0.65, subplot_id == "A03W" ~ 0.65,
    subplot_id == "A04E" ~ 0.00, subplot_id == "A04W" ~ 0.00, subplot_id == "B01E" ~ 0.00,
    subplot_id == "B01W" ~ 0.00, subplot_id == "B02E" ~ 0.45, subplot_id == "B02W" ~ 0.45,
    subplot_id == "B03E" ~ 0.85, subplot_id == "B03W" ~ 0.85, subplot_id == "B04E" ~ 0.65,
    subplot_id == "B04W" ~ 0.65, subplot_id == "C01E" ~ 0.00, subplot_id == "C01W" ~ 0.00,
    subplot_id == "C02E" ~ 0.65, subplot_id == "C02W" ~ 0.65, subplot_id == "C03E" ~ 0.85,
    subplot_id == "C03W" ~ 0.85, subplot_id == "C04E" ~ 0.45, subplot_id == "C04W" ~ 0.45, 
    subplot_id == "D01E" ~ 0.00, subplot_id == "D01W" ~ 0.00, subplot_id == "D02E" ~ 0.85,
    subplot_id == "D02W" ~ 0.85, subplot_id == "D03E" ~ 0.45, subplot_id == "D03W" ~ 0.45,
    subplot_id == "D04E" ~ 0.65, subplot_id == "D04W" ~ 0.65
  )) %>% 
  mutate(treatment = case_when(
    subplot_id == "A01E" ~ "bottom", subplot_id == "A01W" ~ "top", subplot_id == "A02E" ~ "top",
    subplot_id == "A02W" ~ "bottom", subplot_id == "A03E" ~ "bottom", subplot_id == "A03W" ~ "top",
    subplot_id == "A04E" ~ "bottom", subplot_id == "A04W" ~ "top", subplot_id == "B01E" ~ "bottom",
    subplot_id == "B01W" ~ "top", subplot_id == "B02E" ~ "top", subplot_id == "B02W" ~ "bottom",
    subplot_id == "B03E" ~ "bottom", subplot_id == "B03W" ~ "top", subplot_id == "B04E" ~ "top",
    subplot_id == "B04W" ~ "bottom", subplot_id == "C01E" ~ "top", subplot_id == "C01W" ~ "bottom",
    subplot_id == "C02E" ~ "bottom", subplot_id == "C02W" ~ "top", subplot_id == "C03E" ~ "bottom",
    subplot_id == "C03W" ~ "top", subplot_id == "C04E" ~ "top", subplot_id == "C04W" ~ "bottom", 
    subplot_id == "D01E" ~ "bottom", subplot_id == "D01W" ~ "top", subplot_id == "D02E" ~ "bottom",
    subplot_id == "D02W" ~ "top", subplot_id == "D03E" ~ "bottom", subplot_id == "D03W" ~ "top",
    subplot_id == "D04E" ~ "top", subplot_id == "D04W" ~ "bottom"
  )) %>% 
  ungroup()

# Now I want to compare the 2020 and 2021 increments so I'm gonna combine dataframes
two_yr_NPP <- merge(NPP_sc_2020,NPP_sc_2021,by="subplot_id") %>% 
  select(subplot_id, NPP_2020, replicate.x, severity.x, treatment.x, NPP_2021) %>% 
  rename(replicate = replicate.x) %>% 
  rename(severity = severity.x) %>% 
  rename(treatment = treatment.x)
two_yr_NPP<- two_yr_NPP[c("subplot_id", "NPP_2020", "NPP_2021", "replicate", "severity", 
                     "treatment")]

##### Trying to make some figures?????

