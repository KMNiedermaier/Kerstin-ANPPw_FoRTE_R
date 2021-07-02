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
library(tidyverse)

# importing tree count data 
tree_counts <- read.csv("old_data/subcanopy_stemcounts.csv")
# select columns of interest
tree_counts <-  select(tree_counts, subplot, species, count)

# Goal #1: Look at subcanopy NPP over the years w/ variance and compare

#### 2019 first:
sc_2019 <- read.csv("old_data/fd_subcanopy_diameter_2019.csv")


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
subcanopy_select_2019 <- sc_2019 %>% 
  select(subplot_id,  species, tag, dbh_mm, date) %>% 
  filter(!is.na(dbh_mm)) %>% 
  mutate(dbh_cm = dbh_mm/10) %>%
  mutate(biomass_kg = biomass_a(species, dbh_cm)) 

# create a weeks column based on the date 
subcanopy_select_2019$date <- as.Date(subcanopy_select_2019$date,"%Y-%m-%d")
subcanopy_select_2019$week <- as.Date(cut(subcanopy_select_2019$date, breaks = "week", 
                                     start.on.monday = TRUE))

# create a DOY column
subcanopy_select_2019$DOY <- yday(subcanopy_select_2019$date)

# filter subcanopy select for for first and last DBH measurement for every sample
# here is the first measurement of every sampled tree (2019)
may_2019 <- subcanopy_select_2019 %>% 
  group_by(tag) %>% 
  filter(row_number() == 1) %>% 
  mutate(biomass_kg = biomass_a(species, dbh_cm)) 

# here is the final measurement of every sampled tree (2019)
aug_2019 <- subcanopy_select_2019 %>% 
  group_by(tag) %>% 
  filter(row_number() == n()) %>% 
  mutate(biomass_kg = biomass_a(species, dbh_cm)) 


# now I bind these two dfs together, arrange and group for organization. Next, I do
# some vector math to calculat the biomass increment (biomass_inc). Then, I filter 
# for the end of season measuremnt (no increment for first measurement), and again, 
# select and arrange for organization
sub_annual_inc_2019 <- bind_rows(may_2019, aug_2019) %>% 
  arrange(tag) %>% group_by(tag) %>% 
  mutate(biomass_inc = biomass_kg - lag(biomass_kg, default = first(biomass_kg))) %>% 
  filter(row_number() == n()) %>% 
  select(subplot_id, species, tag, biomass_kg, biomass_inc) %>% 
  arrange(subplot_id)


# get rid of the negatives in biomass increments. 
# Negatives likely result from inherent error in measuring diameter with 
# digital calipers 
sub_annual_inc_2019$biomass_inc[sub_annual_inc_2019$biomass_inc < 0] <- 0


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
NPP_sc_2019 <- sub_annual_inc_2019 %>% 
  group_by(subplot_id) %>% 
  summarize(mean_subplot = mean(biomass_inc)) %>% 
  left_join(sub_annual_inc_2019) %>% 
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
  rename(NPP_subcan = kgC_ha_yr) %>% 
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
  mutate(year = 2019)

### Now 2020 data
sc_2020 <- read.csv("old_data/fd_subcanopy_diameter_2020.csv")

# select important columns from subcanopy data. Then remove records where DBH_mm
# was NA (some stems are missing week 1 and 2 of data collection). Next I convert
# DBH in mm to cm, and use my biomass_a function to estimate biomass(kg)
# (DBH allometry) 
subcanopy_select_2020 <- sc_2020 %>% 
  select(subplot_id,  species, tag, dbh_mm, date) %>% 
  filter(!is.na(dbh_mm)) %>% 
  mutate(dbh_cm = dbh_mm/10) %>%
  mutate(biomass_kg = biomass_a(species, dbh_cm)) 


# create a weeks column based on the date 
subcanopy_select_2020$date <- as.Date(subcanopy_select_2020$date,"%Y-%m-%d")
subcanopy_select_2020$week <- as.Date(cut(subcanopy_select_2020$date, breaks = "week", 
                                          start.on.monday = TRUE))

# create a DOY column
subcanopy_select_2020$DOY <- yday(subcanopy_select_2020$date)

# filter subcanopy select for for first and last DBH measurement for every sample
# here is the first measurement of every sampled tree (2020)
jul_2020 <- subcanopy_select_2020 %>% 
  group_by(tag) %>% 
  filter(row_number() == 1) %>% 
  mutate(biomass_kg = biomass_a(species, dbh_cm)) 

# here is the final measurement of every sampled tree (2020)
nov_2020 <- subcanopy_select_2020 %>% 
  group_by(tag) %>% 
  filter(row_number() == n()) %>% 
  mutate(biomass_kg = biomass_a(species, dbh_cm)) 

# now I bind these two dfs together, arrange and group for organization. Next, I do
# some vector math to calculat the biomass increment (biomass_inc). Then, I filter 
# for the end of season measuremnt (no increment for first measurement), and again, 
# select and arrange for organization
sub_annual_inc_2020 <- bind_rows(aug_2019, nov_2020) %>% 
  arrange(tag) %>% group_by(tag) %>% 
  mutate(biomass_inc = biomass_kg - lag(biomass_kg, default = first(biomass_kg))) %>% 
  filter(row_number() == n()) %>% 
  select(subplot_id, species, tag, biomass_kg, biomass_inc) %>% 
  arrange(subplot_id)

# get rid of the negatives in biomass increments. 
# Negatives likely result from inherent error in measuring diameter with 
# digital calipers 
sub_annual_inc_2020$biomass_inc[sub_annual_inc_2020$biomass_inc < 0] <- 0


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
NPP_sc_2020 <- sub_annual_inc_2020 %>% 
  group_by(subplot_id) %>% 
  summarize(mean_subplot = mean(biomass_inc)) %>% 
  left_join(sub_annual_inc_2020) %>% 
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
  rename(NPP_subcan = kgC_ha_yr) %>% 
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
  mutate(year = 2020)


####### Now 2021 Data!
sc_2021 <- read.csv("2021_data/subcanopy_diameter_2021.csv")

### Quickly clean up 2021 data

sc_2021$nested_subplot <- substr(sc_2021$uniqueID, 5, 5)
sc_2021$notes <- "NA"

# this deletes the nested subplot ID from the uniqueID column leaving only the subplot_id
sc_2021$uniqueID <- gsub('.$', '', sc_2021$uniqueID)

# make all column names lowercase 
names(sc_2021) <- tolower(names(sc_2021))

# now rename columns to match other forte data products
names(sc_2021)[names(sc_2021) == "uniqueid"] <- "subplot_id"

# now drop unwanted columns and reorder 
sc_2021 <- sc_2021[c("subplot_id", "nested_subplot", "tag", "species", "dbh_mm", 
                                             "date", "notes")]

# select important columns from subcanopy data. Then remove records where DBH_mm
# was NA (some stems are missing week 1 and 2 of data collection). Next I convert
# DBH in mm to cm, and use my biomass_a function to estimate biomass(kg)
# (DBH allometry) 
subcanopy_select_2021 <- sc_2021 %>% 
  select(subplot_id,  species, tag, dbh_mm, date) %>% 
  filter(!is.na(dbh_mm)) %>% 
  mutate(dbh_cm = dbh_mm/10) %>%
  mutate(biomass_kg = biomass_a(species, dbh_cm)) 


# create a weeks column based on the date 
subcanopy_select_2021$date <- as.Date(subcanopy_select_2021$date,"%Y-%m-%d")
subcanopy_select_2021$week <- as.Date(cut(subcanopy_select_2021$date, breaks = "week", 
                                          start.on.monday = TRUE))

# create a DOY column
subcanopy_select_2021$DOY <- yday(subcanopy_select_2021$date)

# filter subcanopy select for for first and last DBH measurement for every sample
# here is the first measurement of every sampled tree (2021)
jun_2021 <- subcanopy_select_2021 %>% 
  group_by(tag) %>% 
  filter(row_number() == 1) %>% 
  mutate(biomass_kg = biomass_a(species, dbh_cm)) 


# now I bind these two dfs together, arrange and group for organization. Next, I do
# some vector math to calculat the biomass increment (biomass_inc). Then, I filter 
# for the end of season measuremnt (no increment for first measurement), and again, 
# select and arrange for organization
sub_annual_inc_2021 <- bind_rows(nov_2020, jun_2021) %>% 
  arrange(tag) %>% group_by(tag) %>% 
  mutate(biomass_inc = biomass_kg - lag(biomass_kg, default = first(biomass_kg))) %>% 
  filter(row_number() == n()) %>% 
  select(subplot_id, species, tag, biomass_kg, biomass_inc) %>% 
  arrange(subplot_id)

# get rid of the negatives in biomass increments. 
# Negatives likely result from inherent error in measuring diameter with 
# digital calipers 
sub_annual_inc_2021$biomass_inc[sub_annual_inc_2021$biomass_inc < 0] <- 0


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
NPP_sc_2021 <- sub_annual_inc_2021 %>% 
  group_by(subplot_id) %>% 
  summarize(mean_subplot = mean(biomass_inc)) %>% 
  left_join(sub_annual_inc_2021) %>% 
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
  rename(NPP_subcan = kgC_ha_yr) %>% 
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
  mutate(year = 2021)


#############################

# Now I want to compare the 2020 and 2021 increments so I'm gonna combine dataframes
three_yr_NPP <- rbind(NPP_sc_2019, NPP_sc_2020, NPP_sc_2021)

# make sure the year variable is not numerical
three_yr_NPP$year <- as.factor(three_yr_NPP$year)

##### Goal #2: Trying to make some figures?????

# Summary stats
sc_NPP_summary <- three_yr_NPP %>%
  group_by(year) %>%
  summarize(
    mean_sc_NPP = mean(NPP_subcan, na.rm = TRUE),
    median_sc_NPP = median(NPP_subcan, na.rm = TRUE),
    variance_sc_NPP = var(NPP_subcan, na.rm = TRUE),
    sd_sc_NPP = sd(NPP_subcan, na.rm = TRUE),
    n = n(),
    se_NPP = sd_sc_NPP / sqrt(n))
    
    
sc_NPP_bar <-
  ggplot(sc_NPP_summary, aes(x = year, y = mean_sc_NPP)) +
  geom_col() +
  labs(y = "Mean Subcanopy NPP") +
  theme_classic() +
  geom_errorbar(aes(
    ymax = mean_sc_NPP + se_NPP, 
    ymin = mean_sc_NPP - se_NPP),
    width = 0.5) +
  scale_y_continuous(expand = c(0,0))

sc_NPP_bar

####### Using Max's figures to make my own
three_yr_NPP$severity <- recode(three_yr_NPP$severity, "0.00 = 0;
                                  0.45 = 45; 0.65 = 65;
                                  0.85 = 85")

# Plottin'
ggplot(three_yr_NPP, aes(x = factor(severity), y = NPP_subcan, fill = factor(year))) +
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x = element_text(size = 17), axis.text.y = 
          element_text(size = 17), axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17), legend.text = element_text(size = 15),
        legend.title = element_blank()) +
  labs(x = "Disturbance Severity (%)",
       y = expression(paste("ANPP" [w], " ( ",kgC," ",ha^-1," ",day^-1,")"))) 

ggplot(three_yr_NPP, aes(x = factor(treatment), y = NPP_subcan, fill = factor(year))) +
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x = element_text(size = 17), axis.text.y = 
          element_text(size = 17), axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17), legend.text = element_text(size = 15),
        legend.title = element_blank()) +
  labs(x = "Disturbance Type",
       y = expression(paste("ANPP" [w], " ( ",kgC," ",ha^-1," ",day^-1,")")))

#### ANOVA
library(car)

attach(three_yr_NPP)
hist(NPP_subcan)
qqnorm(NPP_subcan)
qqline(NPP_subcan)
shapiro.test(NPP_subcan)
# abnormal data

leveneTest(NPP_subcan_lm)
# equal variance

NPP_subcan_lm <- lm(NPP_subcan ~ year)
anova(NPP_subcan_lm)
# no significant differences between years

################################################################################
# playing around with forestr

if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("FoRTExperiment/fortedata", dependencies = TRUE, build_vignettes = FALSE)

library(fortedata)

fd_canopy_structure_summary()
fd_dendro()
fd_subcanopy_diameter()
