### NPP 2021
# First attempted 08/01/2021
# Kerstin Niedermaier, using Maxim Grigri's NPP_2020 script

# This script is the pipeline from dendrometer bands measurements of ~700 trees
# across the entire FoRTE experiment to ANPPw. First, it runs the subcanopy, and 
# seedling and sapling scripts using the source function. Then it begins the pipeline
# To add new data to the pipeline:
#     1. Import new .csv file (see dendro_2019, dendro_2020)
#     2. Ensure dates are in YYYY-MM-DD format
#     3. Add new df to the bind_rows function
# After ANPPw is calculated, it summarize data and prepares df's for ggplot2 in the 
# "figures.R" script, as well as, the "stats.R" script. 

# load packages I use 
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotrix)
library(gridExtra)
library(grid)
library(ggpubr)
library(lubridate)
library(tidyverse)

# use source to run subcanopy scripts 
# source("Rscripts/subcanopy.R")
# source("Rscripts/seedling_saplings.R")

#importing csv *********ADD NEW DATA HERE
dendro_2019 <- read.csv("old_data/canopy_dendrobands_2019.csv", na.strings = c("", "NA"))
dendro_2020 <- read.csv("old_data/canopy_dendrobands_2020.csv", na.strings = c("", "NA"))
dendro_2021 <- read.csv("2021_data/canopy_dendrobands_2021.csv", na.strings = c("", "NA"))

#change name of bands column for 2021 to match others and a week column
dendro_2021 <- dendro_2021 %>% 
  rename(
    bands_in = bottom_in) %>% 
  mutate(week = 20)


#combine 2019 and 2020 and 2021
dendro_data <- bind_rows(dendro_2019, dendro_2020, dendro_2021)

source("Rscripts/pipeline.R")

#####################################################################################
#This chunk scales D growth to daily NPP (kgC/ha/day) and annual NPP (kgC/ha/yr)

#create a df that has the # of days (inc_days) between measurements to join with 
#summarized verion of df below
inc_days <- select(df, subplot, week, inc_days) %>% unique()

# filter for Nov. 2019 and Nov. 2018, convert to biomass using allometric parameters,
# sum biomass for each subplot and date, find difference between the two measurements,
# scale to area and C. Lastly, manually code in severity and type. 
NPP_2020 <- df %>% 
  filter(week == "2019-11-10" | week == "2020-11-15") %>% 
  mutate(biomass_new = case_when(
    species == "AMEL" ~ (a*(DBH_cm_new*10)^b)/1000,
    species != "AMEL" ~ a*DBH_cm_new^b
  )) %>%
  group_by(subplot, date, week) %>% 
  summarise(subplot_biomass_kg = sum(biomass_new)) %>%
  group_by(subplot) %>% 
  mutate(biomass_diff = subplot_biomass_kg - lag(subplot_biomass_kg, 
                                                 default = first(subplot_biomass_kg))) %>% 
  filter(week == "2020-11-15") %>% 
  mutate(biomass_per_ha = biomass_diff*10) %>% 
  mutate(NPP_2020 = biomass_per_ha*0.48) %>% # kgC/ha/yr
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
  select(subplot, NPP_2020, replicate, severity, treatment)

# make severity %'s for x axis
# NPP_2020$severity <- recode(NPP_2020$severity, "0.00" = "0",
#                              "0.45" = "45", "0.65" = "65",
#                              "0.85" = "85")

# same as NPP 2020 except change the dates 
NPP_2019 <- df %>% 
  filter(week == "2018-10-28" | week == "2019-11-10") %>% 
  mutate(biomass_new = case_when(
    species == "AMEL" ~ (a*(DBH_cm_new*10)^b)/1000,
    species != "AMEL" ~ a*DBH_cm_new^b
  )) %>%
  group_by(subplot, date, week) %>% 
  summarise(subplot_biomass_kg = sum(biomass_new)) %>%
  group_by(subplot) %>% 
  mutate(biomass_diff = subplot_biomass_kg - lag(subplot_biomass_kg, 
                                                 default = first(subplot_biomass_kg))) %>% 
  filter(week == "2019-11-10") %>% 
  mutate(biomass_per_ha = biomass_diff*10) %>% 
  mutate(NPP_2019 = biomass_per_ha*0.48) %>% # kgC/ha/yr
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
  select(subplot, NPP_2019, replicate, severity, treatment) %>% 
  ungroup()

# repeat for 2021 (with different dates)
NPP_2021 <- df %>% 
  filter(week == "2020-11-15" | week == "2021-06-13") %>% 
  mutate(biomass_new = case_when(
    species == "AMEL" ~ (a*(DBH_cm_new*10)^b)/1000,
    species != "AMEL" ~ a*DBH_cm_new^b
  )) %>%
  group_by(subplot, date, week) %>% 
  summarise(subplot_biomass_kg = sum(biomass_new)) %>%
  group_by(subplot) %>% 
  mutate(biomass_diff = subplot_biomass_kg - lag(subplot_biomass_kg, 
                                                 default = first(subplot_biomass_kg))) %>% 
  filter(week == "2021-06-13") %>% 
  mutate(biomass_per_ha = biomass_diff*10) %>% 
  mutate(NPP_2021 = biomass_per_ha*0.48) %>% # kgC/ha/yr
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
  select(subplot, NPP_2021, replicate, severity, treatment) %>% 
  ungroup()

# join 2019 and 2020 and 2021
NPP_alltime <- left_join(NPP_2021,NPP_2019, NPP_2020)

# reformat for plots
NPP_alltime <- NPP_alltime %>% 
  gather(year, kgC_ha_yr, NPP_2019, NPP_2020, NPP_2021) #%>% 
#recode year column to "2019" and "2020"
NPP_alltime$year <- recode(NPP_alltime$year, NPP_2019 = "2019", NPP_2020 = "2020")

# make severity %'s for x axis
NPP_alltime$severity <- recode(NPP_alltime$severity, "0.00" = "0",
                               "0.45" = "45", "0.65" = "65",
                               "0.85" = "85")

# plot the 2020 canopy NPP
ggplot(NPP_alltime, aes(x = factor(severity), y = kgC_ha_yr, fill = factor(year))) +
  geom_boxplot() +
  theme_bw()+
  theme(axis.text.x = element_text(size = 17), axis.text.y = 
          element_text(size = 17), axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17), legend.text = element_text(size = 15),
        legend.title = element_blank()) +
  labs(x = "Disturbance Severity (%)",
       y = expression(paste("ANPP" [w], " ( ",kgC," ",ha^-1," ",day^-1,")")))

# plot by disturbance type
ggplot(NPP_alltime, aes(x = factor(treatment), y = kgC_ha_yr, fill = factor(year))) +
  geom_boxplot() +
  theme_bw()+
  theme(axis.text.x = element_text(size = 17), axis.text.y = 
          element_text(size = 17), axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17), legend.text = element_text(size = 15),
        legend.title = element_blank()) +
  labs(x = "Disturbance Type",
       y = expression(paste("ANPP" [w], " ( ",kgC," ",ha^-1," ",day^-1,")")))

###################################################################################
# break out NPP by fate (Live/Kill)

# starting from df, use same pipe from NPP_2020 but add fate and grouping variable
NPP_2020_fate <- df %>% 
  filter(week == "2019-11-10" | week == "2020-11-15") %>% 
  mutate(biomass_new = case_when(
    species == "AMEL" ~ (a*(DBH_cm_new*10)^b)/1000,
    species != "AMEL" ~ a*DBH_cm_new^b
  )) %>%
  group_by(subplot, date, week, fate) %>% 
  summarise(subplot_biomass_kg = sum(biomass_new)) %>%
  group_by(subplot, fate) %>% 
  mutate(biomass_diff = subplot_biomass_kg - lag(subplot_biomass_kg, 
                                                 default = first(subplot_biomass_kg))) %>% 
  filter(week == "2020-11-15") %>% 
  mutate(biomass_per_ha = biomass_diff*10) %>% 
  mutate(NPP_2020 = biomass_per_ha*0.48) %>% # kgC/ha/yr
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
  select(subplot, NPP_2020, fate, replicate, severity, treatment) %>% 
  ungroup() 

# before summarizing NPP_fate, I need to add 8 cases (rows) to fit the 0% severity-
# kill category (no kill trees on control plots). I Will do this by adding 8 rows 
# to NPP_fate and manually entering the subplot, fate, and NPP
# (0 NPP because there are no kill trees) nd then filling severity, treatmnet (type),
# rep
NPP_2020_fate <- NPP_2020_fate %>% 
  add_row(subplot = "A04E", fate = "kill", NPP_2020 = 0) %>% 
  add_row(subplot = "A04W", fate = "kill", NPP_2020 = 0) %>% 
  add_row(subplot = "B01E", fate = "kill", NPP_2020 = 0) %>% 
  add_row(subplot = "B01W", fate = "kill", NPP_2020 = 0) %>% 
  add_row(subplot = "C01E", fate = "kill", NPP_2020 = 0) %>% 
  add_row(subplot = "C01W", fate = "kill", NPP_2020 = 0) %>% 
  add_row(subplot = "D01E", fate = "kill", NPP_2020 = 0) %>% 
  add_row(subplot = "D01W", fate = "kill", NPP_2020 = 0) %>% 
  group_by(subplot) %>% fill(severity, .direction = "updown") %>% 
  fill(treatment, .direction = "updown") %>% 
  fill(replicate, .direction = "updown") %>% 
  arrange(subplot)

# same as NPP 2020_fate except change the dates 
NPP_2019_fate <- df %>% 
  filter(week == "2018-10-28" | week == "2019-11-10") %>% 
  mutate(biomass_new = case_when(
    species == "AMEL" ~ (a*(DBH_cm_new*10)^b)/1000,
    species != "AMEL" ~ a*DBH_cm_new^b
  )) %>%
  group_by(subplot, date, week, fate) %>% 
  summarise(subplot_biomass_kg = sum(biomass_new)) %>%
  group_by(subplot, fate) %>% 
  mutate(biomass_diff = subplot_biomass_kg - lag(subplot_biomass_kg, 
                                                 default = first(subplot_biomass_kg))) %>% 
  filter(week == "2019-11-10") %>% 
  mutate(biomass_per_ha = biomass_diff*10) %>% 
  mutate(NPP_2019 = biomass_per_ha*0.48) %>% # kgC/ha/yr
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
  select(subplot, NPP_2019, fate, replicate, severity, treatment) %>% 
  ungroup() 

# before summarizing MPP_fate, I need to add 8 cases (rows) to fit the 0% severity-
# kill category (no kill trees on control plots). I Will do this by adding 8 rows 
# to NPP_fate and manually entering the subplot, fate, and NPP
# (0 NPP because there are no kill trees) nd then filling severity, treatmnet (type),
# rep
NPP_2019_fate <- NPP_2019_fate %>% 
  add_row(subplot = "A04E", fate = "kill", NPP_2019 = 0) %>% 
  add_row(subplot = "A04W", fate = "kill", NPP_2019 = 0) %>% 
  add_row(subplot = "B01E", fate = "kill", NPP_2019 = 0) %>% 
  add_row(subplot = "B01W", fate = "kill", NPP_2019 = 0) %>% 
  add_row(subplot = "C01E", fate = "kill", NPP_2019 = 0) %>% 
  add_row(subplot = "C01W", fate = "kill", NPP_2019 = 0) %>% 
  add_row(subplot = "D01E", fate = "kill", NPP_2019 = 0) %>% 
  add_row(subplot = "D01W", fate = "kill", NPP_2019 = 0) %>% 
  group_by(subplot) %>% fill(severity, .direction = "updown") %>% 
  fill(treatment, .direction = "updown") %>% 
  fill(replicate, .direction = "updown") %>% 
  arrange(subplot)

# join 2019 and 2020 
NPP_alltime_fate <- left_join(NPP_2019_fate, NPP_2020_fate)

# reformat for plots with a year column
NPP_alltime_fate <- NPP_alltime_fate %>% 
  gather(year, kgC_ha_yr, NPP_2019, NPP_2020) #%>% 
#recode year column to "2019" and "2020"
NPP_alltime_fate$year <- recode(NPP_alltime_fate$year, NPP_2019 = "2019", NPP_2020 = "2020")

# make severity %'s for x axis
NPP_alltime_fate$severity <- recode(NPP_alltime_fate$severity, "0.00" = "0",
                                    "0.45" = "45", "0.65" = "65",
                                    "0.85" = "85")
# create two new columns in NPP_alltime_fate: 1) with total ANPP, and 2) with % of total
# ANPP. there are 54% more kill trees than live trees, thus I adjust kill tree totals
# by 54%
NPP_alltime_fate <- NPP_alltime_fate %>% 
  group_by(subplot, treatment, year) %>% 
  mutate(subplot_total = kgC_ha_yr + lag(kgC_ha_yr)) %>% 
  fill(subplot_total, .direction = "up") %>%
  mutate(perc_total = (kgC_ha_yr / subplot_total)*100) #%>% 
# #mutate(perc_total_adj = case_when(
#   fate == "kill" ~ perc_total * 0.54,
#   fate == "live" ~ perc_total
# ))

# plot the 2020 canopy NPP
ggplot(NPP_alltime_fate, aes(x = factor(severity), y = kgC_ha_yr, fill = factor(fate))) +
  geom_boxplot() +
  theme_bw()+
  facet_wrap(factor(year) ~.)+
  theme(axis.text.x = element_text(size = 17), axis.text.y = 
          element_text(size = 17), axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17), legend.text = element_text(size = 15),
        legend.title = element_blank()) +
  labs(x = "Disturbance Severity (%)",
       y = expression(paste("ANPP" [w], " ( ",kgC," ",ha^-1," ",day^-1,")")))

# plot by disturbance type
ggplot(NPP_alltime_fate, aes(x = factor(treatment), y = kgC_ha_yr, fill = factor(fate))) +
  geom_boxplot() +
  theme_bw()+
  facet_wrap(factor(year)~.)+
  theme(axis.text.x = element_text(size = 17), axis.text.y = 
          element_text(size = 17), axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17), legend.text = element_text(size = 15),
        legend.title = element_blank()) +
  labs(x = "Disturbance Type",
       y = expression(paste("ANPP" [w], " ( ",kgC," ",ha^-1," ",day^-1,")")))

# plot year:fate 
ggplot(NPP_alltime_fate, aes(x = factor(year), y = kgC_ha_yr, fill = factor(fate))) +
  geom_boxplot() +
  theme_bw()+
  theme(axis.text.x = element_text(size = 17), axis.text.y = 
          element_text(size = 17), axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17), legend.text = element_text(size = 15),
        legend.title = element_blank()) +
  labs(x = "fate",
       y = expression(paste("ANPP" [w], " ( ",kgC," ",ha^-1," ",day^-1,")")))

# make a df with the text I want to annotate 
ann_text <- data.frame(severity = c("45", "65", "85", "65", "85"), 
                       treatment = c("bottom"),
                       lab = "*", 
                       fate = "kill",
                       perc_total = c(0), 
                       year = factor(c(2020, 2020, 2020, 2019, 2019),
                                     levels = c(2019,2020)))

#remove control plots when presenting experiment wide live-kill data 
no_control <- filter(NPP_alltime_fate, severity != "0")

# plot live-kill by severity with % total anpp as the y
ggplot(no_control, aes(x = factor(severity), y = perc_total, fill = factor(fate))) +
  geom_boxplot() +
  theme_bw()+
  facet_wrap(factor(year) ~.)+
  theme(axis.text.x = element_text(size = 17), axis.text.y = 
          element_text(size = 17), axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17), legend.text = element_text(size = 15),
        legend.title = element_blank()) +
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100), limits = c(0,110))+
  labs(x = "Disturbance Severity (%)",
       y = expression(paste(" % of total wood NPP")))+
  geom_text(data = ann_text, label = "*", size = 10, y = 100)+
  scale_fill_manual(values = c("orange3","palegreen3"), labels = c("Dying", "Healthy"))


# plot year:fate with perc_total on the y
# first take out the control plots 
perc_fate <- filter(NPP_alltime_fate, severity != "0")

# plot percent contribution to total ANPPw from live and dead trees agreggated by 
# live and dead trees experiment wide 
ggplot(perc_fate, aes(x = factor(year), y = perc_total, fill = factor(fate))) +
  geom_boxplot() +
  theme_bw()+
  theme(axis.text.x = element_text(size = 17), axis.text.y = 
          element_text(size = 17), axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17), legend.text = element_text(size = 15),
        legend.title = element_blank()) +
  labs(x = "fate",
       y = expression(paste("Adj. % of total ANPP" [w], " ( ",kgC," ",ha^-1," ",day^-1,")")))

# create a summary table for no_control
perc_fate %>% 
  group_by(year, fate, severity) %>% 
  #filter(severity == "45" & year == "2020" & fate == "kill") %>% 
  summarise(mean(perc_total
  ))

