#set library to dyplr
library(dplyr)
library(tidyr)
library(lubridate)

# importing tree coutn data 
tree_counts <- read.csv("data/subcanopy_stemcounts.csv")

#importing csv of subcanopy diameter measurements 
subcanopy_data <- read.csv("data/subcanopy_D.csv")

# calculates subcanopy stem density
subcan_stem_density <- tree_counts %>% 
  mutate(replicate = substr(subplot, 1, 1)) %>%
  group_by(subplot, replicate) %>% 
  summarize(quart_subplot = sum(count)) %>% 
  mutate(stems_per_ha = quart_subplot*40) %>% #coutns were done in 0.025 ha
  group_by(replicate) %>% 
  summarize(rep_stem_density = mean(stems_per_ha), SE = std.error(stems_per_ha))

# calculates subcanopy species comp
subcan_comp <- tree_counts %>% 
  mutate(replicate = substr(subplot, 1, 1)) %>% 
  select(subplot, replicate, species, count) %>% 
  group_by(replicate, species) %>% 
  summarise(total_stems = sum(count)) #%>% # this is in 0.025 ha plots (1/4 subplots)

subcan_comp <- subcan_comp %>% 
  summarise(total_stems_rep = sum(total_stems)) %>% 
  right_join(subcan_comp) %>% 
  mutate(perc_comp = total_stems/total_stems_rep) 

  

# select wanted columns, group by subplot, and sum total stem counts
tree_counts <- tree_counts %>% 
  select(subplot, species, count) %>% 
  group_by(subplot) %>% 
  summarise(total_stems = sum(count))

# use dplyr to get an increment 
subcanopy_select <- subcanopy_data %>% 
  select(subplot, uniqueID, species, tag, DBH_mm, date) %>% 
  filter(!is.na(DBH_mm)) %>% 
  mutate(DBH_cm = DBH_mm/10) %>% 
  group_by(tag) %>%
  mutate(increment_cm = DBH_cm-lag(DBH_cm, default = first(DBH_cm))) %>% 
  #filter(!is.na(increment_cm)) %>% 
  mutate(RGR_cm = increment_cm / lag(DBH_cm, default = first(DBH_cm)))

#get rid of the negaives
subcanopy_select[subcanopy_select <0] <- 0 

# create a weeks column based on the date 
subcanopy_select$date <- as.Date(subcanopy_select$date,"%Y-%m-%d")
subcanopy_select$week <- as.Date(cut(subcanopy_select$date, breaks = "week", 
                                     start.on.monday = FALSE))

# create a DOY column
subcanopy_select$DOY <- yday(subcanopy_select$date)

# create function for initial biomass based on a's and b's for each species 

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

# simply create accumulating biomass by DBH changes 
# get increment_kg
df_sub <- subcanopy_select %>% 
  mutate(biomass_kg = biomass_a(species, DBH_cm)) %>% 
  group_by(tag) %>% 
  mutate(inc_days = abs(DOY-lag(DOY, default = first(DOY)))) %>% 
  mutate(increment_kg = case_when(
    biomass_kg < lag(biomass_kg, default = first(biomass_kg)) ~ 0,
    biomass_kg >= lag(biomass_kg, default = first(biomass_kg)) ~ biomass_kg-lag(biomass_kg, default = first(biomass_kg))
  )) %>% 
  mutate(RGR_biomass = increment_kg/biomass_kg)

df_summary <- df_sub %>% 
  group_by(subplot, week) %>% 
  summarize(mean_delta_B = mean(increment_kg), inc_days_mean = mean(inc_days)) %>% 
  right_join(tree_counts) %>% 
  mutate(subplot_delta_B = mean_delta_B * total_stems) %>% 
  mutate(B_per_ha = subplot_delta_B * 40) %>% # *40 b/c coutns were 0.025ha plots
  mutate(kgC_per_ha = B_per_ha * 0.48) %>% 
  mutate(kgC_ha_day = kgC_per_ha / inc_days_mean) %>% 
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
  )) 
  
# create df to bind with NPP_total
NPP_subcanopy <- df_sub %>% 
  filter(!is.na(increment_kg)) %>% #removes week 1 measurement w/o an increment
  group_by(subplot, week) %>% 
  summarize(mean_delta_B = mean(increment_kg), inc_days_mean = mean(inc_days)) %>% 
  right_join(tree_counts) %>% 
  mutate(subplot_delta_B = mean_delta_B * total_stems * 4) %>% 
  group_by(subplot) %>% 
  summarise(annual_inc = sum(subplot_delta_B)) %>%
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
  )) %>% 
  select(subplot, NPP_subcan = annual_inc)