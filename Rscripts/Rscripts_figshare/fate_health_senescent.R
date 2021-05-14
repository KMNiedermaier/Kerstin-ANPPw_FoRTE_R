# continue from the grow loop part of NPP script; includes fate as a grouping 
# variable. Fate refers to the live (healthy) and kill (senescent) trees. This script 
# creates df with which to analyse the effect of fate on NPP, creates df's for figures,
# and runs an ANOVA to assess within group significance between live/kill

# install/load packages 
library(agricolae)

# using similar pipeline to NPP_final, creat df with fate as a grouping variable 
NPP_fate <- df %>%
  group_by(tag) %>% # NEW WIH inc_cm
  mutate(biomass_new = case_when(
    species == "AMEL" ~ (a*(DBH_cm_new*10)^b)/1000,
    species != "AMEL" ~ a*DBH_cm_new^b
  )) %>%  
  mutate(biomass_diff = biomass_new - lag(biomass_new, default = first(biomass_new))) %>%
  group_by(subplot, fate) %>% # can ADD SPECIES in here if i want to look at effect of species
  summarise(subplot_biomass_kg = sum(biomass_diff)) %>%
  mutate(biomass_per_ha = subplot_biomass_kg*10) %>%
  mutate(kgC_ha_yr = biomass_per_ha*0.48) %>%
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
    subplot == "D01E" ~ "bottom", subplot == "D01W" ~ "top", subplot == "D02E" ~ "top",
    subplot == "D02W" ~ "bottom", subplot == "D03E" ~ "bottom", subplot == "D03W" ~ "top",
    subplot == "D04E" ~ "bottom", subplot == "D04W" ~ "top"
  )) %>% 
  group_by(subplot) %>% 
  mutate(total = kgC_ha_yr + lag(kgC_ha_yr, default = last(kgC_ha_yr))) %>% #subplot total kgC_ha_yr
  mutate(perc_total = kgC_ha_yr/total) %>% 
  ungroup()

# make severity %'s for x axis
NPP_fate$severity <- recode(NPP_fate$severity, "0.00" = "0",
                              "0.45" = "45", "0.65" = "65",
                              "0.85" = "85")

# summarise NPP_fate by severity and prepare a df for plots 
NPP_fate_mean <- NPP_fate %>%
  group_by(severity, fate) %>% 
  summarise(NPP_mean = mean(kgC_ha_yr), SE = std.error(kgC_ha_yr)) %>% 
  group_by(severity) %>% 
  mutate(errorbars = case_when(
    fate == "live" ~ NPP_mean,
    fate == "kill" ~ sum(NPP_mean)
  )) %>% 
  ungroup()

#run some stats on fate contributions to anppw
NPP_fate$severity <- as.factor(NPP_fate$severity)
fate_stats_sev <- aov(perc_total ~ fate + severity + fate*severity, data = NPP_fate)
summary(fate_stats_sev)
TukeyHSD(fate_stats_sev)

# run stats on disturbance type differences 
fate_stats_tx<- aov(perc_total ~ fate + treatment + fate*treatment, data = NPP_fate)
summary(fate_stats_tx)
TukeyHSD(fate_stats_tx)

# df for treatment split by live-kill for figures; EXCLUDES CONTROL TREES
NPP_fate_tx <- NPP_fate %>%
  filter(severity != 0) %>% 
  group_by(treatment, fate) %>% 
  summarise(NPP_mean = mean(kgC_ha_yr), SE = std.error(kgC_ha_yr)) %>% 
  group_by(treatment) %>% 
  mutate(errorbars = case_when(
    fate == "live" ~ NPP_mean,
    fate == "kill" ~ sum(NPP_mean)
  )) 

NPP_fate_tx$treatment <- recode(NPP_fate_tx$treatment, bottom = "Bottom-Up", top = "Top-Down")

