#### continue from NPP script; this script is designed to breakdown NPP by early and 
#### late successional plant functional type (PFT)

# use similar pipe as NPP_final, but break down by PFT as well
# this df includes week as a grouping variable and thus is best suited for a time series 
NPP_PFT <- df %>% 
  group_by(tag) %>% # NEW WIH inc_cm
  mutate(PFT = case_when(
    species == "ACPE" ~ "late",
    species == "ACRU" ~ "late",
    species == "ACSA" ~ "late",
    species == "AMEL" ~ "late",
    species == "BEPA" ~ "early",
    species == "FAGR" ~ "late",
    species == "PIRE" ~ "late",
    species == "PIST" ~ "late",
    species == "POGR" ~ "early",
    species == "POTR" ~ "late",
    species == "QURU" ~ "late",
    species == "ABBA" ~ "late",
    species == "TSCA" ~ "late"
  )) %>% 
  mutate(biomass_new = case_when(
    species == "AMEL" ~ (a*(DBH_cm_new*10)^b)/1000,
    species != "AMEL" ~ a*DBH_cm_new^b)) %>%  
  mutate(biomass_diff = biomass_new - lag(biomass_new, default = first(biomass_new)))%>% 
  group_by(subplot, week, PFT) %>% 
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
    subplot == "D01E" ~ "bottom", subplot == "D01W" ~ "top", subplot == "D02E" ~ "bottom",
    subplot == "D02W" ~ "top", subplot == "D03E" ~ "bottom", subplot == "D03W" ~ "top",
    subplot == "D04E" ~ "top", subplot == "D04W" ~ "bottom"
  )) 

# getting DOY into dates and correcting the november 2018 date 
dates <- data.frame(as.factor(as.Date(NPP_PFT$DOY, "2019-01-01")))
names(dates) <- c("date")
dates <- data.frame(recode_factor(dates$date, '2019-11-03' = "2018-11-15"))
names(dates) <- c("date")
dates$date <- as.Date(dates$date,"%Y-%m-%d")

NPP_PFT <- NPP_PFT %>% 
  bind_cols(dates) 

#create week and month columns to look at different temporal scales 
NPP_PFT$month <- as.Date(cut(NPP_PFT$date, breaks = "month"))
NPP_PFT$week <- as.Date(cut(NPP_PFT$date, breaks = "week", start.on.monday = FALSE))

##############################################################################
######## PFT Annual ANPP; this does not include week as a grouping variable and is 
######## thus meant for visualizing and analyzing annual (cumulative) NPP among PFT's

# continue from NPP script; creating df that has annual NPP by PFT
annual_PFT <- df %>% 
  group_by(tag) %>% # NEW WIH inc_cm
  mutate(PFT = case_when(
    species == "ACPE" ~ "late",
    species == "ACRU" ~ "late",
    species == "ACSA" ~ "late",
    species == "AMEL" ~ "late",
    species == "BEPA" ~ "early",
    species == "FAGR" ~ "late",
    species == "PIRE" ~ "late",
    species == "PIST" ~ "late",
    species == "POGR" ~ "early",
    species == "POTR" ~ "late",
    species == "QURU" ~ "late",
    species == "ABBA" ~ "late",
    species == "TSCA" ~ "late"
  )) %>% 
  mutate(biomass_new = case_when(
    species == "AMEL" ~ (a*(DBH_cm_new*10)^b)/1000,
    species != "AMEL" ~ a*DBH_cm_new^b
  )) %>%  # NEW WITH inc_cm  # NEW WITH inc_cm
  mutate(biomass_diff = biomass_new - lag(biomass_new, default = first(biomass_new))) %>% 
  group_by(subplot, PFT) %>% 
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
  group_by(severity, PFT) %>% 
  summarise(NPP_canopy = mean(kgC_ha_yr), SE = std.error(kgC_ha_yr)) %>% 
  group_by(severity) %>% 
  mutate(errorbars = case_when(
    PFT == "late" ~ NPP_canopy,
    PFT == "early" ~ sum(NPP_canopy)
  ))

# make severity %'s for x axis
annual_PFT$severity <- recode(annual_PFT$severity, "0.00" = "0",
                              "0.45" = "45", "0.65" = "65",
                              "0.85" = "85")
