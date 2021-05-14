# To calc stems/subplot: 
# # select wanted columns, group by subplot, and sum total stem counts
# stems_subplot <- tree_counts %>% 
#   rename(total_stems = count) %>% 
#   slice(rep(1:n(), each = 7)) %>% 
#   group_by(subplot, species) %>% 
#   mutate(week_num = row_number()) %>% 
#   mutate(week = case_when(
#     week_num == 1 ~ as.Date("2019-05-20"),
#     week_num == 2 ~ as.Date("2019-05-27"),
#     week_num == 3 ~ as.Date("2019-06-03"),
#     week_num == 4 ~ as.Date("2019-06-17"),
#     week_num == 5 ~ as.Date("2019-07-01"),
#     week_num == 6 ~ as.Date("2019-07-15"),
#     week_num == 7 ~ as.Date("2019-07-29")
#   )) %>% 
#   select(-week_num)


# group_by(tag) %>%
# mutate(increment_cm = DBH_cm-lag(DBH_cm, default = first(DBH_cm))) %>% 
# #filter(!is.na(increment_cm)) %>% 
# mutate(RGR_cm = increment_cm / lag(DBH_cm, default = first(DBH_cm)))

# Convert diameter increments to biomass, 
# get increment_kg
# df_sub <- subcanopy_select %>% 
#   mutate(biomass_kg = biomass_a(species, DBH_cm)) %>% 
# group_by(tag) %>% 
#mutate(inc_days = abs(DOY-lag(DOY, default = first(DOY)))) %>% 
# #mutate(increment_kg = case_when(
#   biomass_kg < lag(biomass_kg, default = first(biomass_kg)) ~ 0,
#   biomass_kg >= lag(biomass_kg, default = first(biomass_kg)) ~ biomass_kg-lag(biomass_kg, default = first(biomass_kg))
# )) #%>% 
#mutate(RGR_biomass = increment_kg/biomass_kg)


create df of mean species delta biomass
# mean_sp <- sub_annual_inc %>% 
#   group_by(subplot, species) %>% 
#   summarize(mean_sp = mean(biomass_inc)) 
# 
# # create df of mean subplot delta biomass
# mean_subplot <- sub_annual_inc %>% 
#   group_by(subplot) %>% 
#   summarize(mean_subplot = mean(biomass_inc)) 

# join these two df's back to df_sub so there is now 2 new columns and every record 
# has either mean delta sp or mean_deta_subplot
# annual_sub1 <- right_join(mean_sp, sub_annual_inc)
# annual_sub_mean <- right_join(mean_subplot, annual_sub1)

##################################
#######ALright this is garbage; NEED TO: 
######## 1) pretty sure I can remove the subcanopy select chunk Since i generate RGR
########    from biomass increment in the above chunk (why do it twice?)
######## 2) summarize mean subplot and species biomass_inc (per tree) & mean subplot
########    biomass_inc, then apply them to the unsampled population; DON"T SHORTCUT!
######## 3) bring in inc_days cleaner; why are some non-whole numbers? 

##############getting close, i think I have issues with week 1 being either the week
##############of 05/20 or 05/29 - thus I don't have a record for every tree on both
##############of those dates 

# create mean delta biomass for each measured species and an overal sublpot mean
# delta biomass for unmeasured species (my sample did not capture every subcanopy)
# species in the subplot) 

# create df of mean species delta biomass
# delta_B_sp <- df_sub %>% 
#   group_by(subplot, week, species) %>% 
#   summarize(mean_delta_sp = mean(increment_kg)) 
# 
# # create df of mean subplot delta biomass
# delta_B_subplot <- df_sub %>% 
#   group_by(subplot, week) %>% 
#   summarize(mean_delta_subplot = mean(increment_kg)) 
# 
# # join these two df's back to df_sub so there is now 2 new columns and every record 
# # has either mean delta sp or mean_deta_subplot
# df_sub1 <- right_join(delta_B_sp, df_sub)
# df_sub_mean <- right_join(delta_B_subplot, df_sub1)
# 
# # now bring in the stems_subplot so we can scale to the subplot level from our sample
# # population. Then group by subplot, week and fill in mean_delta_subplot (biomass 
# # increment) so that every subplot and week has either a measured or mean delta B.
# # Do not need to fill mean_delta_sp because if there are species that do not have 
# # a increment_kg (measured increment), then there is no mean value for that species 
# # anyway. Then create delta_B_obs_est (measured and mean biomass increment) so that
# # every record has a biomass increment from which to scale. Note: where there is no 
# # measured or mean biomass increment is the initial week of measurement. Thus the 
# # increment for that week is coerced to 0. 
# 
# df_summary <- df_sub_mean %>% 
#   right_join(stems_subplot) %>% arrange(subplot, week, species) %>% 
#   group_by(subplot, week) %>% fill(mean_delta_subplot, .direction = "updown") %>% 
#   mutate(delta_B_obs_est = case_when(
#     !is.na(increment_kg) ~ increment_kg,
#     is.na(increment_kg) & !is.na(mean_delta_sp) ~ mean_delta_sp,
#     is.na(increment_kg) & is.na(mean_delta_sp) ~ mean_delta_subplot,
#     is.na(increment_kg) & is.na(mean_delta_sp) & is.na(mean_delta_subplot) ~ 0
#   )) %>% 
#   group_by(subplot, week, species, total_stems) %>% 
#   summarize(mean_biomass_inc = mean(delta_B_obs_est)) %>% 
#   mutate(stemsxbiomass = mean_biomass_inc * total_stems) %>% 
#   group_by(subplot, week) %>% 
#   summarise(subplot_biomass_ha = sum(stemsxbiomass) * 40) %>% 
#   # *40 b/c coutns were 0.025ha plots
#   mutate(kgC_ha = subplot_biomass_ha * 0.48) %>% 
#   group_by(subplot) %>% 
#   summarize(kgC_ha_yr = sum(kgC_ha))
#            
#            mean_delta_B * total_stems)# %>% 
#   
# mutate(B_per_ha = subplot_delta_B * 40) %>% # *40 b/c coutns were 0.025ha plots
#   mutate(kgC_per_ha = B_per_ha * 0.48) %>% 
#   mutate(kgC_ha_day = kgC_per_ha / inc_days_mean) %>% 
#   mutate(replicate = substr(subplot, 1, 1)) %>% 
#   mutate(severity = case_when(
#     subplot == "A01E" ~ 0.85, subplot == "A01W" ~ 0.85, subplot == "A02E" ~ 0.45,
#     subplot == "A02W" ~ 0.45, subplot == "A03E" ~ 0.65, subplot == "A03W" ~ 0.65,
#     subplot == "A04E" ~ 0.00, subplot == "A04W" ~ 0.00, subplot == "B01E" ~ 0.00,
#     subplot == "B01W" ~ 0.00, subplot == "B02E" ~ 0.45, subplot == "B02W" ~ 0.45,
#     subplot == "B03E" ~ 0.85, subplot == "B03W" ~ 0.85, subplot == "B04E" ~ 0.65,
#     subplot == "B04W" ~ 0.65, subplot == "C01E" ~ 0.00, subplot == "C01W" ~ 0.00,
#     subplot == "C02E" ~ 0.65, subplot == "C02W" ~ 0.65, subplot == "C03E" ~ 0.85,
#     subplot == "C03W" ~ 0.85, subplot == "C04E" ~ 0.45, subplot == "C04W" ~ 0.45, 
#     subplot == "D01E" ~ 0.00, subplot == "D01W" ~ 0.00, subplot == "D02E" ~ 0.85,
#     subplot == "D02W" ~ 0.85, subplot == "D03E" ~ 0.45, subplot == "D03W" ~ 0.45,
#     subplot == "D04E" ~ 0.65, subplot == "D04W" ~ 0.65
#   )) %>% 
#   mutate(treatment = case_when(
#     subplot == "A01E" ~ "bottom", subplot == "A01W" ~ "top", subplot == "A02E" ~ "top",
#     subplot == "A02W" ~ "bottom", subplot == "A03E" ~ "bottom", subplot == "A03W" ~ "top",
#     subplot == "A04E" ~ "bottom", subplot == "A04W" ~ "top", subplot == "B01E" ~ "bottom",
#     subplot == "B01W" ~ "top", subplot == "B02E" ~ "top", subplot == "B02W" ~ "bottom",
#     subplot == "B03E" ~ "bottom", subplot == "B03W" ~ "top", subplot == "B04E" ~ "top",
#     subplot == "B04W" ~ "bottom", subplot == "C01E" ~ "top", subplot == "C01W" ~ "bottom",
#     subplot == "C02E" ~ "bottom", subplot == "C02W" ~ "top", subplot == "C03E" ~ "bottom",
#     subplot == "C03W" ~ "top", subplot == "C04E" ~ "top", subplot == "C04W" ~ "bottom", 
#     subplot == "D01E" ~ "bottom", subplot == "D01W" ~ "top", subplot == "D02E" ~ "bottom",
#     subplot == "D02W" ~ "top", subplot == "D03E" ~ "bottom", subplot == "D03W" ~ "top",
#     subplot == "D04E" ~ "top", subplot == "D04W" ~ "bottom"
#   )) 
#   
# 
# # create df to bind with NPP_total
# NPP_subcanopy <- df_sub %>% 
#   filter(!is.na(increment_kg)) %>% #removes week 1 measurement w/o an increment
#   group_by(subplot, week) %>% 
#   summarize(mean_delta_B = mean(increment_kg), inc_days_mean = mean(inc_days)) %>% 
#   right_join(tree_counts) %>% 
#   mutate(subplot_delta_B = mean_delta_B * total_stems * 4) %>% 
#   group_by(subplot) %>% 
#   summarise(annual_inc = sum(subplot_delta_B)) %>%
#   mutate(severity = case_when(
#     subplot == "A01E" ~ 0.85, subplot == "A01W" ~ 0.85, subplot == "A02E" ~ 0.45,
#     subplot == "A02W" ~ 0.45, subplot == "A03E" ~ 0.65, subplot == "A03W" ~ 0.65,
#     subplot == "A04E" ~ 0.00, subplot == "A04W" ~ 0.00, subplot == "B01E" ~ 0.00,
#     subplot == "B01W" ~ 0.00, subplot == "B02E" ~ 0.45, subplot == "B02W" ~ 0.45,
#     subplot == "B03E" ~ 0.85, subplot == "B03W" ~ 0.85, subplot == "B04E" ~ 0.65,
#     subplot == "B04W" ~ 0.65, subplot == "C01E" ~ 0.00, subplot == "C01W" ~ 0.00,
#     subplot == "C02E" ~ 0.65, subplot == "C02W" ~ 0.65, subplot == "C03E" ~ 0.85,
#     subplot == "C03W" ~ 0.85, subplot == "C04E" ~ 0.45, subplot == "C04W" ~ 0.45, 
#     subplot == "D01E" ~ 0.00, subplot == "D01W" ~ 0.00, subplot == "D02E" ~ 0.85,
#     subplot == "D02W" ~ 0.85, subplot == "D03E" ~ 0.45, subplot == "D03W" ~ 0.45,
#     subplot == "D04E" ~ 0.65, subplot == "D04W" ~ 0.65
#   )) %>% 
#   mutate(treatment = case_when(
#     subplot == "A01E" ~ "bottom", subplot == "A01W" ~ "top", subplot == "A02E" ~ "top",
#     subplot == "A02W" ~ "bottom", subplot == "A03E" ~ "bottom", subplot == "A03W" ~ "top",
#     subplot == "A04E" ~ "bottom", subplot == "A04W" ~ "top", subplot == "B01E" ~ "bottom",
#     subplot == "B01W" ~ "top", subplot == "B02E" ~ "top", subplot == "B02W" ~ "bottom",
#     subplot == "B03E" ~ "bottom", subplot == "B03W" ~ "top", subplot == "B04E" ~ "top",
#     subplot == "B04W" ~ "bottom", subplot == "C01E" ~ "top", subplot == "C01W" ~ "bottom",
#     subplot == "C02E" ~ "bottom", subplot == "C02W" ~ "top", subplot == "C03E" ~ "bottom",
#     subplot == "C03W" ~ "top", subplot == "C04E" ~ "top", subplot == "C04W" ~ "bottom", 
#     subplot == "D01E" ~ "bottom", subplot == "D01W" ~ "top", subplot == "D02E" ~ "top",
#     subplot == "D02W" ~ "bottom", subplot == "D03E" ~ "bottom", subplot == "D03W" ~ "top",
#     subplot == "D04E" ~ "bottom", subplot == "D04W" ~ "top"
#   )) %>% 
#   select(subplot, NPP_subcan = annual_inc) 
# 
# # create df with a PFT column for PFT analysis; 
# # NO POGR (PFT = early) IN THIS DATA SET! so no need to sort by species 
# PFT_NPP_subcanopy <- NPP_subcanopy %>% 
#   mutate(PFT = case_when(
#     !is.na(NPP_subcan) ~ "late",
#     is.na(NPP_subcan) ~ "early"
#   ))

# this summarizes subcanopy NPP by severity   
# subcan_2019_sev <- NPP_subcanopy %>% 
#   group_by(severity) %>% 
#   summarise(NPP = mean(annual_inc), SE = std.error(annual_inc))

#####################################################################################
# #plot this shiiiiit
# ggplot(df_summary, aes(x=week, y=kgC_ha_day)) +
#   stat_summary(fun.y = mean, geom="line")+
#   stat_summary(fun.y = mean, geom="point")
#   #facet_grid(replicate ~ .)
# 
# 
# ggplot(subcanopy_inc, aes(x=DOY, y=increment, color=subplot)) +
#   stat_summary(fun.y = mean, geom="line")+
#   stat_summary(fun.y = mean, geom="point") +
#   ylim(-2.5,2.5)

# # create a function to assign a to each species 
# a <- function(species){
#   if (species == "ACRU"){
#     a <- 0.03117 
#   } else if (species == "ACPE"){
#     a <-  0.2040 
#   } else if (species == "ACSA"){
#     a <- 0.1693 
#   } else if (species == "AMEL"){
#     a <- 0.1630 
#   } else if (species == "FAGR"){
#     a <- 0.1892
#   } else if (species == "PIRE"){
#     a <- 0.0526
#   } else if (species == "PIST"){
#     a <- 0.0408 
#   } else if (species == "POGR"){
#     a <- 0.1387 
#   } else if (species == "QURU"){
#     a <- 0.0398
#   }
#   return(a)
# }
# 
# # write a function for b 
# 
# b <- function(species){
#   if (species == "ACRU"){
#     b <- 2.7780
#   } else if (species == "ACPE"){
#     b <- 2.2524 
#   } else if (species == "ACSA"){
#     b <- 2.3436
#   } else if (species == "AMEL"){
#     b <- 2.4940
#   } else if (species == "FAGR"){
#     b <- 2.3097
#   } else if (species == "PIRE"){
#     b <- 2.5258
#   } else if (species == "PIST"){
#     b <- 2.5735
#   } else if (species == "POGR"){
#     b <- 2.3498
#   } else if (species == "QURU"){
#     b <- 2.7734
#   }
#   return(b)
# }

# run my biomass function through the df 

# df <- subcanopy_select %>% 
#   mutate(a = a(species)) %>% 
#   mutate(b = b(species)) %>% 
#   group_by(tag) %>% 
#   filter(row_number()==1) %>% 
#   mutate(biomass_a = a * DBH_cm ^ b) %>% 
#   right_join(subcanopy_select) %>% 
#   fill(a) %>% fill(b)