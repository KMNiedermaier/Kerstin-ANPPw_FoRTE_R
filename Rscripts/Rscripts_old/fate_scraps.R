#########################################################################
######### NPP Fate by week 
# NPP_fate_wk <- df %>% 
#   group_by(tag) %>% # NEW WIH inc_cm
#   #filter(!is.na(biomass_a)) %>% #NA's are from species that are not present in a subplot; need to remove them 
#   mutate(biomass_new = case_when(
#     species == "AMEL" ~ (a*(DBH_cm_new*10)^b)/1000,
#     species != "AMEL" ~ a*DBH_cm_new^b
#   )) %>%  # NEW WITH inc_cm 
#   group_by(subplot, week, fate) %>% # can ADD SPECIES in here if i want to look at effect of species 
#   summarise(subplot_biomass_kg = sum(biomass_new), DOY = mean(DOY)) %>% 
#   group_by(subplot, fate) %>% 
#   arrange(subplot, fate) %>% 
#   mutate(biomass_diff = subplot_biomass_kg - lag(subplot_biomass_kg, 
#                                                  default = first(subplot_biomass_kg))) %>% 
#   mutate(biomass_per_ha = biomass_diff*10) %>% 
#   mutate(kgC_per_ha = biomass_per_ha*0.48) %>% 
#   right_join(inc_days) %>% 
#   mutate(kgC_per_ha_day = kgC_per_ha / inc_days) %>% # this is kgC/ha/day; call NPP_can for the gather function
#   mutate(replicate = substr(subplot, 1, 1)) %>% 
#   filter(!is.na(kgC_per_ha_day)) %>% 
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
# # getting DOY into dates and correcting the november 2018 date 
# dates <- data.frame(as.factor(as.Date(NPP_fate_wk$DOY, "2019-01-01")))
# names(dates) <- c("date")
# dates <- data.frame(recode_factor(dates$date, '2019-11-03' = "2018-11-15"))
# names(dates) <- c("date")
# dates$date <- as.Date(dates$date,"%Y-%m-%d")
# 
# NPP_fate_wk <- NPP_fate_wk %>% 
#   bind_cols(dates) 
# 
# #create week and month columns to look at different temporal scales 
# NPP_fate_wk$month <- as.Date(cut(NPP_fate_wk$date, breaks = "month"))
# NPP_fate_wk$week <- as.Date(cut(NPP_fate_wk$date, breaks = "week", start.on.monday = FALSE))
#   
# #####################################################################################
# ######## NPP_fate by tag (tree) to keep all annual data in order to make box plots 
# 
# NPP_fate_boxplot <- df %>%
#   group_by(tag) %>% # NEW WIH inc_cm
#   mutate(biomass_new = case_when(
#     species == "AMEL" ~ (a*(DBH_cm_new*10)^b)/1000,
#     species != "AMEL" ~ a*DBH_cm_new^b
#   )) %>%  # NEW WITH inc_cm
#   mutate(biomass_diff = biomass_new - lag(biomass_new, default = first(biomass_new))) %>%
#   group_by(subplot, fate, tag) %>% # can ADD SPECIES in here if i want to look at effect of species
#   summarise(annual_biomass_kg = sum(biomass_diff)) %>%
#   mutate(biomass_per_ha = annual_biomass_kg*10) %>%
#   mutate(kgC_ha_yr = biomass_per_ha*0.48) %>%
#   #mutate(replicate = substr(subplot, 1, 1)) %>%
#   #filter(!is.na(kgC_per_ha_day)) %>%
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
#   )) 
# # # make NA's into 0's so we can sum canopy and subcanopy 
# # NPP_fate_total[is.na(NPP_fate_total)] <- 0
# 
# #################################################################################
# #### Plot df's with summary data and SE 
# 
# # severity time series df; can make it with or without seerity, but if it's without 
# # severity, control should be filtered out 
# # can filter for different severities  or fates to produce different plots
# fate_wk <- NPP_fate_wk %>% 
#   #filter(severity == 0.85) %>% 
#   group_by(week, fate) %>% 
#   summarize(NPP = mean(kgC_per_ha_day), SE = std.error(kgC_per_ha_day)) %>% 
#   mutate(DOY = yday(week))
# 
# #severity_fate_time$fate_sev <- paste(severity_fate_time$fate, 
#                                #severity_fate_time$severity, sep = "_")
# 
# ##################################################################################
# # PLOTS!!
# 
# # forte colors decided by Jeff 
# forte <- c("#000000", "#009E73", "#0072B2", "#D55E00")
# kill_colors <- c("#009E73", "#0072B2", "#D55E00")
# live_kill_colors <- c("brown4", "forestgreen")
# 
# ## LIVE-KILL time series; includes CANOPY only; excludes control plots 
# plot45 <- ggplot(severity_fate_time, aes(x=week, y=NPP, color = fate)) +
#   geom_point() +
#   geom_line() +
#   geom_errorbar(aes(ymin = NPP-SE, ymax = NPP+SE), width = 0.8) +
#   scale_colour_manual(values = live_kill_colors, name = "Fate", 
#                       labels = c("girdled", "live")) +
#   labs(x = "" , 
#        y = expression(paste("NPP (",kgC," ",ha^-1," ",day^-1,")"))) +
#   ylim(0,110) +
#   theme_classic() +
#   theme(legend.position = c(0.8, 0.6)) +
#   geom_vline(xintercept = as.numeric(as.Date("2019-05-21")), linetype = 2, 
#              color = "firebrick3") +
#   ggtitle("45% Defoliation") +
#   theme(plot.title = element_text(hjust=0.5)) +
#   annotate(geom="text", x = as.Date("2019-04-20", "%Y-%m-%d"), y = 92, 
#            label = c("p < 0.05"),
#            color="black", size = 3.5) 
# 
# # filter for 65%
# severity_fate_time <- NPP_fate_total %>% 
#   filter(severity == 0.65) %>% 
#   group_by(week, fate) %>% 
#   summarize(NPP = mean(NPP_canopy), SE = std.error(NPP_canopy)) 
# 
# # plot 65% 0
# 
# #################################################################################
# ######### fate figures along severity gradient with RGR (NOT NPP!) #########
# 
# 
# RGR_fate <- df %>% 
#   mutate(replicate = substr(subplot, 1, 1)) %>%
#   group_by(replicate, subplot, week, fate) %>% # can ADD SPECIES in here if i want to look at effect of species 
#   summarise(RGR = mean(RGR_obs_est), DOY = mean(DOY)) %>% 
#   right_join(inc_days) %>% 
#   mutate(RGR = RGR/inc_days) %>% # get the per day RGR's 
#   filter(!is.na(RGR)) %>% 
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
#   )) 
# 
# # getting DOY into dates and correcting the november 2018 date 
# dates <- data.frame(as.factor(as.Date(RGR_fate$DOY, "2019-01-01")))
# names(dates) <- c("date")
# dates <- data.frame(recode_factor(dates$date, '2019-11-03' = "2018-11-15"))
# names(dates) <- c("date")
# dates$date <- as.Date(dates$date,"%Y-%m-%d")
# 
# RGR_fate <- RGR_fate %>% 
#   bind_cols(dates) 
# 
# #create week and month columns to look at different temporal scales 
# RGR_fate$month <- as.Date(cut(RGR_fate$date, breaks = "month"))
# RGR_fate$week <- as.Date(cut(RGR_fate$date, breaks = "week", start.on.monday = FALSE))
# 
# # create df for live vs kill histogram
# live_kill_hist <- RGR_fate %>% 
#   group_by(fate) %>% 
#   summarise(mean_RGR = mean(RGR)*10, SE = std.error(RGR)) # units are mm/cm
# 
# # histogram of live_kill RGR 
# ggplot(live_kill_hist, aes(x = fate, y = mean_RGR))+
#   geom_bar(aes(fill = fate), stat = "identity", width = 0.5)+
#   theme_classic()+
#   geom_errorbar(aes(ymin = mean_RGR-SE, ymax = mean_RGR+SE), width = 0.1) +
#   scale_fill_manual(values = c("firebrick2","darkgreen"), limits = c("kill", "live"),
#                     name = "Fate") +
#   theme(legend.position = "none")+
#   labs(x = "",
#        y = expression(paste("Mean RGR ( ",mm," ",cm^-1," ",day^-1,")")))+
#   annotate(geom="text", x = 2, y = 0.0065, label = c("p = 0.021"),
#            color="black", size = 5)
# 
# ggsave(filename = "figures/live_kill_RGR.jpeg", plot = last_plot())
# 
# # histogram for live-kill mean NPP
# ggplot(NPP_fate_mean, aes(x = factor(severity), y = NPP_mean, fill = fate)) +
#   geom_bar(stat = "identity", position = "dodge", width = 0.4) +
#   theme_classic() +
#   scale_fill_manual(values = c("firebrick2","darkgreen"), name = "fate")+
#   labs(x = "Disturbance Severity",
#        y = expression(paste("NPP (",kgC," ",ha^-1," ",day^-1,")"))) +
#   theme(axis.title = element_text(size = 10))+
#   geom_errorbar(position = "dodge", aes(ymin = NPP_mean-SE, ymax = NPP_mean+SE), width = 0.4) 
#   
# 
# 
# # create df for figures
# 
# RGR_severity_fate <- RGR_fate %>% 
#   filter(severity == 0.45, week != "2018-11-11") %>% 
#   mutate(RGR = RGR*10) %>% # units are now in mm per cm
#   group_by(week, fate) %>% 
#   summarize(RGR_mean = mean(RGR), SE = std.error(RGR)) 
# 
# ## LIVE-KILL (RGR) time series; includes CANOPY only; excludes control plots 
# RGR_plot45 <- ggplot(RGR_severity_fate, aes(x=week, y=RGR_mean, color = fate)) +
#   geom_point() +
#   geom_line() +
#   geom_errorbar(aes(ymin = RGR_mean-SE, ymax = RGR_mean+SE), width = 0.8) +
#   scale_colour_manual(values = live_kill_colors, name = "Fate", 
#                       labels = c("girdled", "live")) +
#   labs(x = "" , 
#        y = expression(paste("RGR (",mm," ",cm^-1," ",day^-1,")"))) +
#   #ylim(0, 0.018) +
#   theme_classic() +
#   theme(legend.position = c(0.8, 0.6)) +
#   geom_vline(xintercept = as.numeric(as.Date("2019-05-21")), linetype = 2, 
#              color = "firebrick3") +
#   ggtitle("45% Defoliation") +
#   theme(plot.title = element_text(hjust=0.5)) 
# 
# # LIVE-KILL (RGR) 65%
# 
# RGR_severity_fate <- RGR_fate %>% 
#   filter(severity == 0.65, week != "2018-11-11") %>% 
#   mutate(RGR = RGR*10) %>%
#   group_by(week, fate) %>% 
#   summarize(RGR_mean = mean(RGR), SE = std.error(RGR)) 
# 
# RGR_plot65 <- ggplot(RGR_severity_fate, aes(x=week, y=RGR_mean, color = fate)) +
#   geom_point() +
#   geom_line() +
#   geom_errorbar(aes(ymin = RGR_mean-SE, ymax = RGR_mean+SE), width = 0.8) +
#   scale_colour_manual(values = live_kill_colors, name = "Fate", 
#                       labels = c("girdled", "live")) +
#   labs(x = "" , 
#        y = expression(paste("RGR (",mm," ",cm^-1," ",day^-1,")"))) +
#   #ylim(0, 0.015) +
#   theme_classic() +
#   theme(legend.position = c(0.8, 0.6)) +
#   geom_vline(xintercept = as.numeric(as.Date("2019-05-21")), linetype = 2, 
#              color = "firebrick3") +
#   ggtitle("65% Defoliation") +
#   theme(plot.title = element_text(hjust=0.5)) 
# 
# # LIVE-KILL (RGR)  85%
# 
# RGR_severity_fate <- RGR_fate %>% 
#   filter(severity == 0.85, week != "2018-11-11") %>% 
#   mutate(RGR = RGR*10) %>%
#   group_by(week, fate) %>% 
#   summarize(RGR_mean = mean(RGR), SE = std.error(RGR)) 
# 
# RGR_plot85 <- ggplot(RGR_severity_fate, aes(x=week, y=RGR_mean, color = fate)) +
#   geom_point() +
#   geom_line() +
#   geom_errorbar(aes(ymin = RGR_mean-SE, ymax = RGR_mean+SE), width = 0.8) +
#   scale_colour_manual(values = live_kill_colors, name = "Fate", 
#                       labels = c("girdled", "live")) +
#   labs(x = "" , 
#        y = expression(paste("RGR (",mm," ",cm^-1," ",day^-1,")"))) +
#  # ylim(0, 0.018) +
#   theme_classic() +
#   theme(legend.position = c(0.8, 0.6)) +
#   geom_vline(xintercept = as.numeric(as.Date("2019-05-21")), linetype = 2, 
#              color = "firebrick3") +
#   ggtitle("85% Defoliation") +
#   theme(plot.title = element_text(hjust=0.5)) 
# 
# # save multipanel figure with all three severities 
# library(gridExtra)
# timeseries <- grid.arrange(RGR_plot45, RGR_plot65, RGR_plot85, nrow=1)
# g <- arrangeGrob(RGR_plot45, RGR_plot65, RGR_plot85, nrow=1)
# ggsave("figures/RGR_fates.jpeg",height = 2.5, width = 10, units = "in", g)