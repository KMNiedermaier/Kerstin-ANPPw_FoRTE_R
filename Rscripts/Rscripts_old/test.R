#set library to dyplr
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotrix)

# use source to run subcanopy scripts 
source("Rscripts/1-8cm_increment&graph.R")

#importing csv
dendro_data <- read.csv("data/FoRTE_dendro_data_outliers_4.csv")

all_reps_all_trees <- read.csv("data/allreps_dist_treatment_output_update.csv")

#selected only for columns that I will use 
all_reps <- select(dendro_data, week, subplot, tag, 
                   fate, species, DBH_cm, a, b, bands_in, DOY, date, notes) %>%
  mutate(bands_cm= bands_in*2.54) %>% 
  filter(!is.na(bands_cm))

#select my good good columns
all_trees <- all_reps_all_trees %>% 
  select(SubplotID, Species, Tag, fate, dbh)

#rename column names to match all_reps
names(all_trees) <- c("subplot", "species", "tag", "fate", "DBH_cm")

#creating 16 of each tag to build a week column and prep for join; then join and fill in DOY 
all_trees_weeks <- all_trees %>% 
  slice(rep(1:n(), each = 17)) %>% 
  group_by(tag) %>% 
  mutate(week = row_number()) %>% 
  left_join(all_reps) %>% 
  group_by(subplot, week) %>% 
  fill(DOY) %>%
  fill(DOY, .direction = c("up")) %>% 
  ungroup()

all_trees_weeks$uniqueID <- paste(all_trees_weeks$subplot, all_trees_weeks$week,
                                  all_trees_weeks$species, sep = "_")
#building df with all trees sampled and unsampled
#moving from raw increment to both delta biomass AND RGR 
RGR_plot <- all_trees_weeks %>%
  group_by(species) %>% 
  fill(a) %>% fill(a, .direction = "up") %>% 
  fill(b) %>% fill(b, .direction = "up") %>% ungroup() %>%
  filter(week == 1 | notes == "new band" | notes == "new band; 182.2=after adjustment") %>%
  mutate(biomass_a= a*DBH_cm^b) %>% #initial biomass
  mutate(bandread_i= bands_cm) %>% #november band read 
  right_join(all_trees_weeks, by = c("subplot", "species", "tag", "DBH_cm", "week", 
                                     "bands_in", "DOY", "notes", "fate", "bands_cm", "uniqueID"),
             suffix = c("", ".y")) %>% #joining filtered table w/new biomass_a col
  group_by(tag) %>%
  fill(bandread_i) %>% fill(a) %>% fill(b)  #fill in missing values
#mutate(band_diff= bands_cm-bandread_i) %>% # maybe add in swell adj here 

# create function for band_diff (how much it has grown since initial bandread)
banddiff <- function(bands_cm, bandread_i, species, fate){
  if (species == "ACRU" & fate == "kill"){
    band_diff <- (bands_cm - bandread_i) * 0.24
  } else if (species == "ACSA" & fate == "kill"){
    band_diff <- (bands_cm - bandread_i) * 0.17
  } else {
    band_diff <- bands_cm - bandread_i
  }
  return(band_diff)
}

# loop the function through RGR_plot to create new band_diff column adjusted for 
# swelling in the "kill" ACRU and ACSA
for (i in 1:nrow(RGR_plot)){
  RGR_plot$band_diff[i] <- banddiff(RGR_plot$bands_cm[i], RGR_plot$bandread_i[i], 
                                    RGR_plot$species[i], RGR_plot$fate[i])
}

# continue pipe to get RGR 
RGR_plot <- RGR_plot %>% 
  mutate(biomass_b_kg= a*(DBH_cm+band_diff)^b) %>%
  mutate(increment_kg= biomass_b_kg-lag(biomass_b_kg, default = first(biomass_b_kg))) %>%
  mutate(RGR= increment_kg/biomass_b_kg) %>% 
  mutate(inc_days= abs(DOY-lag(DOY, default = first(DOY))))

###commented out previous pipeline, not accounting for swell ######
# RGR_plot <- all_trees_weeks %>%
#   group_by(species) %>% 
#   fill(a) %>% fill(a, .direction = "up") %>% 
#   fill(b) %>% fill(b, .direction = "up") %>% ungroup() %>%
#   filter(week == 1 | notes == "new band" | notes == "new band; 182.2=after adjustment") %>%
#   mutate(biomass_a= a*DBH_cm^b) %>% #initial biomass
#   mutate(bandread_i= bands_cm) %>% #november band read 
#   right_join(all_trees_weeks, by = c("subplot", "species", "tag", "DBH_cm", "week", 
#                                               "bands_in", "DOY", "notes", "bands_cm", "uniqueID"),
#                       suffix = c("", ".y")) %>% #joining filtered table w/new biomass_a col
#   group_by(tag) %>%
#   fill(bandread_i) %>% fill(a) %>% fill(b) %>% #fill in missing values 
#   mutate(band_diff= bands_cm-bandread_i) %>% # maybe add in swell adj here 
#   mutate(biomass_b_kg= a*(DBH_cm+band_diff)^b) %>%
#   mutate(increment_kg= biomass_b_kg-lag(biomass_b_kg, default = first(biomass_b_kg))) %>%
#   mutate(RGR= increment_kg/biomass_b_kg) %>% 
#   mutate(inc_days= abs(DOY-lag(DOY, default = first(DOY))))
  
    
subplot_wk_sp_RGR <- RGR_plot %>%
  filter(!is.na(RGR)) %>% 
  group_by(subplot, week, species) %>% 
  #summarise(mean_delta_biomass= mean(increment_kg)) %>% 
  summarise(mean_RGR= mean(RGR))

RGR_plot_mean <- full_join(subplot_wk_sp_RGR, RGR_plot)

##### sampled trees to generate regression equations to estimate growth ############

#selected only for columns that I will use 
all_reps <- select(dendro_data, severity, fate, treatment, week, subplot, tag, replicate, 
                   species, DBH_cm, a, b, bands_in, DOY, notes) %>%
  mutate(bands_cm= bands_in*2.54) %>% 
  filter(!is.na(bands_cm))

#starting biomass and band read; 
all_reps_filtered <- all_reps %>%
  group_by(tag) %>%
  filter(week == 1 | notes == "new band" | notes == "new band; 182.2=after adjustment") %>%
  mutate(biomass_a= a*DBH_cm^b) %>% #initial biomass
  mutate(bandread_i= bands_cm) %>% #november band read 
  full_join(all_reps) %>% #joining filtered table w/new biomass col
  group_by(tag) %>%
  fill(biomass_a) %>%
  fill(bandread_i) 
  #mutate(band_diff= bands_cm-bandread_i) 
  
# loop the function through RGR_plot to create new band_diff column adjusted for 
# swelling in the "kill" ACRU and ACSA
  for (i in 1:nrow(all_reps_filtered)){
    all_reps_filtered$band_diff[i] <- banddiff(all_reps_filtered$bands_cm[i], 
                                      all_reps_filtered$bandread_i[i], 
                                      all_reps_filtered$species[i], RGR_plot$fate[i])
  }  

all_reps_filtered1 <- all_reps_filtered %>% #calc inc percents 
  group_by(tag) %>%
  mutate(biomass_b_kg= a*(DBH_cm+band_diff)^b) %>%
  mutate(increment_kg= biomass_b_kg-lag(biomass_b_kg, default = first(biomass_b_kg))) %>%
  mutate(RGR= increment_kg/biomass_b_kg)

#select for RGR, tag, date, subplot, species, DBH
RGR_plot_regression <- all_reps_filtered1 %>%
  select(subplot, week, tag, species, 
         DBH_cm, DOY, biomass_a, RGR, increment_kg) %>%
  group_by(tag) %>%
  mutate(inc_days= abs(DOY-lag(DOY, default = first(DOY))))  
#week two is not the right number of days!!

RGR_plot$uniqueID <- paste(RGR_plot$subplot, RGR_plot$week,
                           RGR_plot$species, sep = "_")

#scale to NPP by subplot, week, species, AND DBH (where necessary)
#runs regression for each subplot, week, and species AND returns object with 
#coefficients
regressions <- RGR_plot_regression %>% 
  group_by(subplot, week, species) %>% 
  do(mod = summary(lm(RGR ~ DBH_cm, data = .))) #%>% #-1 removes the intercept

regressions$uniqueID <- paste(regressions$subplot, regressions$week, 
                              regressions$species, sep = "_")

final_output <- matrix(nrow= 2466, ncol=4, dimnames = list(c(), c("uniqueID", "pvalue", "slope", "intercept")))
regressions$uniqueID <- paste(regressions$subplot, regressions$week, 
                              regressions$species, sep = "_")

for(i in 1:nrow(regressions)) {
  
  final_output[i, 4] <- regressions$mod[[i]][4]$coefficients[1]
  final_output[i, 3] <- regressions$mod[[i]][4]$coefficients[2]
  final_output[i, 2] <- regressions$mod[[i]][4]$coefficients[8]
  final_output[i, 1] <- regressions$uniqueID[[i]] #may need [,i] or [i,]
  
} 

#convert matrix to data frame 
final_output_df <- as.data.frame(final_output) 

#remove NA's: NEED TO FIGURE OUT WHY NA'S VS NaN'S 
final_output_df_NAomit <-   filter(final_output_df, pvalue != "NaN", !is.na(pvalue)) 

#coerce pvalue, slope, and intercept into character and then numeric
final_output_df_NAomit$pvalue <-  as.numeric(as.character(final_output_df_NAomit$pvalue))
final_output_df_NAomit$slope <-  as.numeric(as.character(final_output_df_NAomit$slope))
final_output_df_NAomit$intercept <-  as.numeric(as.character(final_output_df_NAomit$intercept))

#filter for pvalue's < 0.05; subplot, weeks, and species where DBH is having a sig effect 
sig_DBH_effect <- filter(final_output_df_NAomit, pvalue < 0.05)

#create new RGR_obs_est column for calculated RGR (sampled trees), avg_RGR (unsampled trees), 
#adj_RGR (for unique id's where DBH had significant effect on RGR)
RGR_obs_est_df <- RGR_plot_mean %>% 
  left_join(sig_DBH_effect) %>% 
  mutate(RGR_obs_est = case_when(
    !is.na(RGR) & is.na(slope) ~ RGR,
    !is.na(RGR) & !is.na(slope) ~ DBH_cm*slope+intercept,
    is.na(RGR) & is.na(slope) ~ mean_RGR,
    is.na(RGR) & !is.na(slope) ~ DBH_cm*slope+intercept)) %>% 
  mutate(biomass_b_kg_new = biomass_a + biomass_a*RGR_obs_est) %>% 
  select(uniqueID, tag, species, DBH_cm, biomass_a, biomass_b_kg_new, fate, 
         RGR_obs_est, DOY, inc_days) %>% 
  arrange(tag)

df <- data.frame(RGR_obs_est_df)
# this should set the RGR to 0 if negative
df$RGR_obs_est[df$RGR_obs_est < 0] <- 0
#split into list of data frames
df.list <- split(df, df$tag)
df.list <- lapply(df.list, function(x){
  for (i in 1:nrow(x)) {
    for(j in 2:nrow(x)){
      x$biomass_a[j] = x$biomass_b_kg_new[j - 1]
      x$biomass_b_kg_new[j] = (x$biomass_a[j] * x$RGR_obs_est[j]) + x$biomass_a[j]
    }
    return(x)
  }
})
df <- plyr::ldply(df.list, data.frame)

# let's get some NPP bebe
inc_days <- RGR_plot %>% 
  group_by(subplot, week) %>% # MUST ADD SPECIES HERE TOO!!! 
  summarise(inc_days= mean(inc_days))

NPP_final <- df %>% 
  filter(!is.na(biomass_a)) %>% #NA's are from species that are not present in a subplot; need to remove them 
  mutate(biomass_diff = biomass_b_kg_new - biomass_a) %>% 
  group_by(subplot, week) %>% # can ADD SPECIES in here if i want to look at effect of species 
  summarise(subplot_biomass_kg = sum(biomass_diff), DOY = mean(DOY)) %>% 
  mutate(biomass_per_ha = subplot_biomass_kg*10) %>% 
  mutate(kgC_per_ha = biomass_per_ha*0.48) %>% 
  right_join(inc_days) %>% 
  mutate(kgC_per_ha_day = kgC_per_ha / inc_days) %>% 
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
  )) 


# getting DOY into dates and correcting the november 2018 date 
dates <- data.frame(as.factor(as.Date(NPP_final$DOY, "2019-01-01")))
names(dates) <- c("date")
dates <- data.frame(recode_factor(dates$date, '2019-11-03' = "2018-11-15"))
names(dates) <- c("date")
dates$date <- as.Date(dates$date,"%Y-%m-%d")

NPP_final <- NPP_final %>% 
  bind_cols(dates) 

#create week and month columns to look at different temporal scales 
NPP_final$month <- as.Date(cut(NPP_final$date, breaks = "month"))
NPP_final$week <- as.Date(cut(NPP_final$date, breaks = "week", start.on.monday = FALSE))

# sum all the species in each week, subplot
NPP_total <- NPP_final %>% 
  filter(!is.na(kgC_per_ha_day)) %>% 
  group_by(week, subplot, severity, treatment) %>% 
  summarise(NPP_canopy = sum(kgC_per_ha_day)) %>% 
  left_join(NPP_subcanopy) 
  
# make NA's into 0's so we can sum canopy and subcanopy 
NPP_total[is.na(NPP_total)] <- 0

# sum canopy andsubcanopy for total NPP 
NPP_total <- mutate(NPP_total, NPP_total = NPP_canopy + NPP_subcan)

# sum all weeks for cumulative NPP 
NPP_cumulative <- NPP_total %>% 
  group_by(subplot, severity, treatment) %>% 
  summarise(NPP_season = sum(NPP_total))

#### Plot df's with summary data and SE 

# severity time series df 
severity_time <- NPP_total %>% 
  group_by(week, severity) %>% 
  summarize(NPP = mean(NPP_canopy), SE = std.error(NPP_canopy))

# Cumulative NPP df for figures
histogram <- NPP_total %>% 
  group_by(subplot, severity) %>%
  summarise(NPP_subplot_can = sum(NPP_canopy),
              NPP_subplot_sub = sum(NPP_subcan)) %>% 
  group_by(severity) %>%
  summarise(canopy_sum = sum(NPP_subplot_can),
            subcan_sum = sum(NPP_subplot_sub)) %>% 
  group_by(severity) %>% 
  gather(NPP_type, kgC_ha_yr, canopy_sum, subcan_sum) %>% 
  mutate(SE = case_when(
    NPP_type == "canopy_sum" ~ kgC_ha_yr * 0.07,
    NPP_type == "subcan_sum" ~ kgC_ha_yr * 0.07
  ))
  
# make severity %'s for x axis
histogram$severity <- recode(histogram$severity, "0.00" = "0%", 
                                    "0.45" = "45%", "0.65" = "65%",
                                    "0.85" = "85%")

# treatement df for figures 
treatment <- NPP_total %>% 
  filter(severity == 0.45) %>% 
  group_by(week, treatment) %>% 
  summarise(NPP = mean(NPP_canopy), SE = std.error(NPP_canopy))


##################################################################################
# PLOTS!!

# forte colors decided by Jeff 
forte <- c("#000000", "#009E73", "#0072B2", "#D55E00")

## Disturbance severity time series; includes CANOPY only 
ggplot(severity_time, aes(x=week, y=NPP, color = factor(severity))) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = NPP-SE, ymax = NPP+SE), width = 0.8) +
  scale_colour_manual(values = forte, name = "Disturbance\nSeverity", labels = 
                        c("0%", "45%", "65%", "85%")) +
  labs(x = "" , 
       y = expression(paste("NPP (",kgC," ",ha^-1," ",day^-1,")"))) +
  theme_classic() +
  theme(legend.position = c(0.8, 0.6)) +
  geom_vline(xintercept = as.numeric(as.Date("2019-05-21")), linetype = 2, 
             color = "firebrick3")+
  annotate(geom="text", x = as.Date("2019-06-02", "%Y-%m-%d"), y = 60, label = c("*"),
           color="black", size = 5) +
  annotate(geom="text", x = as.Date("2019-06-09", "%Y-%m-%d"), y = 100, label = c("*"),
         color="black", size = 5) +
  annotate(geom="text", x = as.Date("2019-07-21", "%Y-%m-%d"), y = 120, label = c("*"),
           color="black", size = 5) +
  annotate(geom="text", x = as.Date("2019-07-28", "%Y-%m-%d"), y = 110, label = c("*"),
           color="black", size = 5) +
  annotate(geom="text", x = as.Date("2019-08-04", "%Y-%m-%d"), y = 70, label = c("*"),
           color="black", size = 5) +
  annotate(geom="text", x = as.Date("2019-08-11", "%Y-%m-%d"), y = 60, label = c("*"),
           color="black", size = 5) +
  annotate(geom="text", x = as.Date("2019-04-20", "%Y-%m-%d"), y = 130, 
           label = c("p = 0.053"),
           color="black", size = 3.5) 
             
ggsave("figures/NPP_Severity_timeseries.png", plot = last_plot())

## Cumulative NPP histogram 
ggplot(histogram, aes(x = factor(severity), y = kgC_ha_yr, fill = NPP_type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5) +
  theme_classic() +
  scale_fill_manual(values = c("lightgreen", "goldenrod1"), name = "", 
                    labels = c("Canopy", "Subcanopy")) +
  labs(x = "Disturbance Severity",
       y = expression(paste("Cumulative NPP ( ",kgC," ",ha^-1," ",yr^-1,")"))) +
  theme(axis.title = element_text(size = 10)) +
  geom_errorbar(aes(ymin = kgC_ha_yr-SE, ymax = kgC_ha_yr+SE), width = 0.1) +
  annotate(geom="text", x = "0%", y = 5000, label = c("a"),
           color="black", size = 3) +
  annotate(geom="text", x = "45%", y = 5000, label = c("a"),
           color="black", size = 3) +
  annotate(geom="text", x = "65%", y = 5000, label = c("a"),
           color="black", size = 3) +
  annotate(geom="text", x = "85%", y = 5000, label = c("a"),
           color="black", size = 3) 

ggsave("figures/cumulative_NPP.jpeg", plot = last_plot())
  
forte <- c("#000000", "#009E73", "#0072B2", "#D55E00")

################# Treatment by severity df's and figures 

# treatement df for figures 
treatment <- NPP_total %>% 
  filter(severity == 0.45) %>% 
  group_by(week, treatment) %>% 
  summarise(NPP = mean(NPP_canopy), SE = std.error(NPP_canopy))

## Treatment-severity time series 
plot_tx45 <- ggplot(treatment, aes(x = week, y = NPP, linetype = treatment)) +
  geom_point(color = "#009E73") +
  geom_line(color = "#009E73") +
  geom_errorbar(aes(ymin = NPP-SE, ymax = NPP+SE), width = 0.8, color = "#009E73") +
  scale_linetype_manual(values = c(1, 6), name = "Treatment", 
                        labels = c("Bottom Up", "Top Down")) +
  labs(x = "" , 
       y = expression(paste("NPP (",kgC," ",ha^-1," ",day^-1,")"))) +
  ylim(0, 170) +
  theme_classic() +
  theme(legend.position = c(0.8, 0.6)) +
  geom_vline(xintercept = as.numeric(as.Date("2019-05-21")), linetype = 2, 
             color = "firebrick3") +
  ggtitle("45% Defoliation") +
  theme(plot.title = element_text(hjust=0.5)) +
  annotate(geom="text", x = as.Date("2019-04-20", "%Y-%m-%d"), y = 130, 
           label = c("N.S."),
           color="black", size = 3.5) 

# treatement df for figures 
treatment <- NPP_total %>% 
  filter(severity == 0.65) %>% 
  group_by(week, treatment) %>% 
  summarise(NPP = mean(NPP_canopy), SE = std.error(NPP_canopy))

## Treatment-severity time series 
plot_tx65 <- ggplot(treatment, aes(x = week, y = NPP, linetype = treatment)) +
  geom_point(color = "#0072B2") +
  geom_line(color = "#0072B2") +
  geom_errorbar(aes(ymin = NPP-SE, ymax = NPP+SE), width = 0.8, color = "#0072B2") +
  scale_linetype_manual(values = c(1, 6), name = "Treatment", 
                        labels = c("Bottom Up", "Top Down")) +
  labs(x = "" , 
       y = expression(paste("NPP (",kgC," ",ha^-1," ",day^-1,")"))) +
  ylim(0,170) +
  theme_classic() +
  theme(legend.position = "none") +
  geom_vline(xintercept = as.numeric(as.Date("2019-05-21")), linetype = 2, 
             color = "firebrick3") +
  ggtitle("65% Defoliation") +
  theme(plot.title = element_text(hjust=0.5)) +
  annotate(geom="text", x = as.Date("2019-04-20", "%Y-%m-%d"), y = 130, 
           label = c("p < 0.05"),
           color="black", size = 3.5) 

# treatement df for figures 
treatment <- NPP_total %>% 
  filter(severity == 0.85) %>% 
  group_by(week, treatment) %>% 
  summarise(NPP = mean(NPP_canopy), SE = std.error(NPP_canopy))

## Treatment-severity time series 
plot_tx85 <- ggplot(treatment, aes(x = week, y = NPP, linetype = treatment)) +
  geom_point(color = "#D55E00") +
  geom_line(color = "#D55E00") +
  geom_errorbar(aes(ymin = NPP-SE, ymax = NPP+SE), width = 0.8, color = "#D55E00") +
  scale_linetype_manual(values = c(1, 6), name = "Treatment", 
                        labels = c("Bottom Up", "Top Down")) +
  labs(x = "" , 
       y = expression(paste("NPP (",kgC," ",ha^-1," ",day^-1,")"))) +
  ylim(0,170) +
  theme_classic() +
  theme(legend.position = "none") +
  geom_vline(xintercept = as.numeric(as.Date("2019-05-21")), linetype = 2, 
             color = "firebrick3") +
  ggtitle("85% Defoliation") +
  theme(plot.title = element_text(hjust=0.5)) +
  annotate(geom="text", x = as.Date("2019-04-20", "%Y-%m-%d"), y = 130, 
           label = c("p < 0.05"),
           color="black", size = 3.5) 

# paste all 3 tx figures onto 1 
library(gridExtra)
timeseries <- grid.arrange(plot_tx45, plot_tx65, plot_tx85, nrow=1)
g <- arrangeGrob(plot_tx45, plot_tx65, plot_tx85, nrow=1)
ggsave("figures/treatment.jpeg", height = 2.5, width = 10, units = "in", g)
###################################################################





# ggplot(NPP_final, aes(x=week, y=kgC_per_ha_day, color=subplot)) +
#   geom_point()+
#   ggtitle("NPP All Plots")
# 
# 
# ggplot(NPP_total, aes(x=week, y=NPP_canopy, color = factor(severity))) +
#   stat_summary(fun.y = mean, geom="line")+
#   stat_summary(fun.y = mean, geom="point")+
#   #facet_grid(replicate ~ .) +
#   ggtitle("NPP Severity Gradient")
# 
# ggplot(NPP_final, aes(x=week, y=NPP_subplot, color=treatment)) +
#   stat_summary(fun.y = mean, geom="line")+
#   stat_summary(fun.y = mean, geom="point")+
#   #facet_grid(replicate ~ .) +
#   ggtitle("NPP Treatments")
# 
# ggplot(NPP_final, aes(x=week, y=NPP_subplot, color=fate)) +
#   stat_summary(fun.y = mean, geom="line")+
#   stat_summary(fun.y = mean, geom="point")+
#   #facet_grid(replicate ~ .) +
#   ggtitle("NPP Treatments")
# 
# NPP_final$tx_sev <- paste(NPP_final$treatment, NPP_final$severity, sep = "_")
# 
# ggplot(NPP_final, aes(x=week, y=NPP_subplot, color=tx_sev)) +
#   stat_summary(fun.y = mean, geom="line")+
#   stat_summary(fun.y = mean, geom="point")+  
#   #facet_grid(treatment ~ .) +
#   ggtitle("NPP Treatments/severity")
