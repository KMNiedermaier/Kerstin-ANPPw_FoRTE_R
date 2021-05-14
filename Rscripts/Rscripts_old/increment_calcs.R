####THIS IS OLD AND OUTDATED 


#importing csv
dendro_data <- read.csv("data/FoRTE_dendro_data_outliers_4.csv")

#set library to dyplr
library(dplyr)
library(tidyr)
library(ggplot2)

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
  fill(bandread_i) %>%
  mutate(band_diff= bands_cm-bandread_i) 
  
#all_reps_filtered[all_reps_filtered < 0] <- 0   #take out (-)'s 

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
#RGR_plot[RGR_plot < 0] <- 0

#removing outliers strategy
RGR_plot[order(-RGR_plot$RGR), , drop = FALSE] #returns 10 most + rgr's
RGR_plot[order(RGR_plot$RGR), , drop = FALSE] #returns 10 most - rgr's
most_neg <- RGR_plot[order(RGR_plot$RGR), , drop = FALSE]
most_pos <- RGR_plot[order(-RGR_plot$RGR), , drop = FALSE]
################################################################
#PLOTS 

#plot RGR; April-Present; ALL READINGS
ggplot(RGR_plot, aes(x=DOY, y=increment_percent, color=subplot)) +
  geom_point()+
  xlim(75,250)+
  ggtitle("% RGR All Individuals")

#plot RGR: SUBPLOT MEAN per rep
ggplot(RGR_plot, aes(x=DOY, y=increment_percent, color=subplot)) +
  stat_summary(fun.y = mean, geom="line")+
  stat_summary(fun.y = mean, geom="point")+
  xlim(75,220)+
  facet_grid(replicate ~ .)+
  ggtitle("SUBPLOT mean RGR")

#plot RGR: SUBPLOT MEAN overall
ggplot(RGR_plot, aes(x=DOY, y=increment_percent, color=subplot)) +
  stat_summary(fun.y = mean, geom="line")+
  stat_summary(fun.y = mean, geom="point")+
  facet_grid(fate ~ .)+
  xlim(75,220)+
  ggtitle("SUBPLOT mean RGR")

#plot RGR: OVERALL MEAN
ggplot(RGR_plot, aes(x=week, y=increment_percent)) +
  stat_summary(fun.y = mean, geom="line")+
  stat_summary(fun.y = mean, geom="point")+
  ggtitle("OVERALL mean RGR")

#plot RGR: SEVERITY MEAN 
ggplot(RGR_plot, aes(x=week, y=increment_percent, color=factor(severity))) +
  stat_summary(fun.y = mean, geom="line")+
  stat_summary(fun.y = mean, geom="point")+
  ggtitle("Mean Relative Growth by Severity")+
  ylab("Mean Increment %")+
  scale_x_discrete(name= "Week",
                   limits= c("1","2","3","4","5",
                             "6","7","8","9","10",
                             "11","12","13"))

#plot RGR: FATE-SEVERITY MEAN
ggplot(RGR_plot, aes(x=DOY, y=increment_percent, color=factor(severity))) +
  stat_summary(fun.y = mean, geom="line")+
  stat_summary(fun.y = mean, geom="point")+
  facet_grid(fate ~ .)+
  xlim(85,240)+
  ggtitle("FATE-SEVERITY RGR")

ggplot(RGR_plot, aes(x=week, y=increment_percent, color=factor(severity))) +
  stat_summary(fun.y = mean, geom="line")+
  stat_summary(fun.y = mean, geom="point")+
  facet_grid(replicate ~ .)+
  ggtitle("SEVERITY per rep")

ggplot(RGR_plot, aes(x=DOY, y=increment_percent, color=fate)) +
  stat_summary(fun.y = mean, geom="line")+
  stat_summary(fun.y = mean, geom="point")+
  xlim(75,210)+
  ylab("Mean Increment %")+
  facet_grid(replicate ~ .)+
  ggtitle("Fate per Replicate")

#plot all RGR; Nov data on wrong side 
ggplot(RGR_plot, aes(x=DOY, y=increment_percent, color=subplot)) +
  geom_point()
  
#severity on the x-axis
ggplot(RGR_plot, aes(x=factor(severity), y=increment_percent, color=factor(severity))) +
  stat_summary(fun.y = mean, geom="line")+
  stat_summary(fun.y = mean, geom="point")+
  ggtitle("Mean Relative Growth by Severity")+
  ylab("Mean Increment %")
  
#scatterplot of DBH
ggplot(RGR_plot, aes(x=species, y=DBH_cm)) +
  geom_point()+
  ggtitle("DBH Scatterplot")

#summarize and plot with error bars 
RGR_summ <- summarySE(RGR_plot, measurevar="increment_percent", groupvars=c("severity","week"))