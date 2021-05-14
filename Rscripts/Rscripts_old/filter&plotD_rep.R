#setting working directory
setwd("C:/Users/PC/Documents/ThesisData_FoRTE/G-driveUploads_new")

#importing csv
dendro_data <- read.csv("FoRTE_dendro_data.csv")

#set library to dyplr
library(dplyr)

#filter for the D_rep
d_rep <- filter(dendro_data, replicate == "D") 
#change capital in Replicate if you overwrite last csv

#selected only for columns that I will use 
d_rep_condensed <- select(d_rep, subplot, replicate, tag, species, DBH_cm, 
       increment_.,DOY, bandread_.)

#summarize avg increment by DOY and subplot, species, etc
d_rep_condensed %>%
  group_by(DOY, subplot) %>%
  summarise(increment_mean = mean(increment_.))

#plot this shiiiiit
library(ggplot2)

ggplot(d_rep_condensed, aes(x=DOY, y=increment_., 
      color=subplot))+geom_point()

######original filter-mutate-join for biomass accumulation#########

#selected only for columns that I will use 
all_reps <- select(dendro_data, severity, fate, treatment, week, subplot, tag, replicate, 
                   species, DBH_cm, a, b, bands_cm,DOY, bandread_.)

#starting biomass and band read; 
all_reps_filtered <- all_reps %>%
  filter(bandread_. == 1) %>%
  mutate(biomass_a= a*DBH_cm^b) %>% #initial biomass
  mutate(bandread_i= bands_cm) %>% #november band read 
  right_join(all_reps) %>% #joining filtered table w/new biomass col
  group_by(tag) %>% #group tags and fill starting biomass
  fill(biomass_a) %>%
  fill(bandread_i) %>%
  mutate(biomass_b_kg= a*(DBH_cm+(bands_cm-bandread_i))^b) %>%
  mutate(increment_kg= biomass_b_kg-lag(biomass_b_kg, default = first(biomass_b_kg))) %>%
  mutate(increment_percent= increment_kg/biomass_b_kg)

#select for RGR, tag, date, subplot, species, DBH
RGR_plot <- all_reps_filtered %>%
  select(subplot, severity, tag, fate, replicate, species, DBH_cm, DOY, increment_percent) %>%
  drop_na() %>%
  group_by(tag) %>%
  mutate(bandread=row_number()) %>%
  ungroup()

### model attempts 
library(lmerTest)
fit <- lmer(increment_percent ~ severity * treatment + (1 | subplot), data = RGR_plot)

RGR_plot$week <- as.factor(RGR_plot$week)
RGR_plot$severity <- as.factor(RGR_plot$severity)

library(agricolae)
modeltry1 <- with(RGR_plot, ssp.plot(week, replicate, severity, treatment, increment_percent))

#deleting outlier rows, but must make sure rows do not shift
all_reps <- all_reps[-c(8169),]
