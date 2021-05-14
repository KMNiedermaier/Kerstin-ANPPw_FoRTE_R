#import treatment sheets; all trees in each plot 
a <- read.csv("data/group_a_lai_disturbance_treatment_output.csv")
b <- read.csv("data/group_b_lai_disturbance_treatment_output.csv")
c <- read.csv("data/group_c_lai_disturbance_treatment_output.csv")
d <- read.csv("data/group_d_lai_disturbance_treatment_output.csv")

#bind all rep tables 
all_reps_all_trees <- bind_rows(a, b, c, d)

all_trees <- all_reps_all_trees %>% 
  select(SubplotID, Species, Tag, dbh)

names(all_trees) <- c("subplot", "species", "tag", "DBH_cm")

all_trees_weeks <- all_trees %>% 
  slice(rep(1:n(), each = 16)) %>% 
  group_by(tag) %>% 
  mutate(week = row_number()) %>% 
  full_join(all_reps) %>% 
  group_by(subplot, week) %>% 
  fill(DOY) %>%
  fill(DOY, .direction = c("up")) %>% 
  ungroup()

all_trees_weeks$uniqueID <- paste(all_trees_weeks$subplot, all_trees_weeks$week,
                                  all_trees_weeks$species, sep = "_")


#count species in each subplot 
species_count <- all_reps_all_trees %>%
  group_by(SubplotID) %>%
  count(Species) 

#rename col's to match RGR_plot 
names(species_count) <- c("subplot", "species", "species_count")


