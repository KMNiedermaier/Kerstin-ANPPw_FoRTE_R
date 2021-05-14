# this script is essentially an ingestion script for the raw seedling and sapling 
# data. The ruler used for measuring stem height had a half cm space before tick marks
# began, thus requiring all stems measured with the ruler to have an additional 0.5cm
# added to their height. The criteria for these stems were >0 but =<30 cm height.

# In the bottom half of this script, I take the seedling_sapling_2019.csv and change
# the format to match the 2020 data. the 2019 data did not need the 0.5cm adjustment. 

library(dplyr)
library(fortedata)
library(googledrive)

# import the raw seedling and sapling file 
ss_2020 <- read.csv("data/seedling_sapling_2020.csv", na.strings = c("", "NA"))

# This some renaming used to adjust column headers then using case_when() to add 0.5cm
# to stems that fit the criteria for being measured with the ruler. 
ss_2020 <-  ss_2020 %>% 
  rename(height_budscar = height_budscar_cm, height_total = height_total_cm) %>% 
  mutate(height_budscar_cm = case_when(
    height_budscar <= 30 & height_budscar > 0 ~ height_budscar + 0.5,
    height_budscar > 30 ~ height_budscar,
    height_budscar == 0 ~ height_budscar
  )) %>% 
  mutate(height_total_cm = case_when(
    height_total <= 30 & height_total > 0 ~ height_total + 0.5,
    height_total > 30 ~ height_total,
    height_total == 0 ~ height_total
  )) %>% 
  select(subplot, veg_plot, species, base_D, height_budscar_cm, height_total_cm, date) %>% 
  rename(basal_diameter_cm = base_D, nested_subplot = veg_plot)

# changed subplot to subplot_id
names(ss_2020)[names(ss_2020) == "subplot"] <- "subplot_id"

# write a csv and save it to tdata folder (would like to find a way to save this 
# as a temp file and then upload )
write.csv(ss_2020, "data/fd_seedling_sapling_2020.csv", row.names = FALSE)

# upload new csv to google drive 

# first I need to find the id 
# find the folder pathway on mydrive 
x <- drive_find(type = "folder", pattern = "seedling_sapling_", n_max = 6)
y <- x[[1,2]]

# now upload using the file id we located above 
drive_upload(
  "data/fd_seedling_sapling_2020.csv",
  name = "fd_seedling_sapling_2020",
  path = as_id(y),
  overwrite = TRUE
)

###################################################################################
# ingest the 2019 seedlin_sapling data na reformat to match 2020

ss_2019 <- read.csv("data/seedling_sapling.csv", na.strings = c("","NA"))


# first selct the wanted columns, then rename height_2018/2019 
ss_2019 <- ss_2019 %>% 
  select(subplot, vegplot_direction, species, baseD_cm, height_2018, height_2019, date) %>% 
  rename(nested_subplot = vegplot_direction, basal_diameter_cm = baseD_cm, 
         height_budscar_cm = height_2018, height_total_cm = height_2019) 

# changed subplot to subplot_id
names(ss_2019)[names(ss_2019) == "subplot"] <- "subplot_id"

# write a csv and save it to tdata folder for upload to fortedata/ google drive 
write.csv(ss_2019, "data/fd_seedling_sapling_2019.csv", row.names = FALSE)

drive_upload(
  "data/fd_seedling_sapling_2019.csv",
  name = "fd_seedling_sapling_2019",
  path = as_id(y),
  overwrite = TRUE
)
