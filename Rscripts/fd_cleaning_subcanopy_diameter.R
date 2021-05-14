# this script is an ingestion script for the subcanopy diameter data, and the 
# the subcanopy stem counts.

library(dplyr)
library(fortedata)
library(googledrive)

# Bring in the subcanopy csv's and make the format matching. Then upload them to 
# google drive

# This brings in raw data from the google drive. First I need to id the files I want.
# An easy way to do this is as follows. q = "starred = true" gets me the file I want quickly,
# but first I need to navigate to the file in My Drive, right click, and select "Add to Starred"

# Let's do this for the 2019 data 
file_id <- drive_find(type = "csv", pattern = "subcanopy_D.csv", q = "starred = true")
sc_2019_id <- file_id[[1,2]]

# Now for the 2020
file_id <- drive_find(type = "csv", pattern = "subcanopy_D_2020.csv", q = "starred = true")
sc_2020_id <- file_id[[1,2]]

# now I need to download this file and save it to my local data folder. Since I probably 
# (in my case definitely...) already have a version of this file, I'll name it something
# different. In this case, I'll use gd for google drive

# alright, first for the 2019
drive_download(
  as_id(sc_2019_id), 
  path = "data/gd_subcanopy_D.csv",
  overwrite = TRUE)

# and tnow the 2020
drive_download(
  as_id(sc_2020_id), 
  path = "data/gd_subcanopy_D_2020.csv",
  overwrite = TRUE)

# This brings in the data from my local machine
sc_2019 <- read.csv("data/gd_subcanopy_D.csv", na.strings = c("", "NA"))
sc_2020 <- read.csv("data/gd_subcanopy_D_2020.csv", na.strings = c("", "NA"))

# since 2020 is more tidy, I will start there and then reformat 2019 to match it 
# First step will be to rename column names to match other fortedata products 

# seperate the uniqueID column into subplot and nested_subplot
# this create a new nested_subplot column and creates a notes column to match 2019
sc_2020 <-  sc_2020 %>% 
  mutate(nested_subplot = substr(uniqueID, 5, 5)) %>% 
  mutate(notes = "NA")

# this deletes the nested subplot ID from the uniqueID column
sc_2020$uniqueID <- gsub('.$', '', sc_2020$uniqueID)

# make all column names lowercase 
names(sc_2020) <- tolower(names(sc_2020))

# now rename columns to match other forte data products
names(sc_2020)[names(sc_2020) == "uniqueid"] <- "subplot_id"
names(sc_2020)[names(sc_2020) == "subplotid"] <- "subplot_id"

# now drop unwanted columns and reorder 
sc_2020 <- sc_2020[c("subplot_id", "nested_subplot", "tag", "species", "dbh_mm", 
                     "date", "notes")]

#####################################################################################
# Now lets work on the 2019 data 

# make all column names lowercase 
names(sc_2019) <- tolower(names(sc_2019))

# change veg_plot to nested_subplot and fix subplotid
names(sc_2019)[names(sc_2019) == "veg_plot"] <- "nested_subplot"
names(sc_2019)[names(sc_2019) == "subplot"] <- "subplot_id"

# select the wanted columns and re order 
sc_2019 <- sc_2019[c("subplot_id", "nested_subplot", "tag", "species", "dbh_mm", 
                     "date", "notes")]

# change both 2019 and 2020 to have to same data classes for each vector; data classes 
# for variables already in fortedata metadata should conform to established classes 

# change 2020 nested to int, 2020 and 2019 date to date, 2020 notes to character
sc_2020$nested_subplot <- as.integer(sc_2020$nested_subplot)
sc_2019$date <- as.Date(sc_2019$date, "%Y-%m-%d")
sc_2020$date <- as.Date(sc_2020$date, "%Y-%m-%d")

#####################################################################################
# alright now they are properly formated and ready to go onto fortedata. 
# upload to google drive 

# first, write csv's
write.csv(sc_2019, "data/fd_subcanopy_diameter_2019.csv", row.names = FALSE)
write.csv(sc_2020, "data/fd_subcanopy_diameter_2020.csv", row.names = FALSE)

# find the folder pathway on mydrive 
x <- drive_find(type = "folder", pattern = "subcanopy_diameter", n_max = 6)
y <- x[[1,2]]

# now upload using the file id we located above 
drive_upload(
  "data/fd_subcanopy_diameter_2019.csv",
  name = "fd_subcanopy_diameter_2019",
  path = as_id(y),
  overwrite = TRUE
)

# do the same for 2020
# now upload using the file id we located above 
drive_upload(
  "data/fd_subcanopy_diameter_2020.csv",
  name = "fd_subcanopy_diameter_2020",
  path = as_id(y),
  overwrite = TRUE
)
########################################################################
# subcanopy stem counts reformatting

stem_counts <- read.csv("data/subcanopy_stemcounts.csv", na.strings = c("", "NA"))

# rename subplot column
names(stem_counts)[names(stem_counts) == "subplot"] <- "subplot_id"

#select wanted column
stem_counts <- stem_counts[c("subplot_id", "species", "count", "plot_area_m2", "date")]

# change date to date class
stem_counts$date <- as.Date(stem_counts$date, "%Y-%m-%d")

# save the new csv
write.csv(stem_counts, "data/fd_subcanopy_stem_count_2019.csv", row.names = FALSE)

# upload to drive
drive_upload(
  "data/fd_subcanopy_stem_count_2019.csv",
  name = "fd_subcanopy_stem_count_2019",
  path = as_id(y),
  overwrite = TRUE
)
