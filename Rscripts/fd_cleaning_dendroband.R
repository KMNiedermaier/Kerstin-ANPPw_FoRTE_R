# this script is an ingestion script for the raw dendroband data.

library(googledrive)

# Bring in the csv's and make the format matching. Then upload them to 
# google drive

# This brings in raw data from the google drive. First I need to id the files I want.
# An easy way to do this is as follows. q = "starred = true" gets me the file I want quickly,
# but first I need to navigate to the file in My Drive, right click, and select "Add to Starred"

# Let's do this for the 2019 data 
file_id <- drive_find(type = "csv", pattern = "raw_canopy_dendroband_2019", q = "starred = true")
dendro_2019_id <- file_id[[1,2]]

# Now for the 2020
file_id <- drive_find(type = "csv", pattern = "raw_canopy_dendroband_2020", q = "starred = true")
dendro_2020_id <- file_id[[1,2]]

# now I need to download this file and save it to my local data folder. Since I probably 
# (in my case definitely...) already have a version of this file, I'll name it something
# different. In this case, I'll use gd for google drive

# alright, first for the 2019
drive_download(
  as_id(dendro_2019_id), 
  path = "data/gd_canopy_dendroband_2019.csv",
  overwrite = TRUE)

# and tnow the 2020
drive_download(
  as_id(dendro_2020_id), 
  path = "data/gd_canopy_dendroband_2020.csv",
  overwrite = TRUE)

# This brings in the data from my local machine
dendro_2019 <- read.csv("data/gd_canopy_dendroband_2019.csv", na.strings = c("", "NA"))
dendro_2020 <- read.csv("data/gd_canopy_dendroband_2020.csv", na.strings = c("", "NA"))

# Now it's time to clean 'em up!! 2019 is a mess, so I'll start there 

# now rename columns to match other forte data products
names(dendro_2019)[names(dendro_2019) == "subplot"] <- "subplot_id"
names(dendro_2019)[names(dendro_2019) == "bands_in"] <- "band_in"

names(dendro_2020)[names(dendro_2020) == "subplot"] <- "subplot_id"
names(dendro_2020)[names(dendro_2020) == "bands_in"] <- "band_in"

# quick date reformat for 2019
dendro_2019$date <- as.Date(dendro_2019$date, "%m/%d/%Y")

# now drop unwanted columns and reorder 
dendro_2019 <- dendro_2019[c("subplot_id", "tag", "species", "band_in", 
                     "date", "notes")]

dendro_2020 <- dendro_2020[c("subplot_id", "tag", "species", "band_in", 
                             "date", "notes")]

#####################################################################################
# alright now they are properly formated and ready to go onto fortedata. 
# upload to google drive 

# first, write csv's
write.csv(dendro_2019, "data/fd_canopy_dendroband_2019.csv", row.names = FALSE)
write.csv(dendro_2020, "data/fd_canopy_dendroband_2020.csv", row.names = FALSE)

# find the folder pathway on mydrive 
x <- drive_find(type = "folder", pattern = "canopy_dendroband", n_max = 6)
y <- x[[1,2]]

# now upload using the file id we located above 
drive_upload(
  "data/fd_canopy_dendroband_2019.csv",
  name = "fd_canopy_dendroband_2019",
  path = as_id(y),
  overwrite = TRUE
)

# do the same for 2020
# now upload using the file id we located above 
drive_upload(
  "data/fd_canopy_dendroband_2020.csv",
  name = "fd_canopy_dendroband_2020",
  path = as_id(y),
  overwrite = TRUE
)
