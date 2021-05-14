# This is a cleaning script for all of the double dendrobanded trees in our sample. 
# In an attempt to track potential stem swelling near the girdle injury, we added a 
# second dendrometer band 50cm above the orginal band to a sub sample of trees. 
# These data can be used to calculate % differences in top and bottom dendrobands 
# and estimate potential stem swelling. 

library(dplyr)
library(googledrive)

# import date
dbl_band_2019 <- read.csv("data/DBL_dendrobands_2019.csv", na.strings = c("", "NA"))
dbl_band_2020 <- read.csv("data/DBL_dendrobands_2020.csv", na.strings = c("", "NA"))

# first take the August 2019 measurements off the 2020 sheet
dbl_band_2020 <-  filter(dbl_band_2020, date != "2019-08-07", date != "2019-08-08")

# get rid of unwanted columns 
dbl_band_2019 <- dbl_band_2019[c("tag", "species", "bottom..in.", "top",
                                 "date", "notes")]

dbl_band_2020 <- dbl_band_2020[c("tag", "species", "bottom..in.", "top",
                                 "date", "notes")]

# change the bottom and top band column names 
names(dbl_band_2019)[names(dbl_band_2019) == "bottom..in."] <- "band_in_bottom"
names(dbl_band_2019)[names(dbl_band_2019) == "top"] <- "band_in_top"

names(dbl_band_2020)[names(dbl_band_2020) == "bottom..in."] <- "band_in_bottom"
names(dbl_band_2020)[names(dbl_band_2020) == "top"] <- "band_in_top"

# change the dates to date format 
dbl_band_2019$date <- as.Date(dbl_band_2019$date, "%Y-%m-%d")
dbl_band_2020$date <- as.Date(dbl_band_2020$date, "%Y-%m-%d")

# now arrange 2020 data by tag 
dbl_band_2020 <- arrange(dbl_band_2020, tag)

#####################################################################################
# alright now they are properly formated and ready to go onto fortedata. 
# upload to google drive 

# write the csv
write.csv(dbl_band_2019, "data/fd_double_dendroband_2019.csv", row.names = FALSE, quote = FALSE)
write.csv(dbl_band_2020, "data/fd_double_dendroband_2020.csv", row.names = FALSE, quote = FALSE)

# find the folder pathway on mydrive 
x <- drive_find(type = "folder", pattern = "canopy_dendroband", n_max = 6)
y <- x[[1,2]]

# now upload using the file id we located above 
drive_upload(
  "data/fd_double_dendroband_2019.csv",
  name = "fd_double_dendroband_2019",
  path = as_id(y),
  overwrite = TRUE
)

# do the same for 2020
# now upload using the file id we located above 
drive_upload(
  "data/fd_double_dendroband_2020.csv",
  name = "fd_double_dendroband_2020",
  path = as_id(y),
  overwrite = TRUE
)

