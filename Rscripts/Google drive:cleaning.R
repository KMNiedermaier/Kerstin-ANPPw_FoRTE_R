if (!require("googledrive")) {
  install.packages("googledrive")
}
library(googledrive)

file_id <- drive_find(type = "csv", pattern = "subcanopy_diameter_2021.csv", q = "starred = true")
sc_2021_id <- file_id[[1,2]]

file_id <- drive_find(type = "csv", pattern = "canopy_dendrobands_2021.csv", q = "starred = true")
canopy_2021_id <- file_id[[1,2]]


drive_download(
  as_id(canopy_2021_id), 
  path = "2021_data/gd_canopy_dendrobands_2021.csv",
  overwrite = FALSE)

drive_download(
  as_id(sc_2021_id), 
  path = "2021_data/gd_subcanopy_diameter_2021.csv",
  overwrite = FALSE)

canopy_2021 <- read.csv("2021_data/gd_canopy_dendrobands_2021.csv", na.strings = c("", "NA"))
sc_2021 <- read.csv("2021_data/gd_subcanopy_diameter_2021.csv", na.strings = c("", "NA"))


## cleaning data
# First step will be to rename column names to match other fortedata products 
# seperate the uniqueID column into subplot and nested_subplot
# this creates a new nested_subplot column and creates a notes column to match 2019
# *see my data was kind of messy! My uniqueID column did not match fortedata requirements
#  and I was missing a notes column in 2020
canopy_2021$nested_subplot <- substr(sc_2020$uniqueID, 5, 5)
sc_2020$notes <- "NA"

