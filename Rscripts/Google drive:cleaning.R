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

