#setting working directory
setwd("C:/Users/PC/Documents/ThesisData_FoRTE/G-driveUploads_new")

#importing csv
dendro_data <- read.csv("FoRTE_dendro_data.csv")

#set library to dyplr
library(dplyr)

#filter for the D_rep
d_rep <- filter(dendro_data, Replicate == "D") 
#change capital in Replicate if you overwrite last csv

#set library to ggplot
library(ggplot2)

head(d_rep)

ggplot(data = d_rep, aes(x = Date, y = Subplot))
  geom_point(color = "#00AFBB", size = 2)
  
#plotting a time series 
  data(d_rep_condensed, package="ggplot2")  # init data
  d_rep_condensed <- data.frame(d_rep_condensed)  # convert to dataframe
  ggplot(d_rep_condensed) + geom_point(aes(x=DOY, y=increment_., color="subplot")) + scale_color_discrete(name="Legend") + labs(title="TheDRep")
