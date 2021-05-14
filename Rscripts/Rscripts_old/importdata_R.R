#opening a file
library(readxl) 

dendro_data <- read_excel("FoRTE_plot_dendro_data.xlsx", sheet=2)

#how many rows?
nrow(dendro_data)

#look at first 5 rows
head(dendro_data)

#Evaluate data structure 
str(dendro_data)

#focus on a row, cell, or column
dendro_data[3,1]

dendro_data$`Bands (cm)`[345]

#let's try and plot this shiiit 
stripchart(dendro_data$`Increment (%)`)

# installing/loading the latest installr package:
install.packages("installr"); library(installr) # install+load installr

updateR() # updating R.

#scraps from increment scripts; un-piped version of filter, mutate, join 
all_reps_filtered <- all_reps %>%
  filter(bandread_. == 1) 
biomass_1 <- mutate(all_reps_filtered, biomass= a*DBH_cm^b)
all_reps_biomass <- left_join(all_reps, biomass_1)
                