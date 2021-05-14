#scale to NPP by subplot, week, species, AND DBH (where necessary)
#runs regression for each subplot, week, and species AND returns object with 
#coefficients
regressions <- RGR_plot_regression %>% 
  group_by(subplot, week, species) %>% 
  do(mod = summary(lm(RGR ~ DBH_cm, data = .))) #%>% #-1 removes the intercept
  
final_output <- matrix(nrow= 2321, ncol=4, dimnames = list(c(), c("uniqueID", "pvalue", "slope", "intercept")))
regressions$uniqueID <- paste(regressions$subplot, regressions$week, 
                              regressions$species, sep = "_")

for(i in 1:2321) {
  
  final_output[i, 4] <- regressions$mod[[i]][4]$coefficients[1]
  final_output[i, 3] <- regressions$mod[[i]][4]$coefficients[2]
  final_output[i, 2] <- regressions$mod[[i]][4]$coefficients[8]
  final_output[i, 1] <- regressions$uniqueID[[i]] #may need [,i] or [i,]
  
} 

#convert matrix to data frame 
final_output_df <- as.data.frame(final_output) 

#remove NA's: NEED TO FIGURE OUT WHY NA'S VS NaN'S 
final_output_df_NAomit <-   filter(final_output_df, pvalue != "NaN", !is.na(pvalue)) 

#coerce pvalue, slope, and intercept into character and then numeric
final_output_df_NAomit$pvalue <-  as.numeric(as.character(final_output_df_NAomit$pvalue))
final_output_df_NAomit$slope <-  as.numeric(as.character(final_output_df_NAomit$slope))
final_output_df_NAomit$intercept <-  as.numeric(as.character(final_output_df_NAomit$intercept))

#filter for pvalue's < 0.05; subplot, weeks, and species where DBH is having a sig effect 
sig_DBH_effect <- filter(final_output_df_NAomit, pvalue < 0.05)


#get species count & total subplot delta biomass per week 
inc_days <- RGR_plot %>% 
  group_by(uniqueID) %>% 
  summarise(inc_days= mean(inc_days)) # cheater way to get inc days with subplot and week  
  
weekly_RGR_species <- RGR_plot %>%
  group_by(subplot, week, species) %>% 
  #increment % is not a raw increment; think I need raw delta biomass to scale 
  #to plot level delta biomass 
  #summarise(mean_delta_biomass= mean(increment_kg)) %>% 
  summarise(mean_RGR= mean(RGR), na.rm = TRUE) %>% 
  full_join(species_count) # need to apply overall subplot means to NA values 

#create uniqueID for subplot, week, species 
weekly_RGR_species$uniqueID <- paste(weekly_RGR_species$subplot, weekly_RGR_species$week,
                                      weekly_RGR_species$species, sep = "_")

weekly_RGR_subplot <- weekly_RGR_species %>% 
  summarise(mean_subplot_deltaB= sum(mean_species_deltaB)) %>% 
  mutate(biomass_per_ha= mean_subplot_deltaB*10) %>% 
  mutate(kgC_per_ha= biomass_per_ha*0.48) %>% 
  group_by(subplot, week) %>% 
  right_join(inc_days) %>% 
  mutate(kgC_ha_day= kgC_per_ha/inc_days)
  
###############################################################################
  


#########################################################################
#scraps

#call a pvalue from a certain row (subplot, week, species)
regressions$mod[[12]][4]$coefficients[8] # n is the row number i want to call








summarise(rsq = summary(reg_coeff)$r.squared) #getting error for the r.squared thing 
do(as.data.frame(coef(.$mod)))

ggplot(weekly_RGR_subplot, aes(x=week, y=kgC_ha_day, color=subplot))+
  geom_point()
