# install agricolae (use this for the LSD test, but could maybe just figure out the
# stats version of an LSD); Knitr and broom are just for visualizing and exporting 
# anova tables 
library(agricolae)
library(knitr)
library(broom)

#run some stats on fate contributions to anppw
NPP_alltime_fate$severity <- as.factor(NPP_alltime_fate$severity)
NPP_alltime_fate$treatment <- as.factor(NPP_alltime_fate$treatment)
NPP_alltime_fate$replicate <- as.factor(NPP_alltime_fate$replicate)
NPP_alltime_fate$fate <- as.factor(NPP_alltime_fate$fate)
NPP_alltime_fate$year <- as.factor(NPP_alltime_fate$year)

# run the mixed anova model on ANPPw by severity, type, year, and fate 
aov.fate <- aov(kgC_ha_yr ~ replicate + severity*year*fate +
                  treatment*year*fate +
                  Error(replicate:(treatment*fate*year)), data = NPP_alltime_fate)
summary(aov.fate)

# This chunk prints a nice pretty anova output and saves it as a df
aov.fate <- tidy(aov.fate)
options(knitr.kable.NA = '')
kable(aov.fate[-14, -1], digits = 3, format = "pandoc")
aov.fate.df <- data.frame(aov.fate)

# LSD test for severity:fate
LSD_ouput <- with(NPP_alltime_fate, LSD.test(kgC_ha_yr, severity:fate:year, 84, 477227, 
                                     console = TRUE, alpha = 0.1))

# LSD test for year:fate
LSD_ouput <- with(NPP_alltime_fate, LSD.test(kgC_ha_yr, fate:year, 3, 27644, 
                                             console = TRUE, alpha = 0.05))

# LSD test for fate:treatment(type):year
LSD_ouput <- with(NPP_alltime_fate, LSD.test(kgC_ha_yr, fate:treatment:year, 3, 951343, 
                                             console = TRUE, alpha = 0.05))

# get post-hoc output into a df
output_df <- data.frame(tibble::rownames_to_column(LSD_ouput$groups))

# seperate into severity and date columns 
output_df <- transform(output_df, severity = substr(rowname, 1, 2), 
                       fate = substr(rowname, 3, 7), year = substr(rowname, 8,12))

# recode output_df 
output_df$severity <- recode(output_df$severity, "0:" = "0")
output_df$fate <- recode(output_df$fate, ":kill" = "kill", ":live" = "live", 
                         "kill:" = "kill", "live:" = "live") 
output_df$year <- recode(output_df$year, ":2020" = "2020",":2019" = "2019")

# arrange output df by severity and fate
output_df <- arrange(output_df, severity, year)


#####################################################################################
# write a model for % of total ANPPw 
# essentially just change the response variable to perc_total (% of total ANPPw)
aov.fate <- aov(perc_total ~ replicate + severity*year*fate +
                  treatment*year*fate +
                  Error(replicate:(treatment*fate*year)), data = NPP_alltime_fate)
summary(aov.fate)

# LSD test for severity:fate
LSD_ouput <- with(NPP_alltime_fate, LSD.test(perc_total, severity:fate:year, 84, 0.019, 
                                             console = TRUE, alpha = 0.05))

# get post-hoc output into a df
output_df <- data.frame(tibble::rownames_to_column(LSD_ouput$groups))

# seperate into severity and date columns 
output_df <- transform(output_df, severity = substr(rowname, 1, 2), 
                       fate = substr(rowname, 3, 7), year = substr(rowname, 8,12))

# recode output_df 
output_df$severity <- recode(output_df$severity, "0:" = "0")
output_df$fate <- recode(output_df$fate, ":kill" = "kill", ":live" = "live", 
                         "kill:" = "kill", "live:" = "live") 
output_df$year <- recode(output_df$year, ":2020" = "2020",":2019" = "2019")

# arrange output df by severity and fate
output_df <- arrange(output_df, severity, year)

# LSD test for severity:fate
LSD_ouput <- with(NPP_alltime_fate, LSD.test(perc_total, severity:year, 84, 0.019, 
                                             console = TRUE, alpha = 0.05))
