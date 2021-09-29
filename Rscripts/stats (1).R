#### Model Runs
# install agricolae (use this for the LSD test, but could maybe just figure out the
# stats version of an LSD); Knitr and broom are just for visualizing and exporting 
# anova tables 
library(agricolae)
library(knitr)
library(broom)

# clean up my NPP_final (from NPP_2019) df for model run. Need to run all variables 
# as factors 
NPP_final$severity <- as.factor(NPP_final$severity)
NPP_final$replicate <- as.factor(NPP_final$replicate)
NPP_final$treatment <- as.factor(NPP_final$treatment)
NPP_final$week <- as.factor(NPP_final$week)

# recode 0 to 0.00 to help me with the post-hoc output viualization
NPP_final$severity <- recode_factor(NPP_final$severity, "0" = "0.00")

# Need to first ungoup my df, select the important columns, and filter
# out the 2020 records (this is 2019 NPP!)
model_run <- NPP_final %>% 
  filter(week!="2020-07-26") %>% 
  ungroup() %>% 
  select(replicate, week, severity, treatment, kgC_per_ha_day)

# Let's run an aov, assign it an object, and then print the summary
aov.NPP<- aov(kgC_per_ha_day ~ replicate + week*severity + treatment*week +
                Error(replicate:(week*treatment)),
              data = model_run)
summary(aov.NPP)

# run LSD post hoc with week:severity ; the numbers are error df's and MSE
LSD_ouput <- with(model_run, LSD.test(kgC_per_ha_day, week:severity, 336, 75.7, 
                                      console = TRUE, alpha = 0.1))

# get post-hoc output into a df
output_df <- data.frame(tibble::rownames_to_column(LSD_ouput$groups))

# seperate into severity and date columns 
output_df <- transform(output_df, severity = substr(rowname, 12, 15), 
                       date = substr(rowname, 1, 10))

# sort by date and severity so that we can look which severities differ within a week 
output_df <- output_df %>% 
  arrange(date, severity)


# run LSD post hoc with severity (error df's and MSE)
LSD_ouput <- with(model_run, LSD.test(kgC_per_ha_day, severity, 336, 75.7, 
                                      console = TRUE, alpha = 0.1))

# This chunk prints a nice pretty anova output, and saves it as a df 
aov.NPP <- tidy(aov.NPP)
options(knitr.kable.NA = '')
kable(aov.NPP[-14, -1], digits = 3, format = "pandoc", caption = "Daily ANPPw ANOVA table")
aov.df <- data.frame(aov.NPP)

# now let's write a csv from the df to copy nd paste numbers into word table 
write.csv(aov.df, "figures/aov_dailyANPPw.csv")

#######################################################################################
########### ANNUAL ANPPw Split Plot Model (without time)
# this is basically the same exact code as above, except run on the all_strata df

# clean up my all_strata df for model run. Need to run all variables as factors 
all_strata$severity <- as.factor(all_strata$severity)
all_strata$replicate <- as.factor(all_strata$replicate)
all_strata$treatment <- as.factor(all_strata$treatment)

#aov for all canopy strata ANPPw 
aov.allstrata<- aov(annual_NPP ~ replicate + severity + treatment +
                Error(replicate:treatment),
              data = all_strata)
summary(aov.allstrata)

# Make a nice pretty table to export 
# This chunk prints a nice pretty anova output
aov.allstrata <- tidy(aov.allstrata)
options(knitr.kable.NA = '')
kable(aov.allstrata[-14, -1], digits = 3, format = "pandoc")
aov.allstrata.df <- data.frame(aov.allstrata)

# run LSD post hoc for all strata NPP by severity  (error df's and error MSE)
LSD_ouput <- with(all_strata, LSD.test(annual_NPP, severity, 21, 546153, 
                                       console = TRUE, alpha = 0.1))

# now I can run an aov for each canopy strata seperately 
#aov for canopy  ANPPw
aov.canopyNPP<- aov(NPP_canopy ~ replicate + severity + treatment +
                      Error(replicate:treatment),
                    data = all_strata)
summary(aov.canopyNPP)

# post hoc for canopy ANPPw by disturbance type 
LSD_ouput <- with(all_strata, LSD.test(NPP_canopy, treatment, 3, 1200627, 
                                       console = TRUE, alpha = 0.1))

#aov for subcanopy ANPPw
aov.subcanNPP<- aov(NPP_subcan ~ replicate + severity + treatment +
                      Error(replicate:treatment),
                    data = all_strata)
summary(aov.subcanNPP)

#aov for seedling/saplings  NPP
aov.seedlings<- aov(NPP_seedlings ~ replicate + severity + treatment +
                      Error(replicate:treatment),
                    data = all_strata)
summary(aov.seedlings)

# run LSD post hoc for seedling/sapling replicate (error df's and error MSE)
LSD_ouput <- with(all_strata, LSD.test(NPP_seedlings, replicate, 3, 56.6, 
                                      console = TRUE, alpha = 0.1))

####################################################################################
###########Annual NPP split plot model WITH PFT 
# clean up my PFT_aov df for model run. Need to run all variables as factors 
PFT_aov$severity <- as.factor(PFT_aov$severity)
PFT_aov$replicate <- as.factor(PFT_aov$replicate)
PFT_aov$treatment <- as.factor(PFT_aov$treatment)
PFT_aov$PFT <- as.factor(PFT_aov$PFT)

# run an aov for annual ANPPw (all strata) split by PFT
aov.PFT <- aov(annual_NPP ~ replicate + severity*PFT + PFT*treatment +
               Error(replicate:(treatment*PFT)), data = PFT_aov)
summary(aov.PFT)

# LSD test for severity:PFT
LSD_ouput <- with(PFT_aov, LSD.test(annual_NPP, severity:PFT, 42, 494320, 
                                       console = TRUE, alpha = 0.1))

# This chunk prints a nice pretty anova output
aov.PFT <- tidy(aov.PFT)
options(knitr.kable.NA = '')
kable(aov.PFT[-14, -1], digits = 3, format = "pandoc")
aov.PFT.df <- data.frame(aov.PFT)

# write a csv that can be copy and pasted into a word doc table 
write.csv(aov.PFT.df, "figures/aov_PFT.csv")
####################################################################################
###########Annual NPP split plot model WITH fate (girdled/ungirdled)

#run some stats on fate contributions to anppw
NPP_fate$severity <- as.factor(NPP_fate$severity)
NPP_fate$treatment <- as.factor(NPP_fate$treatment)
NPP_fate$replicate <- as.factor(NPP_fate$relpicate)

aov.fate <- aov(kgC_ha_yr ~ replicate + severity*fate + fate*treatment +
                 Error(replicate:(treatment*fate)), data = NPP_fate)
summary(aov.fate)

# LSD test for severity:fate
LSD_ouput <- with(NPP_fate, LSD.test(kgC_ha_yr, severity:fate, 42, 326376, 
                                     console = TRUE, alpha = 0.1))

# LSD test for treatment:fate
LSD_ouput <- with(NPP_fate, LSD.test(kgC_ha_yr, fate:treatment, 3, 305328, 
                                     console = TRUE, alpha = 0.1))

# This chunk prints a nice pretty anova output and saves it as a df
aov.fate <- tidy(aov.fate)
options(knitr.kable.NA = '')
kable(aov.fate[-14, -1], digits = 3, format = "pandoc")
aov.fate.df <- data.frame(aov.fate)

# export the df to a csv to copy and paste into word doc table 
write.csv(aov.fate.df, "figures/aov_fate.csv")

####################################################################################
########### End of season (August) LAI by severity and type
# run some stats on LAI by type and severity 
LAI_end$severity <- as.factor(LAI_end$severity)
LAI_end$treatment <- as.factor(LAI_end$treatment)
LAI_end$replicate <- as.factor(LAI_end$replicate)


aov.LAI <- aov(mean_LAI ~ replicate + severity + treatment + 
                   Error(replicate:treatment), data = LAI_end)
summary(aov.LAI)
LSD_ouput <- with(LAI_end, LSD.test(mean_LAI, severity, 21, 0.3545, console = TRUE))

#####################################################################################
#########Live-Kill Stats 2020
#run some stats on fate contributions to anppw
NPP_fate_2020$severity <- as.factor(NPP_fate_2020$severity)
NPP_fate_2020$treatment <- as.factor(NPP_fate_2020$treatment)
NPP_fate_2020$replicate <- as.factor(NPP_fate_2020$relpicate)

aov.fate.2020 <- aov(kgC_ha_yr ~ replicate + severity*fate + fate*treatment +
                  Error(replicate:(treatment*fate)), data = NPP_fate_2020)
summary(aov.fate.2020)

# LSD test for severity:fate
LSD_ouput <- with(NPP_fate_2020, LSD.test(kgC_ha_yr, severity:fate, 42, 304840, 
                                     console = TRUE, alpha = 0.05))

# LSD test for treatment:fate
LSD_ouput <- with(NPP_fate_2020, LSD.test(kgC_ha_yr, fate:treatment, 3, 308256, 
                                     console = TRUE, alpha = 0.05))

#####################################################################################
####### Canopy ANPPw by severity, year, type (treatment)
NPP_alltime$severity <- as.factor(NPP_alltime$severity)
NPP_alltime$treatment <- as.factor(NPP_alltime$treatment)
NPP_alltime$replicate <- as.factor(NPP_alltime$replicate)
NPP_alltime$year <- as.factor(NPP_alltime$year)

aov.cNPP.2020 <- aov(kgC_ha_yr ~ replicate + severity*year + year*treatment +
                       Error(replicate:(treatment*year)), data = NPP_alltime)
summary(aov.cNPP.2020)

# LSD test for severity:year
LSD_ouput <- with(NPP_alltime, LSD.test(kgC_ha_yr, severity:year, 42, 811609, 
                                           console = TRUE, alpha = 0.05))

####### Subcanopy ANPPw by severity, year, type (treatment)
NPP_sc_alltime$severity <- as.factor(NPP_sc_alltime$severity)
NPP_sc_alltime$treatment <- as.factor(NPP_sc_alltime$treatment)
NPP_sc_alltime$replicate <- as.factor(NPP_sc_alltime$replicate)
NPP_sc_alltime$year <- as.factor(NPP_sc_alltime$year)

aov.scNPP.2020 <- aov(NPP1 ~ replicate + severity*year + year*treatment +
                      Error(replicate:(treatment*year)), data = NPP_sc_alltime)
summary(aov.scNPP.2020)

# LSD test for severity:year
LSD_ouput <- with(NPP_sc_alltime, LSD.test(NPP1, severity:year, 42, 35920, 
                                          console = TRUE, alpha = 0.05))
