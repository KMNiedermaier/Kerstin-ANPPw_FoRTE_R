fit <- aov(increment_percent ~ species + DBH_cm, data=RGR_plot)
summary(fit)

fit_I <- aov(increment_percent ~ species, data=RGR_plot)
summary(fit_I)

fit_II <- aov(increment_percent ~ species*DBH_cm + DOY, data=RGR_plot)
summary(fit_II)
#FOR TIME AS CATEGORICAL -- run factor(week) instead of DOY
#NEED TO: loop aov for every subplot
library(stargazer) 

analysis <- by(RGR_plot, RGR_plot$subplot,
   function(x) summary(aov(increment_percent ~ species*DBH_cm + factor(week), 
                           data = x)))
analysis[1:32]

stargazer(analysis[1:32], type = "text")  

