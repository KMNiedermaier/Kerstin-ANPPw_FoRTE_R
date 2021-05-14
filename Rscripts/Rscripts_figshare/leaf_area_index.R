# this script takes LAI data collected with hemispherical images, plots the mean 
# end-of-season LAI by disturbance severity and type, and runs an ANOVA to test for
# significant differences 

# load packages 
library(agricolae)
library(gridExtra)
library(grid)

# reformat LAI data and make a histogram of end of season LAI among severities and types 

camera_data <- read.csv("data/FoRTE_LAI.csv")

# import data
# camera_data <- fd_hemi_camera()
# camera_data <- distinct(camera_data)

#select columns of interest and rename columns to match other format
LAI <- select(camera_data, SubplotID, severity, treatment, Date, Day, LAI_cam)
names(LAI) <- c("subplot", "severity", "treatment", "date", "DOY", "LAI")
LAI <- arrange(LAI, subplot)

# filter for end of season measurements 
LAI_end <- LAI %>% 
  filter(DOY == 214 | DOY == 213 | DOY == 215) %>% 
  arrange(subplot) 

# create df for severity plot and for treatment (type) plot 
LAI_severity <- LAI_end %>% 
  group_by(severity) %>% 
  summarize(LAI_mean = mean(LAI), SE = std.error(LAI))

LAI_severity$severity <- as.numeric(LAI_severity$severity) 


LAI_type <- LAI_end %>% 
  group_by(treatment) %>% 
  summarise(LAI_mean = mean(LAI), SE = std.error(LAI))

# recode type to read bottom up and top down 
LAI_type$treatment <- recode(LAI_type$treatment, B = "Bottum-Up")
LAI_type$treatment <- recode(LAI_type$treatment, T = "Top-Down")


# generate some figures for LAI_severity and LAI_type
forte <- c("#000000", "#009E73", "#0072B2", "#D55E00")
#LAI_severity$severity <- as.factor(LAI_severity$severity)

LAI_sev <- ggplot(LAI_severity, aes(x = factor(severity), y = LAI_mean, fill = factor(severity))) +
  geom_bar(stat = "identity", position = "dodge", width = 0.4) +
  theme_bw() +
  theme(axis.text.x= element_text(size = 17), axis.text.y= element_text(size=17), 
        axis.title.x = element_text(size=19), axis.title.y  = element_text(size = 19), 
        legend.title = element_blank(), legend.position = c(0.5, 1.06), 
        legend.text = element_text(size = 17), legend.direction = "horizontal", 
        panel.grid = element_blank(), axis.line = element_line(size = 0.5), 
        plot.margin = margin(1.5,1,0.5,1, "cm")) +
  scale_y_continuous(limits = c(0,3.6), sec.axis = sec_axis(~ .,labels = NULL)) +
  scale_x_discrete(position = "top") +
  scale_x_discrete(position = "bottom")+
  scale_fill_manual(aesthetics = "fill", values = forte)+
  guides(fill = FALSE) +
  labs(x = "Disturbance Severity (%)",
       y = "LAI") +
  theme(axis.title = element_text(size = 10)) +
  geom_errorbar(position = position_dodge(width = 0.4), 
                aes(ymin = LAI_mean-SE, ymax = LAI_mean+SE), width = 0.1, color = "grey32") +
  annotate(geom="text", x = 1, y = 3.5, label = c("a"), color="black", size = 7) +
  annotate(geom="text", x = 2, y = 3.6, label = c("a"), color="black", size = 7) +
  annotate(geom="text", x = 3, y = 3.15, label = c("ab"), color="black", size = 7) +
  annotate(geom="text", x = 4, y = 2.9, label = c("b"), color="black", size = 7) +
  annotate(geom="text", x = 0.6, y = 3.5, label = c("A"), color="black", size = 5) 

# and for disturbance type 

LAI_tx <- ggplot(LAI_type, aes(x = factor(treatment), y = LAI_mean, fill = factor(treatment))) +
  geom_bar(stat = "identity", position = "dodge", width = 0.4) +
  theme_bw() +
  theme(axis.text.x= element_text(size = 17), axis.text.y= element_text(size=17), 
        axis.title.x = element_text(size=19), axis.title.y  = element_text(size = 19), 
        legend.title = element_blank(), legend.position = "blank", 
        legend.text = element_text(size = 17), panel.grid = element_blank(),
        axis.line = element_line(size = 0.5), plot.margin = margin(0,1,1.5,1, "cm")) +
  scale_y_continuous(limits = c(0,3.6), sec.axis = sec_axis(~ .,labels = NULL)) +
  scale_fill_manual(values = c("olivedrab4", "firebrick4"), 
                    labels = c("DOY: 186-190", "DOY: 213-215"))+
  labs(x = "Disturbance Type",
       y = "LAI") +
  theme(axis.title = element_text(size = 10)) +
  geom_errorbar(position = position_dodge(0.4), aes(ymin = LAI_mean-SE, ymax = LAI_mean+SE), 
                width = 0.1, color = "grey32") +
  annotate(geom="text", x = 2.5, y = 3.5, label = c("n.s."),
           color="black", size = 7) +
  annotate(geom="text", x = 0.5, y = 3.5, label = c("B"), color="black", size = 5) 

# paste plots together 
g2 <- ggplotGrob(LAI_sev)
g3 <- ggplotGrob(LAI_tx)
g <- rbind(g2, g3, size = "first")
g$widths <- unit.pmax(g2$widths, g3$widths)
grid.newpage()
grid.draw(g)

ggsave(filename = "figures/LAI_tx_sev.jpeg", height = 11, width = 9, units = "in", g)

# run some stats on LAI by type and severity 
LAI_end$severity <- as.factor(LAI_end$severity)
LAI_end$treatment <- as.factor(LAI_end$treatment)

LAI$severity <- as.factor(LAI$severity)
LAI$treatment <- as.factor(LAI$treatment)

# create df for early and late season LAI
LAI_july <- filter(LAI_end, DOY == 186 | DOY == 189 | DOY == 190)
LAI_aug <- filter(LAI_end, DOY == 213 | DOY == 214 | DOY == 215)

# run an ANOVA to see if LAI differs significantly in august among severities and types 
LAI_stats <- aov(LAI ~ severity + treatment + severity*treatment, data = LAI_end)
summary(LAI_stats)
LSD_ouput <- with(LAI_end, LSD.test(LAI, severity, 157, 0.908, console = TRUE))

