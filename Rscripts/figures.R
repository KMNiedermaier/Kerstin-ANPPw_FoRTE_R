##################################################################################
# PLOTS!!

# forte colors decided by Jeff 
forte <- c("#000000", "#009E73", "#0072B2", "#D55E00")

# create text grobs for months to then annotate them onto our figures 
april <- textGrob("April", gp=gpar(fontsize=22))
may <- textGrob("May", gp=gpar(fontsize=22))
june <- textGrob("June", gp=gpar(fontsize=22))
july <- textGrob("July", gp=gpar(fontsize=22))
aug <- textGrob("Aug", gp=gpar(fontsize=22))
sept <- textGrob("Sept", gp=gpar(fontsize=22))
oct <- textGrob("Oct", gp=gpar(fontsize=22))
nov <- textGrob("Nov", gp=gpar(fontsize=22))

# create a simple histogram of cumulative NPP as of 07/26/2020
ggplot(NPP_july_2020, aes(x = factor(severity), y = kgC_ha))+
  geom_bar(stat = "identity", width = 0.5)+
  theme_bw()+
  theme(axis.text.x= element_text(size = 22), axis.text.y= element_text(size=22), 
        axis.title.x = element_text(size=30), axis.title.y  = element_text(size = 30), 
        legend.title = element_blank(), legend.position = c(0.15, 0.9), 
        legend.text = element_text(size = 22), panel.grid = element_blank(),
        axis.line = element_line(size = 0.5), plot.margin = margin(1.5,1,1.5,1, "cm")) +
  labs(x = "Disturbance Severity (%)",
       y = expression(paste("Annual ANPP" [w], " ( ",kgC," ",ha^-1," ",yr^-1,")")))+
  geom_errorbar(aes(ymin = kgC_ha-SE, ymax = kgC_ha+SE), width = 0.1) 

# save the figure to figures folder
ggsave("figures/NPP_july_2020.jpeg", height = 8, width = 11, units = "in", plot = last_plot())

## Disturbance severity time series; includes CANOPY only 
sev_time <- ggplot(severity_time, aes(x=week, y=NPP, color = factor(severity))) +
  scale_colour_manual(values = forte, name = "Disturbance\nSeverity (%)", labels = 
                        c("0", "45", "65", "85")) +
  theme_bw() +
  theme(axis.text.x= element_blank(), axis.text.y= element_text(size=22), 
        axis.title.x = element_blank(), axis.title.y  = element_blank(), 
        legend.title = element_text("Disturbance\nSeverity", size = 22), legend.position = c(0.8, 0.6), 
        legend.text = element_text(size = 22), panel.grid = element_blank(),
        axis.line = element_line(size = 0.5), plot.margin = margin(0.5,1,0,1, "cm")) +
  geom_point(size = 4) +
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50), limits = c(0,60),
                     sec.axis = sec_axis(~ .,labels = NULL, breaks = c(0, 10, 20, 30, 40, 50)))+
  geom_path(aes(group = factor(severity)), size = .9) +
  geom_errorbar(aes(ymin = NPP-SE, ymax = NPP+SE), width = 1.9, size = 0.9) +
  #coord_cartesian(clip = "off")+
  geom_vline(xintercept = as.numeric(as.Date("2019-05-21")), linetype = 2, 
             color = "firebrick3", size = 1) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%j",
               limits = as.Date(c("2019-04-02", "2019-11-11")))+
  annotate(geom="text", x = as.Date("2019-06-02", "%Y-%m-%d"), y = 19, label = c("*"),
           color="black", size = 12) +
  annotate(geom="text", x = as.Date("2019-06-09", "%Y-%m-%d"), y = 33, label = c("*"),
           color="black", size = 12) +
  annotate(geom="text", x = as.Date("2019-06-30", "%Y-%m-%d"), y = 43, label = c("*"),
           color="black", size = 12) +
  annotate(geom="text", x = as.Date("2019-07-07", "%Y-%m-%d"), y = 42, label = c("*"),
           color="black", size = 12) +
  annotate(geom="text", x = as.Date("2019-07-14", "%Y-%m-%d"), y = 51, label = c("*"),
           color="black", size = 12) +
  annotate(geom="text", x = as.Date("2019-07-21", "%Y-%m-%d"), y = 41.5, label = c("*"),
           color="black", size = 12) +
  annotate(geom="text", x = as.Date("2019-07-28", "%Y-%m-%d"), y = 40, label = c("*"),
           color="black", size = 12) +
  annotate(geom="text", x = as.Date("2019-08-04", "%Y-%m-%d"), y = 28, label = c("*"),
           color="black", size = 12) +
  annotate(geom="text", x = as.Date("2019-08-11", "%Y-%m-%d"), y = 22, label = c("*"),
           color="black", size = 12)+
  annotate(geom="text", as.Date("2019-04-7", "%Y-%m-%d"), y = 55, 
           label = c("A"), color="black", size = 7.5) 

##################################################################################
######## treatment time series 
treatment <- NPP_final %>% 
  group_by(week, treatment) %>% 
  summarize(NPP = mean(kgC_per_ha_day), SE = std.error(kgC_per_ha_day)) %>% 
  mutate(DOY = yday(week))

tx_time <- ggplot(treatment, aes(x = week, y = NPP, color = treatment)) +
  scale_colour_manual(values = c("#C4961A", "#293352"), name = "Disturbance\n  Type", labels = 
                        c("Bottom-Up", "Top-Down")) +
  geom_path(aes(group = treatment), size = .9) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = NPP-SE, ymax = NPP+SE), width = 1.9, size = 0.9) +
  #scale_linetype_manual(name = "Disturbance Type", 
                        # labels = c("Bottom Up", "Top Down")) +
  theme_bw() +
  scale_x_date(date_breaks = "20 days", date_labels = "%j",
               limits = as.Date(c("2019-04-02", "2019-11-11")))+
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50), limits = c(0,60), 
                     sec.axis = sec_axis(~ .,labels = NULL, breaks = c(0, 10, 20, 30, 40, 50)))+
  theme(axis.text.x= element_text(size = 22), axis.text.y= element_text(size=22), 
        axis.title.x = element_blank(), axis.title.y  = element_blank(), 
        legend.title = element_text(size = 22), legend.position = c(0.8, 0.6), 
        legend.text = element_text(size = 22), panel.grid = element_blank(),
        axis.line = element_line(size = 0.5), plot.margin = margin(-0.1,1,1.5,1, "cm")) +
  geom_vline(xintercept = as.numeric(as.Date("2019-05-21")), linetype = 2, 
             color = "firebrick3", size = 0.9) +
  annotate(geom="text", x = as.Date("2019-11-09", "%Y-%m-%d"), y = 60, 
           label = c("n.s."),
           color="black", size = 10) +
  annotation_custom(april, xmin = as.Date("2019-04-15", "%Y-%m-%d"),
                    xmax = as.Date("2019-04-15", "%Y-%m-%d"), ymin=-15,ymax=-15) +
  annotation_custom(may, xmin = as.Date("2019-05-15", "%Y-%m-%d"),
                    xmax = as.Date("2019-05-15", "%Y-%m-%d"), ymin=-15,ymax=-15) +
  annotation_custom(june, xmin = as.Date("2019-06-15", "%Y-%m-%d"),
                    xmax = as.Date("2019-06-15", "%Y-%m-%d"), ymin=-15,ymax=-15) +
  annotation_custom(july, xmin = as.Date("2019-07-15", "%Y-%m-%d"),
                    xmax = as.Date("2019-07-15", "%Y-%m-%d"), ymin=-15,ymax=-15) +
  annotation_custom(aug, xmin = as.Date("2019-08-15", "%Y-%m-%d"),
                    xmax = as.Date("2019-08-15", "%Y-%m-%d"), ymin=-15,ymax=-15) +
  annotation_custom(sept, xmin = as.Date("2019-09-15", "%Y-%m-%d"),
                    xmax = as.Date("2019-09-15", "%Y-%m-%d"), ymin=-15,ymax=-15) +
  annotation_custom(oct, xmin = as.Date("2019-10-15", "%Y-%m-%d"),
                    xmax = as.Date("2019-10-15", "%Y-%m-%d"), ymin=-15,ymax=-15) +
  annotation_custom(nov, xmin = as.Date("2019-11-15", "%Y-%m-%d"),
                    xmax = as.Date("2019-11-15", "%Y-%m-%d"), ymin=-15,ymax=-15) +
  coord_cartesian(clip = "off")+
  annotate(geom="text", as.Date("2019-04-7", "%Y-%m-%d"), y = 55, 
           label = c("B"), color="black", size = 7.5) 

# try other way to move axis labels around 
# g <- arrangeGrob(sev_time, tx_time, ncol=1, 
#                  left = grid::textGrob(expression(paste("Annual ANPP" [w], " ( ",kgC," ",ha^-1," ",yr^-1,")")), 
#                                        x = unit(0.5, "npc"), y = unit(0.5, "npc"), 
#                                        rot=90, gp=gpar(fontsize = 30)))
# grid.newpage()
# grid.draw(g)


####Paste these guys together! 
# g2 <- ggplotGrob(sev_time)
# g3 <- ggplotGrob(tx_time)
# g <- rbind(g2, g3, size = "first")

g <- arrangeGrob(sev_time, tx_time, ncol=1, 
            left = grid::textGrob(expression(paste("Daily ANPP" [w], " ( ",kgC," ",ha^-1," ",day^-1,")")), 
                                  x = unit(0.5, "npc"), y = unit(0.5, "npc"), 
                                  rot=90, gp=gpar(fontsize = 30)))
#g$widths <- unit.pmax(g2$widths, g3$widths)
grid.newpage()
grid.draw(g)
ggsave(filename = "figures/sev_tx_time.jpeg", dpi = 400, height = 8, width = 10, units = "in", g)
# ####

##################################################################################
######## Cumulative NPP histogram 
ggplot(annual_NPP_fig, aes(x = factor(severity), y = annual_NPP, fill = canopy_strata)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5) +
  theme_bw() +
  theme(axis.text.x= element_text(size = 22), axis.text.y= element_text(size=22), 
        axis.title.x = element_text(size=30), axis.title.y  = element_text(size = 30), 
        legend.title = element_blank(), legend.position = c(0.2, 0.9), 
        legend.text = element_text(size = 22), panel.grid = element_blank(),
        axis.line = element_line(size = 0.5), plot.margin = margin(1.5,1,1.5,1, "cm")) +
  scale_fill_manual(values = c("lightgreen", "darksalmon", "goldenrod1"), name = "", 
                    labels = c("Canopy", "Subcanopy", "Seedling/Sapling")) +
  scale_y_continuous(breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000), 
                     sec.axis = sec_axis(~ .,labels = NULL, breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000))) +
  labs(x = "Disturbance Severity (%)",
       y = expression(paste("Annual ANPP" [w], " ( ",kgC," ",ha^-1," ",yr^-1,")"))) +
  #theme(axis.title = element_text(size = 10)) +
  geom_errorbar(aes(ymin = errorbars-SE, ymax = errorbars+SE), width = 0.1) +
  annotate(geom="text", x = "0", y = 1700, label = c("a"),
           color="black", size = 7.5) +
  annotate(geom="text", x = "45", y = 1700, label = c("a"),
           color="black", size = 7.5) +
  annotate(geom="text", x = "65", y = 1700, label = c("a"),
           color="black", size = 7.5) +
  annotate(geom="text", x = "85", y = 1700, label = c("a"),
           color="black", size = 7.5) 

ggsave("figures/cumulative_NPP.jpeg", height = 8, width = 11, units = "in", 
       dpi = 400 ,plot = last_plot())

######### Fate by week time series 
# ggplot(fate_wk, aes(x = week , y = NPP, color = factor(fate)))+
#   geom_path(aes(group = factor(fate)), size = .9) 
  

####################################################################################
############### Live-Kill Stacked Histogram of Annual NPP


sev <- ggplot(NPP_fate_mean, aes(x = factor(severity), y = NPP_mean, fill = fate)) +
  geom_bar(stat = "identity", position = "stack", width = 0.4) +
  geom_line(aes(x= factor(severity), y = annual_NPP), group = 0) +
  geom_point(aes(x = factor(severity), y = annual_NPP), size = 3, show.legend = FALSE) +
  theme_bw() +
  theme(axis.text.x= element_text(size = 22), axis.text.y= element_text(size=22), 
        axis.title.x = element_text(size=30), axis.title.y  = element_blank(), 
        legend.title = element_blank(), legend.position = c(0.24, 0.87), 
        legend.text = element_text(size = 22), panel.grid = element_blank(),
        axis.line = element_line(size = 0.5), plot.margin = margin(1.5,1,.5,1, "cm")) +
  scale_fill_manual(values = c("orange3","palegreen3"), labels = c("Girdled", "Healthy"))+
  scale_y_continuous(breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000), 
                     sec.axis = sec_axis(~ .,labels = NULL, breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000))) +
  labs(x = "Disturbance Severity (%)",
       y = expression(paste("Annual ANPP" [w], " ( ",kgC," ",ha^-1," ",yr^-1,")"))) +
  geom_errorbar(aes(ymin = errorbars-SE, ymax = errorbars+SE), width = 0.1) +
  annotate(geom="text", x = 4.35, y = 3000, label = c("n.s."),
           color="black", size = 10) +
  annotate(geom="text", x = 0.6, y = 2700, 
           label = c("A"), color="black", size = 7.5) 

###################################################################################
############## Treatment stacked histogram 
tx <- ggplot(NPP_fate_tx) +
  geom_bar(aes(x = treatment, y = NPP_mean, fill = fate), stat = "identity", position = "stack", width = 0.4) +
  geom_bracket(xmin = "Bottom-Up",
               xmax = "Top-Down", y.position = 2900, label = "n.s.",
               label.size = 10, tip.length = c(0.09, .02), size = 0.6) +
  theme_bw() +
  theme(axis.text.x= element_text(size = 22), axis.text.y= element_text(size=22), 
        axis.title.x = element_text(size=30), axis.title.y  = element_blank(), 
        legend.title = element_blank(), legend.position = "none", 
        legend.text = element_text(size = 22), panel.grid = element_blank(),
        axis.line = element_line(size = 0.5), plot.margin = margin(0,1,2,1, "cm")) +
  scale_fill_manual(values = c("orange3","palegreen3"), labels = c("Girdled", "Ungirdled"))+
  scale_y_continuous(breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000),limits = c(0,3000),
                     sec.axis = sec_axis(~ .,labels = NULL, breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000))) +
  labs(x = "Disturbance Type",
       y = "") +
  geom_errorbar(aes(x = treatment, y = NPP_mean, ymin = errorbars-SE, ymax = errorbars+SE), 
                width = 0.1) +
  # annotate(geom="text", x = 2.45, y = 3000, label = c("n.s."),
  #          color="black", size = 10) +
  annotate(geom="text", x = 0.51, y = 2700,
           label = c("B"), color="black", size = 7.5) +  # annotate(geom="text", x=1.5, y=500, label = c("["), color = "black",
  #          size = 15, ) 
  annotate(geom="text", x=1.5, y=850, label = c("P<0.05"), color = "black",
           size = 10)

####Paste these guys together! 
# try other way to move axis labels around 
g <- arrangeGrob(sev, tx, ncol=1, 
                 left = grid::textGrob(expression(paste("Annual ANPP" [w], " ( ",kgC," ",ha^-1," ",yr^-1,")")), 
                                       x = unit(0.5, "npc"), y = unit(0.5, "npc"), 
                                       rot=90, gp=gpar(fontsize = 30)))
grid.newpage()
grid.draw(g)
ggsave(filename = "figures/fate_tx_sev.jpeg", dpi = 400, height = 11, width = 8, units = "in", g)


##################################################################################
############## PFT early-late successional side by side figure 
ggplot(annual_PFT, aes(x = factor(severity), y = NPP_canopy, fill = PFT)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.4) +
  theme_bw() +
  theme(axis.text.x= element_text(size = 22), axis.text.y= element_text(size=22), 
        axis.title.x = element_text(size=30), axis.title.y  = element_text(size = 30), 
        legend.title = element_blank(), legend.position = c(0.3, 0.9), 
        legend.text = element_text(size = 22), panel.grid = element_blank(),
        axis.line = element_line(size = 0.5), plot.margin = margin(1.5,1,1.5,1, "cm")) +
  scale_fill_manual(values = c("greenyellow","chartreuse4"), 
                    labels = c("Early Successional", "Middle/Late Succesional"))+
  scale_y_continuous(breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000), 
                     sec.axis = sec_axis(~ .,labels = NULL, breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000))) +
  labs(x = "Disturbance Severity (%)",
       y = expression(paste("Annual ANPP" [w], " ( ",kgC," ",ha^-1," ",yr^-1,")"))) + 
  geom_errorbar(position = "dodge", aes(ymin = NPP_canopy-SE, ymax = NPP_canopy+SE), width = 0.4) +
  annotate(geom="text", x = "65", y = 2130, label = c("*"),
           color="black", size = 12) +
  annotate(geom="text", x = "0", y = 1650, label = c("n.s."),
           color="black", size = 7.5) +
  annotate(geom="text", x = "45", y = 1690, label = c("n.s."),
           color="black", size = 7.5) +
  annotate(geom="text", x = "85", y = 1800, label = c("n.s."),
           color="black", size = 7.5) 

ggsave(filename = "figures/PFT_NPP.jpeg", dpi = 400, height = 8, width = 11, units = "in",plot = last_plot())

##################################################################################
########LAI FIGures by severity and type (treatment)

LAI_sev <- ggplot(LAI_severity, aes(x = factor(severity), y = LAI_mean, fill = factor(severity))) +
  geom_bar(stat = "identity", position = "dodge", width = 0.4) +
  theme_bw() +
  theme(axis.text.x= element_text(size = 22), axis.text.y= element_text(size=22), 
        axis.title.x = element_text(size=30), axis.title.y  = element_blank(), 
        legend.title = element_blank(), legend.position = c(0.5, 1.06), 
        legend.text = element_text(size = 22), legend.direction = "horizontal", 
        panel.grid = element_blank(), axis.line = element_line(size = 0.5), 
        plot.margin = margin(1.5,1,0.5,1, "cm")) +
  scale_y_continuous(limits = c(0,3.6), sec.axis = sec_axis(~ .,labels = NULL)) +
  # scale_x_discrete(position = "top") +
  # scale_x_discrete(position = "bottom")+
  scale_fill_manual(aesthetics = "fill", values = forte)+
  guides(fill = FALSE) +
  labs(x = "Disturbance Severity (%)",
       y = "LAI") +
  #theme(axis.title = element_text(size = 10)) +
  geom_errorbar(position = position_dodge(width = 0.4), 
                aes(ymin = LAI_mean-SE, ymax = LAI_mean+SE), width = 0.1, 
                color = "grey32", size = 1) +
  annotate(geom="text", x = 1, y = 3.5, label = c("a"), color="black", size = 10) +
  annotate(geom="text", x = 2, y = 3.6, label = c("a"), color="black", size = 10) +
  annotate(geom="text", x = 3, y = 3.15, label = c("ab"), color="black", size = 10) +
  annotate(geom="text", x = 4, y = 2.9, label = c("b"), color="black", size = 10) +
  annotate(geom="text", x = 0.6, y = 3.5, label = c("A"), color="black", size = 7.5) 

# and for disturbance type 

LAI_tx <- ggplot(LAI_type, aes(x = factor(treatment), y = LAI_mean, fill = factor(treatment))) +
  geom_bar(stat = "identity", position = "dodge", width = 0.4) +
  theme_bw() +
  theme(axis.text.x= element_text(size = 22), axis.text.y= element_text(size=22), 
        axis.title.x = element_text(size=30), axis.title.y  = element_blank(), 
        legend.title = element_blank(), legend.position = "blank", 
        legend.text = element_text(size = 22), panel.grid = element_blank(),
        axis.line = element_line(size = 0.5), plot.margin = margin(0,1,1.5,1, "cm")) +
  scale_y_continuous(limits = c(0,3.6), sec.axis = sec_axis(~ .,labels = NULL)) +
  scale_fill_manual(values = c("#C4961A", "#293352"), 
                    labels = c("DOY: 186-190", "DOY: 213-215"))+
  labs(x = "Disturbance Type",
       y = "LAI") +
  #theme(axis.title = element_text(size = 10)) +
  geom_errorbar(position = position_dodge(0.4), aes(ymin = LAI_mean-SE, ymax = LAI_mean+SE), 
                width = 0.1, color = "grey32", size = 1) +
  annotate(geom="text", x = 2.4, y = 3.5, label = c("n.s."),
           color="black", size = 10) +
  annotate(geom="text", x = 0.5, y = 3.5, label = c("B"), color="black", size = 7.5) 

# paste plots together 
# two ways but this first one is best. 


g <- arrangeGrob(LAI_sev, LAI_tx, ncol=1, 
                 left = grid::textGrob("LAI (dimensionless)", 
                                       x = unit(0.5, "npc"), y = unit(0.5, "npc"), 
                                       rot=90, gp=gpar(fontsize = 30)))
grid.newpage()
grid.draw(g)
ggsave(filename = "figures/LAI_tx_sev.jpeg", dpi = 400, height = 11, width = 9, units = "in", g)

##################################################################################
################# Treatment by severity df's and figures 

# treatement df for figures 
treatment <- NPP_final %>% 
  group_by(week, treatment) %>% 
  summarize(NPP = mean(kgC_per_ha_day), SE = std.error(kgC_per_ha_day)) %>% 
  mutate(DOY = yday(week))

## Treatment-severity time series 
ggplot(treatment, aes(x = week, y = NPP, linetype = factor(treatment))) +
  geom_point(size = 0.9) +
  geom_line(size = 0.9) +
  geom_errorbar(aes(ymin = NPP-SE, ymax = NPP+SE), width = 0.9, size = 0.9, color = "#009E73") +
  scale_linetype_manual(name = "Disturbance Type", 
                        labels = c("Bottom Up", "Top Down")) +
  labs(x = "" , 
       y = expression(paste("ANPP" [w], " ( ",kgC," ",ha^-1," ",day^-1,")"))) +
  theme_bw() +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%j",
               limits = as.Date(c("2019-04-07", "2019-11-11")))+
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50), limits = c(0,50), 
                     sec.axis = sec_axis(~ .,labels = NULL, breaks = c(0, 10, 20, 30, 40, 50)))+
  theme(axis.text.x= element_blank(), axis.text.y= element_text(size=17), 
        axis.title.x = element_blank(), axis.title.y  = element_blank(), 
        legend.title = element_text("Disturbance\nSeverity", size = 17), legend.position = c(0.8, 0.6), 
        legend.text = element_text(size = 15), panel.grid = element_blank(),
        axis.line = element_line(size = 0.5), plot.margin = margin(.5,1,1,1, "cm"),
        plot.title = element_text(size = 20)) +
  geom_vline(xintercept = as.numeric(as.Date("2019-05-21")), linetype = 2, 
             color = "firebrick3", size = 0.9) +
  theme(plot.title = element_text(hjust=0.5)) +
  annotate(geom="text", x = as.Date("2019-04-20", "%Y-%m-%d"), y = 45, 
           label = c("n.s."),
           color="black", size = 7) 

# treatement df for figures 
treatment65 <- NPP_final%>% 
  filter(severity == 0.65) %>% 
  group_by(week, severity, treatment) %>% 
  summarize(NPP = mean(kgC_per_ha_day), SE = std.error(kgC_per_ha_day)) %>% 
  mutate(DOY = yday(week))

## Treatment-severity time series 
tx_65 <- ggplot(treatment65, aes(x = week, y = NPP, linetype = treatment)) +
  geom_point(color = "#0072B2") +
  geom_line(color = "#0072B2", size = 0.9) +
  geom_errorbar(aes(ymin = NPP-SE, ymax = NPP+SE), width = 0.9, size = 0.9, color = "#0072B2") +
  scale_linetype_manual(values = c(1, 6), name = "Treatment", 
                        labels = c("Bottom Up", "Top Down")) +
  labs(x = "" , 
       y = expression(paste("ANPP" [w], " ( ",kgC," ",ha^-1," ",day^-1,")"))) +
  theme_bw() +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%j",
               limits = as.Date(c("2019-04-07", "2019-11-11")))+
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70), limits = c(0,70),
                     sec.axis = sec_axis(~ .,labels = NULL, breaks = c(0, 10, 20, 30, 40, 50, 60, 70)))+
  theme(axis.text.x= element_blank(), axis.text.y= element_text(size=17), 
        axis.title.x = element_blank(), axis.title.y  = element_blank(), 
        legend.title = element_blank(), legend.position = "none", 
        legend.text = element_text(size = 15), panel.grid = element_blank(),
        axis.line = element_line(size = 0.5), plot.margin = margin(.5,1,0,1, "cm"),
        plot.title = element_text(size = 20)) +
  geom_vline(xintercept = as.numeric(as.Date("2019-05-21")), linetype = 2, 
             color = "firebrick3") +
  ggtitle("Severity - 65%") +
  theme(plot.title = element_text(hjust=0.5)) +
  annotate(geom="text", x = as.Date("2019-04-20", "%Y-%m-%d"), y = 45, 
           label = c("n.s"),
           color="black", size = 7) 

# treatement df for figures 
treatment85 <- NPP_final %>% 
  filter(severity == 0.85) %>% 
  group_by(week, severity, treatment) %>% 
  summarize(NPP = mean(kgC_per_ha_day), SE = std.error(kgC_per_ha_day)) %>% 
  mutate(DOY = yday(week))

## Treatment-severity time series 
tx_85 <- ggplot(treatment85, aes(x = week, y = NPP, linetype = treatment)) +
  geom_point(color = "#D55E00") +
  geom_line(color = "#D55E00", size = 0.9) +
  geom_errorbar(aes(ymin = NPP-SE, ymax = NPP+SE), width = 0.9, size = 0.9, color = "#D55E00") +
  scale_linetype_manual(values = c(1, 6), name = "Treatment", 
                        labels = c("Bottom Up", "Top Down")) +
  labs(x = "" , 
       y = expression(paste("ANPP" [w], " ( ",kgC," ",ha^-1," ",day^-1,")"))) +
  theme_bw() +
  theme(axis.text.x= element_text(size = 17), axis.text.y= element_text(size=17), 
        axis.title.x = element_blank(), axis.title.y  = element_text(size = 17), 
        legend.title = element_blank(), legend.position = "none", 
        legend.text = element_text(size = 15), panel.grid = element_blank(),
        axis.line = element_line(size = 0.5), plot.margin = margin(.5,1,1,1, "cm"), 
        plot.title = element_text(size = 20)) +
  geom_vline(xintercept = as.numeric(as.Date("2019-05-21")), linetype = 2, 
             color = "firebrick3") +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%j",
               limits = as.Date(c("2019-04-07", "2019-11-11")))+
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70), 
                     sec.axis = sec_axis(~ .,labels = NULL, breaks = c(0, 10, 20, 30, 40, 50, 60, 70)))+
  annotation_custom(april, xmin = as.Date("2019-04-15", "%Y-%m-%d"),
                    xmax = as.Date("2019-04-15", "%Y-%m-%d"), ymin=-14,ymax=-14) + 
  annotation_custom(may, xmin = as.Date("2019-05-15", "%Y-%m-%d"),
                    xmax = as.Date("2019-05-15", "%Y-%m-%d"), ymin=-14,ymax=-14) + 
  annotation_custom(june, xmin = as.Date("2019-06-15", "%Y-%m-%d"),
                    xmax = as.Date("2019-06-15", "%Y-%m-%d"), ymin=-14,ymax=-14) + 
  annotation_custom(july, xmin = as.Date("2019-07-15", "%Y-%m-%d"),
                    xmax = as.Date("2019-07-15", "%Y-%m-%d"), ymin=-14,ymax=-14) +
  annotation_custom(aug, xmin = as.Date("2019-08-15", "%Y-%m-%d"),
                    xmax = as.Date("2019-08-15", "%Y-%m-%d"), ymin=-14,ymax=-14) + 
  annotation_custom(sept, xmin = as.Date("2019-09-15", "%Y-%m-%d"),
                    xmax = as.Date("2019-09-15", "%Y-%m-%d"), ymin=-14,ymax=-14) +
  annotation_custom(oct, xmin = as.Date("2019-10-15", "%Y-%m-%d"),
                    xmax = as.Date("2019-10-15", "%Y-%m-%d"), ymin=-14,ymax=-14) + 
  annotation_custom(nov, xmin = as.Date("2019-11-15", "%Y-%m-%d"),
                    xmax = as.Date("2019-11-15", "%Y-%m-%d"), ymin=-14,ymax=-14) +
  coord_cartesian(clip = "off")+
  ggtitle("Severity - 85%") +
  theme(plot.title = element_text(hjust=0.5)) +
  annotate(geom="text", x = as.Date("2019-04-20", "%Y-%m-%d"), y = 45, 
           label = c("n.s"),
           color="black", size = 7) 

g2 <- ggplotGrob(tx_45)
g3 <- ggplotGrob(tx_65)
g4 <- ggplotGrob(tx_85)
g <- rbind(g2, g3, g4, size = "first")
g$widths <- unit.pmax(g2$widths, g3$widths, g4$widths)
grid.newpage()
grid.draw(g)

ggsave("figures/tx_timeseries.jpeg", height = 11, width = 9, units = "in", g)


# paste all 3 tx figures onto 1 
library(gridExtra)
timeseries <- grid.arrange(plot_tx45, plot_tx65, plot_tx85, ncol=1)
g <- arrangeGrob(plot_tx45, plot_tx65, plot_tx85, ncol=1)
ggsave("figures/treatment.jpeg", height = 11, width = 9, units = "in", g)

##############################################################################
# live-dead histrogram for UMBS winter meeting 

fate_hist <- RGR_fate %>% 
  filter(severity != "0.00") %>% 
  group_by(severity, fate) %>% 
  summarize(RGR_mean = mean(RGR)*10, SE = std.error(RGR)) #mm/cm

# Cumulative NPP split by live-dead bar graph for winter meeting 

ggplot(fate_hist, aes(x = factor(severity), y = RGR_mean, fill = fate)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.4) +
  theme_classic() +
  scale_fill_manual(values = c("firebrick2","darkgreen"), name = "Fate")+
  labs(x = "Disturbance Severity",
       y = expression(paste( "Mean RGR ( ",mm," ",cm^-1," ",day^-1,")"))) +
  theme(axis.title = element_text(size = 10))+
  annotate(geom="text", x = "0.45", y = 0.00072, label = c("n.s."),
           color="black", size = 5) +
  annotate(geom="text", x = "0.65", y = 0.00115, label = c("*"),
           color="black", size = 8) +
  annotate(geom="text", x = "0.85", y = 0.00085, label = c("n.s."),
           color="black", size = 5) +
  geom_errorbar(position = "dodge", aes(ymin = RGR_mean-SE, ymax = RGR_mean+SE), width = 0.4) 

ggsave(filename = "figures/RGR_histogram.jpeg", plot = last_plot())

# treatment histogram for UMBS winter meeting 

treatment_hist <- annual_inc %>% 
  filter(severity != "0%") %>% 
  group_by(severity, treatment) %>% 
  summarise(NPP = mean(NPP_canopy), SE = std.error(NPP_canopy))

ggplot(treatment_hist, aes(x = factor(severity), y = NPP, fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.4) +
  theme_classic() +
  scale_fill_manual(values = c("coral4","darkolivegreen4"), name = "Disturbance\nType")+
  labs(x = "Disturbance Severity",
       y = expression(paste("NPP (",kgC," ",ha^-1," ",day^-1,")"))) +
  theme(axis.title = element_text(size = 10))+
  #annotate(geom="text", x = "65%", y = 70, label = c("*"),
           #color="black", size = 8) +
  #annotate(geom="text", x = "85%", y = 65, label = c("*"),
           #color="black", size = 8) +
  #annotate(geom="text", x = "45%", y = 3000, label = c("n.s."),
           #color="black", size = 5) +
  geom_errorbar(position = "dodge", aes(ymin = NPP-SE, ymax = NPP+SE), width = 0.4) 

ggsave(filename = "figures/treatment_histogram.jpeg", plot = last_plot())


###################################################################################
#### ALL time subcanopy and canopy NPP by year and severity 

# make a df with the text I want to annotate 
ann_text <- data.frame(severity = c("0","45", "65", "85", "0","45","65", "85"), 
                       lab = c("ab", "ab", "ab", "ab", "a", "bc", "bc", "c"), 
                       NPP1 = c(600, 600, 600, 600, 600, 600, 950, 750), 
                       year = factor(c(2019, 2019, 2019, 2019, 2020, 2020, 2020, 2020),
                                     levels = c(2019,2020)))

# figure for Chris, subcanopy by severity faceted by year
sc <- ggplot(NPP_sc_alltime, aes(x = factor(severity), y = NPP1, fill = factor(severity))) +
  geom_boxplot()+
  theme_classic()+
  facet_grid(.~factor(year), scales = "free")+
  scale_fill_manual(values = forte)+
  theme(axis.text.x= element_text(size = 20), axis.text.y= element_text(size=20), 
        axis.title.x = element_text(size = 25), axis.title.y  = element_text(size=25), 
        legend.title = element_blank(),  strip.text.x =element_text(size = 25), 
        legend.text = element_blank(), legend.position = "none", 
        plot.margin = margin(1,1,1,1, "cm")) +
  scale_y_continuous(sec.axis = sec_axis(~ .,labels = NULL)) +
  labs(x = "Disturbance Severity (%)",
       y = expression(paste("ANPP" [w], " ( ",kgC," ",ha^-1," ",year^-1,")")))+
  geom_text(data = ann_text, label = c("ab", "ab", "ab", "ab", "a", "bc", "bc", "c"), size = 10)

# figure for Chris, canopy by severity faceted by year
c <- ggplot(NPP_alltime, aes(x = factor(severity), y = kgC_ha_yr, fill = factor(severity))) +
  geom_boxplot()+
  theme_classic()+
  facet_grid(.~factor(year), scales = "free")+
  scale_fill_manual(values = forte)+
  theme(axis.text.x= element_text(size = 20), axis.text.y= element_text(size=20), 
        axis.title.x = element_text(size = 25), axis.title.y  = element_text(size=25), 
        legend.title = element_blank(),  strip.text.x =element_text(size = 25), 
        legend.text = element_blank(), legend.position = "none", 
        plot.margin = margin(1,1,1,1, "cm")) +
  scale_y_continuous(sec.axis = sec_axis(~ .,labels = NULL)) +
  labs(x = "Disturbance Severity (%)",
       y = expression(paste("ANPP" [w], " ( ",kgC," ",ha^-1," ",year^-1,")")))+
  annotate(geom="text", x = 1.5, y = 5500, label = c("n.s."),
           color="black", size = 10) 

# paste these two together 
g <- ggarrange(c, sc, labels = c("A. Canopy", "B. Subcanopy"), 
               font.label = list(size = 25), vjust = 1)

ggsave("figures/sc_c_yr_sev.jpeg", plot = g, units = "in", height = 5, width = 11)

# take sc_alltime and vector math for log(disturbance/control)
sc_log <- NPP_sc_alltime %>% 
  mutate(replicate = substr(subplot, 1, 1)) %>%
  group_by(replicate, year, severity) %>% 
  summarize(NPP_sev = mean(NPP1)) %>% 
  group_by(year, replicate) %>% 
  mutate(log_ratio = case_when(
    year == "2019" ~ log(NPP_sev/NPP_sev[1]),
    year == "2020" ~ log(NPP_sev/NPP_sev[1])
  )) %>% 
  group_by(severity, year) %>% 
  summarize(avg_log_ratio = mean(log_ratio), SE = std.error(log_ratio))
     
# do the same for the canopy 
canopy_log <- NPP_alltime %>% 
  mutate(replicate = substr(subplot, 1, 1)) %>%
  group_by(replicate, year, severity) %>% 
  summarize(NPP_sev = mean(kgC_ha_yr)) %>% 
  group_by(year, replicate) %>% 
  mutate(log_ratio = case_when(
    year == "2019" ~ log(NPP_sev/NPP_sev[1]),
    year == "2020" ~ log(NPP_sev/NPP_sev[1])
  )) %>% 
  group_by(severity, year) %>% 
  summarize(avg_log_ratio = mean(log_ratio), SE = std.error(log_ratio))

#make year a number
sc_log$year <- as.numeric(sc_log$year)
canopy_log$year <- as.numeric(canopy_log$year)

# create time series of log response ratio for the SUBCANOPY
sc_resi <- ggplot(sc_log, aes(x = year, y = avg_log_ratio, color = severity)) +
  scale_color_manual(values = forte)+
  theme_classic()+
  geom_point(size = 3)+
  geom_path(size = 1)+
  scale_x_continuous(breaks = c(2019,2020), sec.axis = sec_axis(~ .,labels = NULL,breaks = c(2018,2019,2020))) +
  theme(axis.text= element_text(size = 20), axis.title = element_text(size = 25), 
        legend.text = element_text(size = 15), legend.title = element_text(size = 15), 
        legend.position = c(0.2, 0.8), plot.margin = margin(1,1,1,1, "cm")) +
  scale_y_continuous(sec.axis = sec_axis(~ .,labels = NULL)) +
  geom_errorbar(aes(ymin=avg_log_ratio - SE, ymax=avg_log_ratio + SE), width = 0.05)+
  labs(y = "Log(ANPPw Disturbance/Control)")

# create time series of log response ratio for the CANOPY
c_resi <- ggplot(canopy_log, aes(x = year, y = avg_log_ratio, color = severity)) +
  scale_color_manual(values = forte)+
  theme_classic()+
  geom_point(size = 3)+
  geom_path(size = 1)+
  scale_x_continuous(breaks = c(2019,2020), sec.axis = sec_axis(~ .,labels = NULL,breaks = c(2018,2019,2020))) +
  theme(axis.text= element_text(size = 20), axis.title = element_text(size = 25), 
        legend.text = element_text(size = 15), legend.title = element_text(size = 15), 
        legend.position = c(0.2, 0.8), plot.margin = margin(1,1,1,1, "cm")) +
  scale_y_continuous(sec.axis = sec_axis(~ .,labels = NULL)) +
  geom_errorbar(aes(ymin=avg_log_ratio - SE, ymax=avg_log_ratio + SE), width = 0.05)+
  labs(y = "Log(ANPPw Disturbance/Control)")

# paste these two together 
g <- ggarrange(sc_resi, c_resi, labels = c("C", "D"), font.label = list(size = 25))



###################################################################################
#### ALL time subcanopy and canopy NPP by year and type

# plot subcanopy npp by type and faceted by year
sc_type <- ggplot(NPP_sc_alltime, aes(x = factor(treatment), y = NPP1, fill = factor(treatment))) +
  geom_boxplot()+
  theme_classic()+
  facet_grid(.~factor(year), scales = "free")+
  scale_fill_manual(values = c("#C4961A", "#293352"))+
  theme(axis.text.x= element_text(size = 20), axis.text.y= element_text(size=20), 
        axis.title.x = element_text(size = 25), axis.title.y  = element_text(size=25), 
        legend.title = element_blank(),  strip.text.x =element_text(size = 25), 
        legend.text = element_blank(), legend.position = "none", 
        plot.margin = margin(1,1,1,1, "cm")) +
  scale_y_continuous(sec.axis = sec_axis(~ .,labels = NULL)) +
  labs(x = "Disturbance Type",
       y = expression(paste("ANPP" [w], " ( ",kgC," ",ha^-1," ",year^-1,")")))+
  annotate(geom="text", x = 1.5, y = 1000, label = c("n.s."),
           color="black", size = 10) 

# plot canopy npp by type and faceted by year
c_type <- ggplot(NPP_alltime, aes(x = factor(treatment), y = kgC_ha_yr, fill = factor(treatment))) +
  geom_boxplot()+
  theme_classic()+
  facet_grid(.~factor(year), scales = "free")+
  scale_fill_manual(values = c("#C4961A", "#293352"))+
  theme(axis.text.x= element_text(size = 20), axis.text.y= element_text(size=20), 
        axis.title.x = element_text(size = 25), axis.title.y  = element_text(size=25), 
        legend.title = element_blank(),  strip.text.x =element_text(size = 25), 
        legend.text = element_blank(), legend.position = "none", 
        plot.margin = margin(1,1,1,1, "cm")) +
  scale_y_continuous(sec.axis = sec_axis(~ .,labels = NULL)) +
  labs(x = "Disturbance Type",
       y = expression(paste("ANPP" [w], " ( ",kgC," ",ha^-1," ",year^-1,")")))+
  annotate(geom="text", x = 1.5, y = 5000, label = c("n.s."),
           color="black", size = 10) 

# paste these two together 
g <- ggarrange(c_type, sc_type, labels = c("A. Canopy", "B. Subcanopy"), 
               font.label = list(size = 25), vjust = 1)

ggsave("figures/sc_c_yr_type.jpeg", plot = g, units = "in", height = 5, width = 11)
