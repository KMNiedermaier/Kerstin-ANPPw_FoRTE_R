

#############################################################################
####### Above plots with alpha shading and early AND late season LAI 
### include early season dates 
# LAI_end <- LAI %>% 
#   filter(DOY == 214 | DOY == 213 | DOY == 215 | DOY == 189 | DOY == 186 | 
#            DOY == 190) %>% 
#   arrange(subplot) %>% 
#   mutate(season = case_when(
#     DOY == 214 | DOY == 213 | DOY == 215 ~ "early",
#     DOY == 186 | DOY == 189 | DOY == 190 ~ "late"
#   ))



# ggplot(LAI_severity, aes(x = factor(severity), y = LAI_mean, fill = factor(severity), alpha = season)) +
#   geom_bar(stat = "identity", 
#            position = "dodge", width = 0.4) +
#   scale_alpha_discrete(range = c(0.5,1), 
#                        labels = c("DOY: 186-190", "DOY: 213-215"))+
#   theme_bw() +
#   theme(axis.text.x= element_text(size = 17), axis.text.y= element_text(size=17), 
#         axis.title.x = element_text(size=19), axis.title.y  = element_text(size = 19), 
#         legend.title = element_blank(), legend.position = c(0.5, 1.06), 
#         legend.text = element_text(size = 17), legend.direction = "horizontal", 
#         panel.grid = element_blank(), axis.line = element_line(size = 0.5), 
#         plot.margin = margin(1.5,1,0.5,1, "cm")) +
#   scale_y_continuous(limits = c(0,3.5)) +
#   scale_fill_manual(aesthetics = "fill", values = forte,
#                     labels = c("DOY: 186-190", "DOY: 213-215"))+
#   guides(fill = FALSE) +
#   labs(x = "Disturbance Severity (%)",
#        y = "LAI") +
#   theme(axis.title = element_text(size = 10)) +
#   geom_errorbar(position = position_dodge(width = 0.4), aes(ymin = LAI_mean-SE, 
#                                                             ymax = LAI_mean+SE), width = 0.2) +
#   annotate(geom="text", x = 3.9, y = 2.9, label = c("*"),
#            color="black", size = 12) +
#   annotate(geom="text", x = 4.1, y = 3.2, label = c("*"),
#            color="black", size = 12) 

# TREATMENTS
# ggplot(LAI_type, aes(x = factor(treatment), y = LAI_mean, alpha = season, fill = factor(treatment))) +
#   geom_bar(stat = "identity", position = "dodge", width = 0.4) +
#   theme_bw() +
#   scale_alpha_discrete(range = c(0.5,1), 
#                        labels = c("DOY: 186-190", "DOY: 213-215"))+
#   theme(axis.text.x= element_text(size = 17), axis.text.y= element_text(size=17), 
#         axis.title.x = element_text(size=19), axis.title.y  = element_text(size = 19), 
#         legend.title = element_blank(), legend.position = "blank", 
#         legend.text = element_text(size = 17), panel.grid = element_blank(),
#         axis.line = element_line(size = 0.5), plot.margin = margin(0,1,1.5,1, "cm")) +
#   scale_y_continuous(limits = c(0,3.5)) +
#   scale_fill_manual(values = c("olivedrab4", "firebrick4"), 
#                     labels = c("DOY: 186-190", "DOY: 213-215"))+
#   labs(x = "Disturbance Type",
#        y = "LAI") +
#   theme(axis.title = element_text(size = 10)) +
#   geom_errorbar(position = position_dodge(0.4), aes(ymin = LAI_mean-SE, ymax = LAI_mean+SE), 
#                 width = 0.1) +
#   annotate(geom="text", x = 2.5, y = 3.5, label = c("n.s."),
#            color="black", size = 7) 