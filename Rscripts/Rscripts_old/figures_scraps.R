#############################################################################
############data visualization: RGR by DBH_cm of all trees with sig DBH effect

###################################################################################
######## Fate-severity box plot
# 
# #annual NPP boxplot including all canopy strata 
# ggplot(annual_inc1, aes(x = factor(severity), y = annual_NPP)) +
#   geom_boxplot()
# 
# # annual NPP boxplot for healthy 
# # filter NPP fate for healthy trees
# NPP_healthy <- NPP_fate %>% 
#   filter(fate == "live")
# 
# ggplot(NPP_healthy, aes(x = factor(severity), y = kgC_ha_yr)) +
#   geom_boxplot()
# 
# # annual NPP boxplot for senescent
# # filter NPP_fate for kill trees
# NPP_senescent <- NPP_fate %>% 
#   filter(fate == "kill")
# 
# ggplot(NPP_senescent, aes(x = factor(severity), y = kgC_ha_yr)) +
#   geom_boxplot()


##############################################################################
####### Exploring replicates 
# create 4 panel figures with each replicate broken out by severity

# df with treatments and replicate to have a look at replicate comparisons
# rep <- NPP_total %>% 
#   group_by(week, severity) %>% 
#   summarize(NPP = mean(NPP_canopy), SE = std.error(NPP_canopy))
# 
# # filter for each rep and plot; then paste plots together using plottrix 
# rep_A <- rep %>% 
#   filter(replicate == "A")
# 
# plot_tx45 <- ggplot(rep_A, aes(x = week, y = NPP, linetype = treatment)) +
#   geom_point(color = "#009E73") +
#   geom_line(color = "#009E73") +
#   geom_errorbar(aes(ymin = NPP-SE, ymax = NPP+SE), width = 0.8, color = "#009E73") +
#   scale_linetype_manual(values = c(1, 6), name = "Treatment", 
#                         labels = c("Bottom Up", "Top Down")) +
#   labs(x = "" , 
#        y = expression(paste("NPP (",kgC," ",ha^-1," ",day^-1,")"))) +
#   ylim(0, 170) +
#   theme_classic() +
#   theme(legend.position = c(0.8, 0.6)) +
#   geom_vline(xintercept = as.numeric(as.Date("2019-05-21")), linetype = 2, 
#              color = "firebrick3") +
#   ggtitle("45% Defoliation") +
#   theme(plot.title = element_text(hjust=0.5)) +
#   annotate(geom="text", x = as.Date("2019-04-20", "%Y-%m-%d"), y = 130, 
#            label = c("N.S."),
#            color="black", size = 3.5) 



# This way works just fine but can't figure out how to manually put in axis labels 
# g2 <- ggplotGrob(sev)
# g3 <- ggplotGrob(tx)
# g <- rbind(g2, g3, size = "first")
# g$widths <- unit.pmax(g2$widths, g3$widths)
# grid.newpage()
# grid.draw(g)
####