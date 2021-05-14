#####################################################################################
##########Species specific annual RGR figure 

# first take df, group by tag and sum weekly RGR_obs_est to get annual RGR
RGR_fig_2020 <- df %>% 
  filter(week == 17 | week == 18) %>% 
  group_by(tag, species, fate) %>% 
  mutate(annual_inc = DBH_cm_new - lag(DBH_cm_new, default = first(DBH_cm_new))) %>%
  mutate(annual_RGR = annual_inc / lag(DBH_cm_new, default = first(DBH_cm_new))) %>% 
  filter(week == 18) %>% 
  group_by(species, fate) %>% 
  summarise(RGR_2020 = mean(annual_RGR), SE = std.error(annual_RGR)) %>% 
  filter(species == "ACRU" | species == "ACSA" | species == "FAGR" | 
           species == "POGR" | species =="PIST" | species == "QURU") %>% 
  ungroup()

# recode species codes as species names 
RGR_fig_2020$species <- recode(RGR_fig_2020$species, "ACRU" = "Acer rubrum", "ACSA" = "Acer saccharum", 
                          "FAGR" = "Fagus grandifolia", "PIST" = "Pinus strobus", 
                          "POGR" = "Populus grandidentata", "QURU" = "Quercus rubra") 

# create x-axis label with species names 
xlabs <- c("Acer\nrubrum", "Acer\nsaccharum", "Fagus\ngrandifolia", "Pinus\nstrobus",
           "Populus\ngrandidentata" , "Quercus\nrubra")

# plot a figure with annual RGR by species/fate
ggplot(RGR_fig_2020, aes(x = species, y = RGR_2020, fill = fate)) +
  geom_bar(stat = "identity", position = "dodge")+
  theme_bw() +
  theme(axis.text.x= element_text(size = 18, face = "italic"), 
        axis.text.y= element_text(size=22), 
        axis.title.x = element_blank(), axis.title.y  = element_text(size=30), 
        legend.title = element_blank(), legend.position = c(0.85, 0.85), 
        legend.text = element_text(size = 22), legend.direction = "vertical", 
        panel.grid = element_blank(), axis.line = element_line(size = 0.5), 
        plot.margin = margin(1.5,1,0.5,1, "cm")) +
  scale_y_continuous(breaks = c(0, 0.01, 0.02, 0.03, 0.04), limits = c(0,0.04))+ 
  labs(y = expression(paste("RGR Nov-July", " ( ",cm," ",cm^-1," ",yr^-1,")")))+
  scale_x_discrete(labels = xlabs) +
  scale_fill_manual(values = c("orange3","palegreen3"), labels = c("Girdled", "Healthy"))+
  geom_errorbar(position = position_dodge2(0.7, padding = 0.75), 
                aes(ymin = RGR_2020-SE, ymax = RGR_2020+SE)) +
  annotate(geom="text", x = "Acer rubrum", y = 0.014, label = c("*"),
           color="black", size = 12)+
  annotate(geom="text", x = "Acer saccharum", y = 0.015, label = c("*"),
           color="black", size = 12) +
  annotate(geom="text", x = "Fagus grandifolia", y = 0.015, label = c("*"),
             color="black", size = 12)+
  # annotate(geom="text", x = "Pinus strobus", y = 0.021, label = c("n.s."),
  #          color="black", size = 12)+
  annotate(geom="text", x = "Populus grandidentata", y = 0.012, label = c("*"),
           color="black", size = 12) +
  annotate(geom="text", x = "Quercus rubra", y = 0.014, label = c("*"),
           color="black", size = 12)

# save to figures folder 
ggsave(filename = "figures/RGR_sp_2020.jpeg", dpi = 400, height = 8, width = 11, 
       units = "in", last_plot())

# create a df to run stats on all the species/fate RGRs
RGR_stats_2020 <- df %>% 
  filter(week == 17 | week == 18) %>% 
  group_by(tag, species, fate) %>% 
  mutate(annual_inc = DBH_cm_new - lag(DBH_cm_new, default = first(DBH_cm_new))) %>%
  mutate(annual_RGR = annual_inc / lag(DBH_cm_new, default = first(DBH_cm_new))) %>% 
  filter(week == 18) %>% 
  filter(species == "ACRU" | species == "ACSA" | species == "FAGR" | 
           species == "POGR" | species =="PIST" | species == "QURU") %>% 
  ungroup()

# create a t.test function that filters RGR_stats for a given species and runs a 
# paired t.test on RGR by fate 
t_test_RGR <- function(DF, SPECIES){
  filtered <- filter(DF, species == SPECIES)
  output <- t.test(annual_RGR ~ fate, data = filtered, paired = FALSE)
  return(output)
}
# run a paired t-test on RGR of girdled and ungirdled trees of each species 
t_test_RGR(RGR_stats_2020,"ACRU")
t_test_RGR(RGR_stats_2020,"ACSA")
t_test_RGR(RGR_stats_2020,"FAGR")
t_test_RGR(RGR_stats_2020,"PIST")
t_test_RGR(RGR_stats_2020,"POGR")
t_test_RGR(RGR_stats_2020,"QURU")
