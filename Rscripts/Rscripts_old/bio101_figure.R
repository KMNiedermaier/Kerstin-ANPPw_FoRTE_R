live <- RGR_plot %>% 
  filter(subplot == "D02E", species == "QURU", fate == "live") %>% 
  summarise(sum(increment_kg)) 

kill <- RGR_plot %>% 
  filter(subplot == "D02E", species == "QURU", fate == "kill") %>% 
  summarise(sum(increment_kg))

single_plot_growth <- full_join(live, kill)
names(single_plot_growth) <- c("growth_kg")

single_plot_growth <- single_plot_growth %>% 
  mutate(tree= c("Oak_live", "Oak_dead"))

single_plot_growth <- single_plot_growth[c("tree", "growth_kg")]

ggplot(single_plot_growth, aes(x=tree, y=growth_kg, fill= tree)) +
  geom_col() +
  scale_fill_manual(values = c("brown4", "forestgreen")) +
  ylab("Tree Growth (kg of biomass)") +
  xlab("") +
  theme_classic()
  
  theme(legend.position= "none") 
  
  
  
