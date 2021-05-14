# ## ACRU; MUST RUN TOGETHER 
# t_test_acru <- DBLband_transform %>% 
#   filter(species == "ACRU" & fate == "live" & week != "2019-07-14" & 
#            week != "2019-07-21") 
# t.test(increment ~ bot_or_top, data = t_test_acru, paired = TRUE)
# t_test_acru %>% 
#   ungroup() %>% 
#   filter(bot_or_top == "bot_increment") %>% 
#   summarize(mean = mean(increment))
# 
# t_test_acru <- DBLband_transform %>% 
#   filter(species == "ACRU" & fate == "kill" & week != "2019-07-14"& week != "2019-07-21") 
# t.test(increment ~ bot_or_top, data = t_test_acru, paired = TRUE)  
# t_test_acru %>% 
#   ungroup() %>% 
#   filter(bot_or_top == "top_increment") %>% 
#   summarize(mean = mean(increment))
# 
# ## ACSA; MUST RUN TOGETHER 
# t_test_acsa <- DBLband_transform %>% 
#   filter(species == "ACSA" & fate == "live" & week != "2019-07-07" & week != "2019-07-14" & week != "2019-11-10") 
# t.test(increment ~ bot_or_top, data = t_test_acsa, paired = TRUE)
# 
# t_test_acsa <- DBLband_transform %>% 
#   filter(species == "ACSA" & fate == "kill" & week != "2019-07-07" & week != "2019-07-14" & week != "2019-11-10") 
# t.test(increment ~ bot_or_top, data = t_test_acsa, paired = TRUE)
# 
# ## PIST: MUST RUN TOGETHER 
# t_test_pist <- DBLband_transform %>% 
#   filter(species == "PIST" & fate == "live") 
# t.test(increment ~ bot_or_top, data = t_test_pist, paired = TRUE)
# 
# t_test_pist <- DBLband_transform %>% 
#   filter(species == "PIST" & fate == "kill") 
# t.test(increment ~ bot_or_top, data = t_test_pist, paired = TRUE)
# 
# ## POGR: MUST RUN TOGETHER 
# t_test_pogr <- DBLband_transform %>% 
#   filter(species == "POGR" & fate == "live") 
# t.test(increment ~ bot_or_top, data = t_test_pogr, paired = TRUE)
# 
# t_test_pogr <- DBLband_transform %>% 
#   filter(species == "POGR" & fate == "kill") 
# t.test(increment ~ bot_or_top, data = t_test_pogr, paired = TRUE)
# 
# ## QURU: MUST RUN TOGETHER 
# t_test_quru <- DBLband_transform %>% 
#   filter(species == "QURU" & fate == "live") 
# t.test(increment ~ bot_or_top, data = t_test_quru, paired = TRUE)
# 
# t_test_quru <- DBLband_transform %>% 
#   filter(species == "QURU" & fate == "kill") 
# t.test(increment ~ bot_or_top, data = t_test_quru, paired = TRUE)




################ not really sure what I was going for down here....############################
# Plots top against bottom increment for both live and kill trees
# filter out start dates 
DBLband_incs_kill <- DBLband_incs %>%
  filter(week != "2019-07-07" & week != "2019-07-14") %>% 
  filter(fate == "kill") 

DBLband_incs_live <- DBLband_incs %>%
  filter(week != "2019-07-07" & week != "2019-07-14") %>%
  filter(fate == "live")

#mean by species_kill
bottom_mean <- DBLband_incs_kill %>%
  group_by(species) %>%
  summarise(bot_increment = mean(bot_increment)) 

top_mean <- DBLband_incs_kill %>%
  group_by(species) %>%
  summarise(top_increment = mean(top_increment))

both_means_kill <- bottom_mean %>%
  right_join(top_mean)

#mean by species_live
bottom_mean <- DBLband_incs_live %>%
  group_by(species) %>%
  summarise(bot_increment = mean(bot_increment)) 

top_mean <- DBLband_incs_live %>%
  group_by(species) %>%
  summarise(top_increment = mean(top_increment))

both_means_live <- bottom_mean %>%
  right_join(top_mean)

# kill plot
#ggplot(both_means_kill, aes(x=top_increment, y=bot_increment, color="red"))+
#geom_point()

ggplot(DBLband_incs_kill, aes(x=top_increment, y=bot_increment, color="kill"))+
  geom_point()

# live plot 
#ggplot(both_means_live, aes(x=top_increment, y=bot_increment, color="green"))+
#geom_point()

ggplot(DBLband_incs_live, aes(x=top_increment, y=bot_increment, color="live"))+
  geom_point() 
