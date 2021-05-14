# this script calcuates the increments of upper (top) and lower (Bottom) dendrometer 
# bands that were meant to capture potential swelling in the stem immediately above 
# the severed phloem (stem-girdle), and runs a series of t-tests to assess 
# significant differences

#set library to dyplr
library(dplyr)
library(tidyr)
library(ggplot2)

#importing csv
# edit path (...) to access data file
DBLband <- read.csv(".../DBL_dendrobands.csv")

DBLband_data <- select(DBLband, tag, fate, species, 
                       DBH_cm, bottom..in., top, date)

# pipeline to raw difference between upper (top) and lower (bottom) bands 
DBLband_incs <- DBLband_data %>%
  group_by(tag) %>%
  filter(!is.na(top)) %>% 
  filter(!is.na(bottom..in.)) %>% 
  mutate(bot_increment= bottom..in.-lag(bottom..in., 
                                    default = first(bottom..in.))) %>%
  mutate(top_increment= top-lag(top, default = first(top))) %>% 
  mutate(t_b_diff = bot_increment - top_increment)
  
  
#all increment_percent < 0 = 0
DBLband_incs[DBLband_incs < 0] <- 0

# need two-sample t-test for each species on each week 

# but first I need to convert date from factor to date and create a week column 
DBLband_incs$date <- as.Date(DBLband_incs$date, "%Y-%m-%d")
DBLband_incs$week <- as.Date(cut(DBLband_incs$date, breaks = "week", start.on.monday = FALSE))

# reformat data to fit a t-test
DBLband_transform <- DBLband_incs %>% 
  gather(bot_or_top, increment, bot_increment, top_increment) %>% 
  ungroup()

#make a t test funcion to compare top and bottom differences 

t_test_diff <- function(SPECIES){
  filtered <- filter(DBLband_incs, species == SPECIES & 
                       week != "2019-07-07" & week != "2019-07-14")
  output <- t.test(t_b_diff ~ fate, data = filtered, paired = FALSE)
  return(output)
}

#t_test runs a t-test comparing top and bottom bands of a given species and fate
t_test <- function(SPECIES, FATE){
  filtered <- filter(DBLband_transform, species == SPECIES & fate == FATE & 
                     week != "2019-07-07" & week != "2019-07-14")
  output <- t.test(increment ~ bot_or_top, data = filtered, paired = TRUE)
  return(output)
}

# swell calculates the percent swelling in girdled trees
swell <- function(SPECIES){
  filter_top <- filter(DBLband_transform, species == SPECIES & fate == "kill" & 
                      bot_or_top == "top_increment" & week != "2019-07-07" & 
                      week != "2019-07-14")
  filter_bot <- filter(DBLband_transform, species == SPECIES & fate == "kill" & 
                       bot_or_top == "bot_increment" & week != "2019-07-07" & 
                       week != "2019-07-14")
  mean_top <- mean(filter_top$increment)
  mean_bot <- mean(filter_bot$increment)
  percent_diff <- 1 - mean_top/mean_bot
  return(percent_diff)
}

# t-tests asking if the top and bottom bands are different
t_test("ACRU", "live"); t_test("ACRU", "kill"); swell("ACRU")

t_test("ACSA", "live"); t_test("ACSA", "kill"); swell("ACSA")

t_test("PIST", "live"); t_test("PIST", "kill"); #NS

t_test("POGR", "live"); t_test("POGR", "kill"); #NS

t_test("QURU", "live"); t_test("QURU", "kill"); #NS

# t-tests asking if the difference in top and bottom bands is different in the live 
# and dead trees 
t_test_diff("ACRU")
t_test_diff("ACSA")
t_test_diff("PIST")
t_test_diff("POGR")
t_test_diff("QURU")
 