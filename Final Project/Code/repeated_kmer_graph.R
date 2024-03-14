# Julia Fairbank, Mia Tarantola, Caroline Cutter
# CSCI 0321A - Bioinformatics
# Final Project - Repeated K-Mer Data
# Monday, May 16th, 2022



# -------------------------------- CODE -----------------------------------

# IMPORTS

library("ggplot2")                                   # ggplot2 library
library("tidyverse")                                 # tidyverse library
library("gridExtra")                                 # gridExtra library
library("dplyr")


setwd("~/Desktop/Final Project")


data <- read_csv("Repeated K-Mer Data.csv")                    # data file




data %>%
  ggplot(aes(x = K, y = Repeats)) +
  geom_col(color = "darkviolet", fill = "darkviolet") +
  theme_bw() +
  scale_x_continuous(breaks=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) + 
  scale_y_continuous(limits = c(0, 5000)) +
  labs(title = "Number of Repeated Sequences of Length K" ) +
  theme(plot.title=element_text(hjust=0.5))



data %>%
  ggplot(aes(x = K, y = Repeats)) +
  geom_line(color = "darkviolet") +
  #geom_smooth(method = "loess", se = F, color = "darkviolet") + 
  theme_bw() +
  scale_x_continuous(breaks=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) + 
  scale_y_continuous(limits = c(0, 5000)) 


data %>%
  ggplot(aes(x = K, y = Repeats)) +
  geom_point(color = "black") +
  geom_smooth(method = "loess", se = F, color = "darkviolet") + 
  theme_bw() +
  scale_x_continuous(breaks=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) + 
  scale_y_continuous(limits = c(0, 5000)) 


