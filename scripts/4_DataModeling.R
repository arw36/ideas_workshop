# Intro to Data Modeling
# Day 2, Module 4
# Anna Willoughby

# load packages 
library(tidyverse)
library(GGally)
library(magrittr)
library(viridis)

# Task 1: load data 
load("lyme.RData")

# Task 2: 
ggpairs(all.dfs, columns=c("prcp", "avtemp", "size", "cases"))

# Task 3: create two new columns with log10(size) and log10(cases+1)
all.dfs$size.log <- log10(all.dfs$size)
all.dfs$cases.log <- log10(all.dfs$cases + 1)
# need to add 1 to cases as some counties have 0 cases, and the log would create errors as non-real number


plot.log <- ggpairs(all.dfs, columns = c("prcp", "avtemp", "size.log", "cases.log"))
plot.log

set.seed(222)
subset <- sample_n(all.dfs, 100)

myPlot <- ggplot(subset, aes(x = prcp, y = avtemp )) + geom_point()
myPlot + geom_smooth(method = "lm")

# create a linear model for my data 
myModel <- lm(avtemp ~ prcp, data = subset)
summary(myModel)

# Task 7: What is the slope and p value? 
summary(myModel)$coefficients[2,1] # slope value
summary(myModel)$coefficients[2,4] # p value is significant

# Task 8: 
all.dfs %>% 
  group_by(year) %>% 
  summarise(total = sum(size)) %>% 
  ggplot(.) + 
  geom_point(aes(x = year, y = total))

# Task 9 + 10: create a data frame called "by_state" from the main data frame, that groups by state, and inspect is 
by_state <- all.dfs %>% group_by(state) %>% nest()

# Task 11: Display the Georgia data 
by_state$data[[10]]

# Task 12: function to return linear model from data frame
linGrowth_model <- function(df){
  lm(size ~ year, data = df)
}

models <- purrr::map(by_state$data, linGrowth_model)

#Task 13: add column with model objects 
by_state %<>% mutate(model = purrr::map(data, linGrowth_model))
library(modelr)
by_state %<>% mutate(resids = map2(data, model, add_residuals))

# Task 14: what is the structure of resids?
# The resids is a tible with a residual value stored for each data point 

# Task 15: 
# Function to extract sum for residuals of by_state, adjusting for absolute values 
sum.resid <- function(x){
  sum(abs(by_state$resids[[x]]$resid))
}

# run the value for each residual 
by_state <- tibble::rowid_to_column(by_state, "ID") # form ID to run through row numbers 
by_state %<>% mutate(totalResid = purrr::map(ID, sum.resid)) # run function on all by_state rows 

  
# Task 16: Write a function that accepts a linear model and returns the slope (model M has slope)

find.slope <- function(M){
  slope <- M$coefficients[[2]]
  return(slope)
}
by_state %<>% mutate(slope = purrr::map(model, find.slope))
by_state$slope %<>% as.character()
by_state$slope %<>% str_extract("(\\d)+")


# Task 17: Plot the growth rate (slope value) for all states 
ggplot(by_state, aes(x = state, y = as.numeric(slope))) + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Task 18: Plot the total residuals for all states
ggplot(by_state, aes(x = state, y = as.numeric(totalResid))) + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Task 19: Repeat tasks 9 & 10 using a diff data frame name, by_state2
by_state2 <- all.dfs %>% group_by(state) %>% nest

# Task 20: Write a function that accepts an element of the by_state2$data list-colmn and 
# returns the spearman correlation coefficient between Lyme disease cases and preciptiation


cortest <- function(df){
  cor.test(df$prcp, df$cases, method = "spearman")$estimate
}

by_state2 %<>% mutate(spearman = purrr::map(data, cortest))
spearmans <- unnest(by_state2, spearman) # unnest data

# plot correlation estimate of Lyme Model
ggplot(spearmans, aes(x = state, y = spearman)) + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
