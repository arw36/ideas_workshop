# Intro to Sci Programming
# Day 1, Module 2
# May 13, 2019
# Anna Willoughby

# load packages 
library(tidyverse)

#load data
wnv <- read.csv("data/wnv.csv")

# Histogram of WNV cases by state by year
g <- ggplot(wnv, aes(Year)) + 
  geom_bar(aes(weight = Total, fill = State)) + 
  facet_wrap(~ State) +
  theme(legend.position = "none")
g

# Due to skewed data, we will look at the logarithm of the case data
# Using natural log
g2 <- ggplot(wnv, aes(Year)) + 
  geom_bar(aes(weight = log(Total), fill = State)) + 
  facet_wrap(~ State) +
  theme(legend.position = "none")
g2

# using log10 
g3 <- ggplot(wnv, aes(Year)) + 
  geom_bar(aes(weight = log10(Total), fill = State)) + 
  facet_wrap(~ State) +
  theme(legend.position = "none")
g3

# Calculate raw case fatality rate by state
wnv$case.fatality <- wnv$Fatal / wnv$Total
g4 <- ggplot(wnv, aes(Year)) + 
  geom_bar(aes(weight = case.fatality, fill = State)) + 
  facet_wrap(~ State) +
  theme(legend.position = "none")
g4

# verify Total is sum of febrile cases, neuroinvasive cases, and other cases 
wnv$sum_total <-  wnv$EncephMen + wnv$Fever + wnv$Other
table(wnv$Total == wnv$sum_total)

# round to nearest dozen
wnv$case.reported <- floor(wnv$Total/12)
wnv$case.reported.error <- wnv$Total - floor(wnv$Total/12)*12
total.error <- sum(wnv$case.reported.error)
total.error

# calculate neurodisease  rate 
wnv$ndr <- wnv$EncephMen / wnv$Total

# write function to get mean and standard error ndr per country 

ndr.by.state.year <- function(state = 'Colorado', years = 1999:2007){
  # computes mean and standard error of neuroinvasive disease rate by state 
  # 
  # Args: 
  #   state: vector of state names 
  #   years: vector of years to be included in the calculation 
  #
  # Returns: 
  #   a dataframe containing a state name, mean neuroinvasive disease rate, and se 
  
  x <- wnv[wnv$State %in% state & wnv$Year %in% years,]
  y <- data.frame(state = x$State, ndr = x$ndr)
  m <- aggregate(y$ndr, by = list(y$state), FUN = mean)
  se <- aggregate(y$ndr, by = list(y$state), FUN = function(x) sd(x)/sqrt(length(x)) )
  out <- merge(m, se, by = 'Group.1')
  names(out) <- c('state', 'mean.ndr', 'se.ndr')
  return(out)
} 

disease <- ndr.by.state.year(state = c('California', 'Colorado', 'New York'))

ggplot(disease, aes(x = state, y = mean.ndr, fill = state)) + 
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin=mean.ndr-se.ndr, ymax = mean.ndr+se.ndr)) + 
  labs(x = 'State', y = 'Neuroinvasive disease rate', 
       title = 'Neuroinvasive disease rate, 1999-2007 (mean +/- se)', caption="Data from: https://disease")

states <- as.data.frame(unique(wnv$State))
names(states) <- c("State")

states[,2] <- ndr.by.state.year(states$State)[2]
states[,3] <- ndr.by.state.year(states$State)[3]

states %>% ggplot(aes(x = State, y = mean.ndr, fill = State)) + 
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin=mean.ndr-se.ndr, ymax = mean.ndr+se.ndr)) + 
  labs(x = 'State', y = 'Neuroinvasive disease rate', 
       title = 'Neuroinvasive disease rate, 1999-2007 (mean +/- se)', caption="Data from: https://disease") +
  theme(axis.text.x = element_text(angle = 45))

# Choose a longitude to get the center of the country 

center <- (max(wnv$Longitude)+min(wnv$Longitude))/2
wnv$coast <- ifelse(wnv$Longitude > center, "East", "West") # assign west and east by longitude
coast <- select(wnv, State, coast) # simplify data frame
states <- full_join(states, coast, by = "State") # merge to get ndr data
# plot 
ggplot(states, aes(x = coast, y = mean.ndr)) + 
  geom_boxplot() + geom_jitter()

# choose a latitude to get the center of the country (this time horizontal)
center_lat <- (max(wnv$Latitude)+min(wnv$Latitude))/2
# assign north or sourth based on center line 
wnv$region <- ifelse(wnv$Latitude > center_lat, "North", "South")
# simplify data frame
region <- select(wnv, State, region)
# merge in ndr data 
states <- full_join(states, region, by = "State")
# plot 
ggplot(states, aes(x = region, y = mean.ndr)) + 
  geom_boxplot() + geom_jitter()
  

# Loop over all the WNV years of data and compute 

# Total number of states reporting cases 


# total number of reported cases 
# total number of fatalities 
# case fatality rate 
# produce plots to explore how thes quanitites change over time and with respect to eachouther 
# What have you learned or suspect 


