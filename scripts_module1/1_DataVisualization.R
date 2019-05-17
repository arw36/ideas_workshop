# Intro to Data Visualization
# Day 1, Module 1
# Anna Willoughby

# load packages 
library(lubridate)
library(tidyverse)
library(ggridges)

#load data
mers <- read.csv('data/cases.csv')


# correct errors
mers$hospitalized[890] <- c('2015-02-20')
mers <- mers[-417,]
mers$onset2 <- ymd(mers$onset)
mers$hospitalized2 <- ymd(mers$hospitalized)

# calculate days from onset of epidemic 
day0 <- min(na.omit(mers$onset2)) 
# Q1 would get errors for na, as many rows without onset date
mers$epi.day <- as.numeric(mers$onset2 - day0) # create epidemic day value 
#Q2 changes the class from diff in time to a numeric

# Plot 
ggplot(data=mers) + 
  geom_bar(aes(x=epi.day, fill = country), position = "fill") + 
  labs(x  = 'Epidemic day', y = 'Case count', 
  title = 'Global count of MERS cases by date of symptom onset', 
  caption = "Data from: https://github.com/rambaut/MERS-Cases/blob/gh-pages/data/cases.csv") + 
  coord_flip() 

# Exercise 1: Position = "fill" stacks so you can see the percent from each country, can see which country the % of cases came from 
# Exercise 2: flips the axes to emphasize the case counts

mers$infectious.period <- mers$hospitalized2-mers$onset2 #calculate "raw" infectious period 
mers$infectious.period <- as.numeric(mers$infectious.period, units = "days") # convert to days
ggplot(data = mers) + 
  geom_histogram(aes(x = infectious.period)) +
  labs(x = 'Infectious period', y = 'Frequency', 
       title = 'Distribution of calculated MERS infectious period', 
       caption="Data from: https://github/rambaut/MERS-cases/blob/gh-pages/data/cases.csv")
# Histogram plot
mers$infectious.period2 <- ifelse(mers$infectious.period < 0 , 0, mers$infectious.period)
ggplot(data = mers)+
  geom_histogram(aes(x=infectious.period2)) +
  labs(x = 'Infectious period', y = 'Frequency', 
       title = 'Distribution of calculated MERS infectious period (positive values only)', 
       caption = "Data from: https://github/rambaut/MERS-cases/blob/gh-pages/data/cases.csv")
# Density plot
ggplot(data = mers)+
  geom_density(aes(x=infectious.period2)) +
  labs(x = 'Infectious period', y = 'Frequency', 
       title = 'Distribution of calculated MERS infectious period (positive values only)', 
       caption = "Data from: https://github/rambaut/MERS-cases/blob/gh-pages/data/cases.csv")
# Area plot
ggplot(data = mers)+
  geom_area(stat = 'bin', aes(x=infectious.period2)) +
  labs(x = 'Infectious period', y = 'Frequency', 
       title = 'Distribution of calculated MERS infectious period (positive values only)', 
       caption = "Data from: https://github/rambaut/MERS-cases/blob/gh-pages/data/cases.csv")
# Exercise 3 
ggplot(mers, aes(x=infectious.period2, fill = country))+
  geom_dotplot(binaxis = "x", stackdir = "up", dotsize = 0.5) +
  labs(x = 'Infectious period', y = 'Frequency', 
       title = 'Distribution of calculated MERS infectious period (positive values only)', 
       caption = "Data from: https://github/rambaut/MERS-cases/blob/gh-pages/data/cases.csv")

# Exercise 4 
ggplot(mers, aes(y = infectious.period2, x = epi.day)) +
         geom_point(aes(colour = country))


# Exercise 5, add a smooth curve fit for the total data
ggplot(mers, aes(y = infectious.period2, x = epi.day)) +
  geom_point(aes(colour = country)) + geom_smooth(method = "loess")

# Exercise 6, plot infect period against time with smooth by country
ggplot(mers, aes(y = infectious.period2, x = epi.day)) +
  geom_point(aes(colour = country)) + geom_smooth(aes(group = country, colour = country), method = "loess")

# Faceting
ggplot(mers, aes( x = epi.day, y = infectious.period2)) +
  geom_point(aes(colour = country)) + 
  facet_wrap(~country) + 
  scale_y_continuous(limits = c(0,50)) +
  labs(x = 'Epidemic Day', y = 'Infectious period', 
       title = 'MERS infectious period (positive values only) over time', 
       caption = "Data from: https://github/rambaut/MERS-cases/blob/gh-pages/data/cases.csv")

ggplot(subset(mers, gender %in% c('M', 'F') & country %in% c('KSA', 'Oman', 'Iran', 'Jordan', 'Qatar', 'South Korea', 'UAE'))) +
  geom_point(aes(x = epi.day, y = infectious.period2, colour = country)) +
  facet_grid(gender ~ country) + 
  scale_y_continuous(limits = c(0,50)) +
  labs(x = 'Epidemic Day', y = 'Infectious period', 
        title = 'MERS infectious period by gender and country', 
        caption = "Data from: https://github/rambaut/MERS-cases/blob/gh-pages/data/cases.csv")

# Exercise 7, study variation in case fatality rate (the fraction of cases that end in death) over time and across countries 
mers$case_died <- ifelse(mers$outcome %in% c("?fatal", "fatal", "fatal?"), "Yes", "No")


cf_plot <- ggplot(mers) +
  geom_bar(aes(x = epi.day, colour = case_died), position = "fill") + 
  facet_wrap(~country) + 
  scale_y_continuous(limits = c(0,1)) +
  scale_colour_manual(values = c( "blue", "red")) +
  labs(x = 'Epidemic Day', y = 'Case fatality', 
       title = 'MERS case fatality over time by country', 
       caption = "Data from: https://github/rambaut/MERS-cases/blob/gh-pages/data/cases.csv")
cf_plot


mers_country <- mers %>% group_by(country, epi.day, case_died) %>% summarise(cases = n())
mers_country <- spread(mers_country, key = case_died, value = cases)
mers_country$Yes <- ifelse(is.na(mers_country$Yes) == TRUE, 0, mers_country$Yes)
mers_country$No <- ifelse(is.na(mers_country$No) == TRUE, 0, mers_country$No)
mers_country$case.fatality <- mers_country$Yes / (mers_country$Yes + mers_country$No)

# Exercise 8: 
cf2_plot <- ggplot(mers_country, aes(x = epi.day, y = country, fill = case.fatality)) + 
  geom_density_ridges() +
  scale_fill_continuous() # values are not currently colored
cf2_plot

library(ggridges)
cf_plot2 <- ggplot(mers, aes(x = epi.day, y = country, fill = case_died)) + 
  geom_density_ridges_gradient() +
  scale_fill_manual(values = c("blue", "red"))
cf_plot2 

# Use plotly
library(plotly)
epi.curve <- ggplot(mers)+ 
  geom_bar(aes(x = epi.day)) +
  labs(x = 'Epidemic day', y = 'Case count', title ='Global count of MERS cases by date of symptoms onset', 
       caption = "Data from: https://github/rambaut/MERS-cases/blob/gh-pages/data/cases.csv")
ggplotly(epi.curve)

# Exercise 9 
ggplotly(cf_plot)
