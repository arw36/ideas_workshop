# Intro to Data Wrangling
# Day 2, Module 3
# Anna Willoughby

# load packages 
library(tidyverse)
library(GGally)
library(maptools)
library(maps)
library(magrittr)

# Task 1: load data 
ld <- read_csv(file = "lyme.csv")
pop <- read_csv(file = "pop.csv")
prism <- read_csv(file = "climate.csv")

# Task 2: which ways does the data fail to conform to tidy? 
# pop is not tidy as there are multiple rows for the same county, with separate date (eg. 1970-2013 population on one row, then 2014+ data on)
# pop has several errors with years being misspelled
# we want each row to have a unique row per county/fips number

# Task 3: annotate code to make pop into a tidy dataframe
pop %<>% select(fips, starts_with("pop2")) # subset columns to include fips code and population levels for the 2000s  
pop %<>% gather(starts_with("pop2"), key = "str_year", value = "size") %>% na.omit # wrangle data from wide years to long year and population size data
pop %<>% mutate(year = str_replace_all(str_year, "pop", "")) # extract years from the str_year column
pop %<>% mutate(year=as.integer(year)) # convert year from character string to an integer
pop %<>% mutate(fips = str_replace_all(fips, "^0", "")) # remove leading 0 value
pop %<>% mutate(fips=as.integer(fips)) # convert fipsom character string to an integer

# If we wanted to remove state-level data: 
# try filtering by a list method 
# states <- 1:50*1000 # create a vector of state codes
# states <- as.integer(states) # convert to integer
# counties <- 1:max(pop$fips) # create a vector of all fips code
# counties <- counties[! counties %in% states] # filter out the state codes
# pop %<>% subset(fips %in% counties) # subset by only county codes


# Task 4: Format lyme data into tidy data
# create fips code via nested conditionals
ld$fips <- ifelse(str_length(ld$CTYCODE) == 2, paste(ld$STCODE, "0", ld$CTYCODE, sep = ""), 
                                    ifelse(str_length(ld$CTYCODE) == 1, paste(ld$STCODE, "00", ld$CTYCODE, sep = ""), 
                                                      paste(ld$STCODE, ld$CTYCODE, sep = "")))# add county part of fips code 


ld %<>% gather(starts_with("Cases"), key = "str_year", value = "cases")  # wrangle data from wide years to long year and population size data
ld %<>% mutate(year = str_replace_all(str_year, "Cases", "")) # extract years from the str_year column
ld %<>% mutate(year = as.integer(year))
ld %<>% rename(state=STNAME, county = CTYNAME)

# Alternative Method via a function
# fips.builder <- function(st, ct){
#  if (str_length(ct)==3){
#    fips <- paste(as.character(st), as.character(ct), sep ="") %>% as.integer
#  }
#  else if (str_length(ct)==2){
#    fips <- paste(as.character(st), "0", as.character(ct), sep="") %>% as.integer
#  }
#  else {
#    fips <- paste(as.character(st), "00", as.character(ct), sep = "") %>% as.integer
#  }
#  return(fips)
# }

# ld %<>% rowwise() %>% mutate(fips=fips.builder(STCODE, CTYCODE)) # takes about 10 seconds 

ld %<>% select(-c(STCODE, CTYCODE, str_year))

# Task 5: join ld and prism data frame
ld.prism <- inner_join(ld, prism, by = c("year", "fips")) 

# join merged data with demographic data
all.dfs <- inner_join(ld.prism, pop, by = c("year", "fips"))

# Determine how many cases of lyme were reported each year 
cases_by_year <- ld %>% ungroup %>% group_by(year) %>% 
  summarize(total = sum(cases)) %>% arrange(desc(total))
# average number of cases in each state, averaged across county and year 
cases_by_state_year <- ld %>% ungroup %>% group_by(year, state) %>% 
  summarize(mean_cases = mean(cases)) %>% arrange(desc(mean_cases))
# What is the worst year? 
# 2009 was the worst year

# Which three states have been impacted on average 
# Conneticut, Massachusetts, New Jersey

# save data 
# save(ld, prism, pop, all.dfs, cases_by_state_year, cases_by_year, file = "lyme.RData")
# write_csv(all.dfs,"lyme_combined.csv")


# make a map 
# get map data for US counties and states 
county_map <- map_data("county")
state_map <- map_data("state")

ag.fips <- group_by(all.dfs, fips) # rename dataframe
ld.16y <- summarize(ag.fips, all.cases = sum(cases)) # sum cases by county
ld.16y <- left_join(select(all.dfs, c(state, county, fips)), ld.16y) # add state and county name
ld.16y <- distinct(ld.16y) # remove duplicates
ld.16y %<>% rename(region=state, subregion=county) # rename columnes
ld.16y$subregion <- str_replace_all(ld.16y$subregion, " County", "") # simplify county name
ld.16y$region <- tolower(ld.16y$region) # translate upper characters to lower 
ld.16y$subregion <- tolower(ld.16y$subregion) # translate upper characters to lower 
ld.16y %<>% mutate(log10cases=log10(1+all.cases)) # transform to log cases
map.ld.16y <- left_join(county_map, ld.16y) # merge to spatial data
# plot 
ggplot(map.ld.16y) + 
  geom_polygon(aes(long, lat, group = group, fill = log10cases), color= "gray", lwd = 0.2) +
  scale_colour_gradientn(colours=rev(heat.colors(10)))

