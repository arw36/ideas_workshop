# load library
library(tidyverse)
# load data
cases <- read.csv(file = "cases.csv")

# make a plot 
ggplot(cases[grepl(cases$age, "?") == FALSE, ], aes(x = as.numeric(age))) + 
  geom_histogram() + 
  labs(x = "Patient age")