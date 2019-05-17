### R code from vignette source 'stochastic-simulation.rnw'
### Encoding: UTF-8
### Anna Willoughby
### May 16, 2019

###################################################
### code chunk number 1: stochastic-simulation.rnw:70-84
###################################################
SI.onestep <- function (x, params) {     #function for one step of the stochastic SI epidemic
   X <- x[2]                             #the second element of x is number of susceptibles X
   Y <- x[3]                             #the third element of x is number of infecteds Y
   with(
     as.list(params),
     {
       new.Y <- Y+1                      #whenever an event occurs we increase infecteds by 1...
       new.X <- X-1                      #and descrease susceptibles by 1
       
       tau <- -log(runif(1))/(params$beta*X*Y/(X+Y)) #exponential random time to next event
       c(tau=tau,X=new.X,Y=new.Y)        #store result
       }
     )
}


###################################################
### code chunk number 2: stochastic-simulation.rnw:89-98
###################################################
SI.model <- function (x, params, nstep) { #function to iterate the stochastic SI for nstep events
  output <- array(dim=c(nstep+1,3))       #set up an array to store all the results
  colnames(output) <- c("time","X","Y")   #name the variables in the array
  output[1,] <- x                         #the first record of the array is the initial condition
  for (k in 1:nstep) {                    #iterate the model for nstep events
    output[k+1,] <- x <- SI.onestep(x,params) #update x and store result
  }
  output                                  #return output
}


###################################################
### code chunk number 3: stochastic-simulation.rnw:105-122
###################################################
set.seed(38499583)                        #set the random seed so results are repeatable
nsims <- 10                               #number of simulations to run
pop.size <- 200                           #total size of the population
Y0 <- 2                                   #initial number infected
nstep <- pop.size-Y0                      #how many steps to run? until everyone infected
xstart <- c(time=0,X=(pop.size-Y0),Y=Y0)  #initial conditions
params <- list(beta=0.3) #parameters
data <- vector(mode='list',length=nsims)  #create a list called ``data'' to store all runs
for (k in 1:nsims) {                      #simulate k different runs
  data[[k]] <- as.data.frame(SI.model(xstart,params,nstep))  #main simulation step
  data[[k]]$cum.time <- cumsum(data[[k]]$time) #calculates the running sum of inter-event intervals
}
max.y<-max(data[[1]]$cum.time)               #find the maximum time any simulation ran (to set x axis)
plot(c(0,pop.size),c(0,pop.size),type='n',xlab='time',ylab='incidence',xlim=c(0,max.y)) #set up plot
for (k in 1:nsims) {                      #loop over each simulation...
  lines(Y~cum.time,data=data[[k]],col=k,type='o')  #to plot
}

# Exercise 1. Simulate the stochastic SI model using Gillespie’s direct method. 
# Experiment with the initial number of infecteds (Y0) and with the total population size (N). 
# wrap in a function for easy iterations 

plot.SImodel <- function(pop.size, Y0){
  nstep <- pop.size-Y0                      #how many steps to run? until everyone infected
  xstart <- c(time=0,X=(pop.size-Y0),Y=Y0)  #initial conditions
  params <- list(beta=0.3) #parameters
  data <- vector(mode='list',length=nsims)  #create a list called ``data'' to store all runs
  for (k in 1:nsims) {                      #simulate k different runs
    data[[k]] <- as.data.frame(SI.model(xstart,params,nstep))  #main simulation step
    data[[k]]$cum.time <- cumsum(data[[k]]$time) #calculates the running sum of inter-event intervals
  }
  max.y<-max(data[[1]]$cum.time)               #find the maximum time any simulation ran (to set x axis)
  plot(c(0,pop.size),c(0,pop.size),type='n',xlab='time',ylab='incidence',xlim=c(0,max.y)) #set up plot
  for (k in 1:nsims) {                      #loop over each simulation...
    lines(Y~cum.time,data=data[[k]],col=k,type='o')  #to plot
  }
  
}
par(mfrow = c(2,2)) # make output window a grid
plot.SImodel(10000,1) # variance, but still similar trends
plot.SImodel(50,1) #  with low starting infected, high variance and low predicatability
plot.SImodel(10000,100) # least variance, most predictable
plot.SImodel(50,8) #

# What eﬀects do these have on the predictability of the epidemic? 
# Increases on the infected individuals make the epidemic more predictable 
# I think larger population sizes increase the predictability of the 

# What eﬀects do these have on the variability of the ﬁnal outbreak size?
# A small population increases the variability.
# A small number of infected individuals to begin make the epidemic more variable

###################################################
### code chunk number 4: stochastic-simulation.rnw:138-162
###################################################
SIR.onestep <- function (x, params) {  #function to calculate one step of stochastic SIR
  X <- x[2]                            #local variable for susceptibles
  Y <- x[3]                            #local variable for infecteds
  Z <- x[4]                            #local variable for recovereds
  N <- X+Y+Z                           #total population size (subject to demographic change)
  with(                                #use with as in deterministic model to simplify code
       as.list(params), 
       {
         rates <- c(mu*N, beta*X*Y/N, mu*X, mu*Y, gamma*Y, mu*Z)
         changes <- matrix(c( 1, 0, 0,
                             -1, 1, 0,
                             -1, 0, 0,
                              0,-1, 0,
                              0,-1, 1,
                              0, 0,-1),
                           ncol=3, byrow=TRUE)
         tau <- -log(runif(1)) / sum(rates)   # exponential waiting time
         U <- runif(1)       #uniform random deviate
         m <- min(which(cumsum(rates)>=U*sum(rates)))
         x <- x[2:4] + changes[m,]
         return(out <- c(tau, x))
       }
       )
}


###################################################
### code chunk number 5: stochastic-simulation.rnw:167-177
###################################################

SIR.model <- function (x, params, nstep) {  #function to simulate stochastic SIR
  output <- array(dim=c(nstep+1,4))         #set up array to store results
  colnames(output) <- c("time","X","Y","Z") #name variables
  output[1,] <- x                           #first record of output is initial condition
  for (k in 1:nstep) {                      #iterate for nstep steps
    output[k+1,] <- x <- SIR.onestep(x,params)
  }
  output                                    #return output
}


###################################################
### code chunk number 6: stochastic-simulation.rnw:182-201
###################################################
set.seed(38499583)                #set seed
nsims <- 10                       #number of simulations
pop.size <- 100                   #total population size
Y0 <- 8                           #initial number infected
X0 <- round(0.9*pop.size)        #initial number suscepitlble (~90% of population)
nstep <- 1600                     #number of events to simulate
xstart <- c(time=0,X=X0,Y=Y0,Z=pop.size-X0-Y0) #initial conditions
params <- list(mu=0.00001,beta=60,gamma=365/13) #parameters
data <- vector(mode='list',length=nsims) #initialize list to store the output
for (k in 1:nsims) {              #simulate nsims times
  data[[k]] <- as.data.frame(SIR.model(xstart,params,nstep))
  data[[k]]$cum.time <- cumsum(data[[k]]$time)
}
max.time<-data[[1]]$cum.time[max(which(data[[1]]$Y>0))] #maximum time in first simulation
max.y<-1.8*max(data[[1]]$Y)       #find max infected in run 1 and increase by 80% for plot
plot(Y~cum.time,data=data[[1]],xlab='Time',ylab='Incidence',col=1,xlim=c(0,max.time),ylim=c(0,max.y))
for (k in 1:nsims) {              #add multiple epidemics to plot
  lines(Y~cum.time,data=data[[k]],col=k,type='o')
}

# Exercise 2. Simulate the stochastic SIR model using Gillespie’s direct method. 
# As before, experiment with the initial number of infecteds (Y0) and with the 
# total population size (N). 

# make a function that wraps the plot 
plot.smodel <- function(Y0, pop.size){
  X0 <- round(0.9*pop.size)        #initial number suscepitlble (~90% of population)
  nstep <- 1600                     #number of events to simulate
  xstart <- c(time=0,X=X0,Y=Y0,Z=pop.size-X0-Y0) #initial conditions
  params <- list(mu=0.00001,beta=60,gamma=365/13) #parameters
  data <- vector(mode='list',length=nsims) #initialize list to store the output
  for (k in 1:nsims) {              #simulate nsims times
    data[[k]] <- as.data.frame(SIR.model(xstart,params,nstep))
    data[[k]]$cum.time <- cumsum(data[[k]]$time)
    }
  max.time<-data[[1]]$cum.time[max(which(data[[1]]$Y>0))] #maximum time in first simulation
  max.y<-1.8*max(data[[1]]$Y)       #find max infected in run 1 and increase by 80% for plot
  plot(Y~cum.time,data=data[[1]],xlab='Time',ylab='Incidence',col=1,xlim=c(0,max.time),ylim=c(0,max.y))
  for (k in 1:nsims) {              #add multiple epidemics to plot
    lines(Y~cum.time,data=data[[k]],col=k,type='o')
    }
} 

# experimenting with different variables 
par(mfrow=(c(2,2)))
plot.smodel(200,8)
plot.smodel(50,20)
plot.smodel(50, 1)
plot.smodel(10000,20)
plot.smodle(10000,1)

# What eﬀects do these have on the predictability of the epidemic?
# Decreasing the population size decreases the predictability of the epidemic.
# Decreasing the infected decreaseed the predictability of the epidemic. 


# I personally really hate lists, as I think it adds a layer of complexity and doesn't conform to tidy guidelines
# To remedy this, and have more flexibility by using ggplot, I'd like to apend to a long data frame, with unique 
# identifiers of "time.step" and "simulation". 

library(ggplot2) # format to ggplot for easier faceting
library(tidyr) # for transforming data 
library(dplyr) # for transforming data 
library(stringr) # for text parsing

# wrap the plot into a function for quick iterations
# this funtion doesn't work yet! 
# plot.smodel <- function(pop.size, Y0){
  X0 <- round(0.9*pop.size)        #initial number suscepitlble (~90% of population)
  nstep <- 1600                     #number of events to simulate
  xstart <- c(time=0,X=X0,Y=Y0,Z=pop.size-X0-Y0) #initial conditions
  params <- list(mu=0.00001,beta=60,gamma=365/13) #parameters
  data <- vector(mode='list',length=nsims) #initialize list to store the output
  for (k in 1:nsims) {              #simulate nsims times
    data[[k]] <- as.data.frame(SIR.model(xstart,params,nstep))
    data[[k]]$cum.time <- cumsum(data[[k]]$time)
    }
  max.time<-data[[1]]$cum.time[max(which(data[[1]]$Y>0))] #maximum time in first simulation
  max.y<-1.8*max(data[[1]]$Y)       #find max infected in run 1 and increase by 80% for plot
  # create column for each simulation and time step
  for (k in data[[k]]) {
    data[[k]]$time.step <- seq(0, nstep, by = 1)
  }
  # create column for simulation number 
  for (k in data[[k]]) {
    data[[k]]$simulation <- str_extract(colnames(data[[k]][1]), "\\d+")
  }
  for (k in 1:nsims)  
    ggplot(data[[1]], aes(x = cum.time, y = Y)) + # set up ggplot 
            geom_point() + 
            labs(x = 'Time',y ='Incidence',col=1,xlim=c(0,max.time),ylim=c(0,max.y))
  {              #add multiple epidemics to plot
  lines(Y~cum.time,data=data[[k]],col=k,type='o')
    }
}
