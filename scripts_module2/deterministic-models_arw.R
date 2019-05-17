### R code from vignette source 'deterministic-models.rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: deterministic-models.rnw:104-105
###################################################
require(deSolve)                          #deSolve library needed for this computing session


###################################################
### code chunk number 2: deterministic-models.rnw:109-124
###################################################
sir.model.closed <- function (t, x, params) {    #here we begin a function with three arguments
  S <- x[1]                               #create local variable S, the first element of x
  I <- x[2]                               #create local variable I
  R <- x[3]                               #create local variable R
  with(                                   #we can simplify code using "with"
       as.list(params),                   #this argument to "with" lets us use the variable names
       {                                  #the system of rate equations
         dS <- -beta*S*I
         dI <- beta*S*I-gamma*I
         dR <- gamma*I
         dx <- c(dS,dI,dR)                #combine results into a single vector dx
         list(dx)                         #return result as a list
       }
       )
}


###################################################
### code chunk number 3: deterministic-models.rnw:133-136
###################################################
times <- seq(0,120,by=5)                    #function seq returns a sequence
params <- c(beta=0.3,gamma=1/7)             #function "c" combines values into a vector. these are per day 
xstart <- c(S=9999/10000,I=1/10000,R=0)     #initial conditions


###################################################
### code chunk number 4: deterministic-models.rnw:140-141
###################################################
out <- as.data.frame(ode(xstart,times,sir.model.closed,params))  #result stored in dataframe


###################################################
### code chunk number 5: deterministic-models.rnw:144-149
###################################################
op <- par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))                  #set graphical parameters
plot(I~time,data=out,type='b')                              #plot the I variable against time
par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)                  #re-set graphical parameters
plot(I~S,data=out,type='b',yaxt='n',xlab='S')               #plot phase portrait
par(op)                                                     #re-set graphical parameters

# Exercise 1. Explore the dynamics of the system for different values of the 
# beta and lambda parameters by simulating and plotting trajectories as time series and in phase space (e.g., I vs. S). 


# created a function that produces plot based on beta and gam
plot_I <- function(b, g){
  params <- c(beta=b,gamma=g)  
  out1 <- as.data.frame(ode(xstart, times, sir.model.closed, params))
  op <- par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))
  plot(I~time,data=out1,type='b')               #plot the I variable against time
  par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)                  #re-set graphical parameters
  plot(I~S,data=out1,type='b',yaxt='n',xlab='S')               #plot phase portrait
  par(op)                                           
}


# looping for a variety of diff beta and gamma values
# this would likely be better as a matrix of values then facetted or plotted on top of eachother
b_seq <- seq(0.1,0.9, by =0.1)
g_seq <- seq(1/7, 1, by = 0.1)
library(tidyverse)
I.plot <- purrr::map2(b_seq, g_seq, plot_I)

# Explore with ggplot 

ggplot(out, aes(x = time,y = I)) + geom_line(colour = viridis(params[[1]])) 


# Exercise 2. Explore the dynamics of the system for one set of beta and gamma at different initial conditions.

# beta and gamma params of choice 

params <- c(beta=0.7,gamma=1/7)   

plot_diffstart <- function(population, infected){
  xstart <- c(S=(population - infected)/population,I= infected/population,R=0) 
  out1 <- as.data.frame(ode(xstart, times, sir.model.closed, params))
  op <- par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))
  plot(I~time,data=out1,type='b')                     #plot the I variable against time
  par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)          #re-set graphical parameters
  plot(I~S,data=out1,type='b',yaxt='n',xlab='S')      #plot phase portrait
  par(op)                                           
} 

#initial conditions
plot_diffstart(10000,1) # this should match the original example
plot_diffstart(100000,1) # a much bigger pop 
# in this case the epidemic takes off a bit later 
plot_diffstart(100,1) # a much smaller pop
# in this case the epidemic starts much earlier 

# Exercise 3. Modify the codes given to study a demographically open SIR model

# add in new parameters for birth and death rates
params <- c(beta=0.3,gamma=1/7, b = 0.02, d = 0.015)
sir.model.open <- function (t, x, params) {    #here we begin a function with three arguments
  S <- x[1]                               #create local variable S, the first element of x
  I <- x[2]                               #create local variable I
  R <- x[3]                               #create local variable R
  with(                                   #we can simplify code using "with"
       as.list(params),                   #this argument to "with" lets us use the variable names
       {                                  #the system of rate equations
         dS <- -beta*S*I + b - d
         dI <- beta*S*I-gamma*I - d
         dR <- gamma*I - d
         dx <- c(dS,dI,dR)                #combine results into a single vector dx
         list(dx)                         #return result as a list
       }
       )
}

# plot with the open SIR system
xstart <- c(S=9999/10000,I=1/10000,R=0)     #initial conditions# remember to redo the initial conditions 
plot.bd <- function(b, d){
  params <- c(beta=0.3,gamma=1/7, b = b, d = d)  
  out1 <- as.data.frame(ode(xstart, times, sir.model.open, params))
  op <- par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))
  plot(I~time,data=out1,type='b')               #plot the I variable against time
  par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)                  #re-set graphical parameters
  plot(I~S,data=out1,type='b',yaxt='n',xlab='S')               #plot phase portrait
  par(op)                                           
}
plot.bd(0.02, 0.015)

# Exercise 4. Modify the codes given to study a SEIR model 

seir.model.closed <- function (t, x, params) {    #here we begin a function with three arguments
  S <- x[1]                               #create local variable S, the first element of x
  E <- x[2]                               #create local variable E
  I <- x[3]                               #create local variable I
  R <- x[4]                               #create local variable R
  with(                                   #we can simplify code using "with"
       as.list(params),                   #this argument to "with" lets us use the variable names
       {                                  #the system of rate equations
         dS <- -beta*S*I
         dE <- beta*S*I - sigma*E
         dI <- sigma*E-gamma*I
         dR <- gamma*I
         dx <- c(dS,dE,dI,dR)                #combine results into a single vector dx
         list(dx)                         #return result as a list
       }
       )
}




