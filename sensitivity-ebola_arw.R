### R code from vignette source 'sensitivity-ebola.rnw'
### Encoding: UTF-8
### Anna Willoughby 
### May 17, 2019


###################################################
### code chunk number 1: sensitivity-ebola.rnw:86-121
###################################################
# Exercise 1: write a function to return the rates of change for this model. 

legrand <- function(t, x, params){

  # states
  S <- x[1]
  E <- x[2]
  I <- x[3]
  H <- x[4]
  F <- x[5]
  R <- x[6]
  N <- S+E+I+H+F+R
  
  # calculated parameters
  gammaih = params$gammai*params$gammah/(params$gammah - params$gammai)
  gammadh = params$gammad*params$gammah/(params$gammah - params$gammad)
  # rates
  r1 <- (params$betaI*S*I + params$betaH*S*H + params$betaF*S*F)/N
  r2 <- params$alpha*E
  r3 <- params$gammah*params$theta1*I
  r4 <- gammadh*params$delta2*H
  r5 <- params$gammaf*F
  r6 <- params$gammai*(1-params$theta1)*(1-params$delta1)*I
  r7 <- params$delta1*(1-params$theta1)*params$gammad*I
  r8 <- gammaih*(1-params$delta2)*H
  
  # derivatives
  dS <- -r1
  dE <- r1 - r2
  dI <- r2 - r3 - r7 - r6
  dH <- r3 - r4 - r8
  dF <- r4 + r7 - r5
  dR <- r5 + r6 + r8
    
  #output
  out <- list(c(dS, dE, dI, dH, dF, dR))
}


require(deSolve)

times <- seq(0, 52, by=1)  #solve for 52 weeks
params <- list(betaI=2.532,
               betaH=0.012,
               betaF=0.462,
               alpha=7/12,
               gammah=7/4.2,
               gammai=7/10,
               gammaf=7/2,
               gammad = 7/8,
               theta1=0.65,
               delta1=0.47,
               delta2=0.42)
pop.size <- 470000
I0 <- 9
xstart <- c(S=pop.size-I0, E=0, I=I0 , H=0, F=0, R=0)
out <- as.data.frame(ode(xstart, times, legrand, params))

# Exercise 2: plot the six state variables 

plot.legrand <- function(out){
  par(mfrow=c(3,2))
  plot(out$time, out$S, xlab='Time', ylab='Susceptible', type='l', col='steelblue4', lwd=3)
  plot(out$time, out$E, xlab='Time', ylab='Exposed', type='l', col='steelblue4', lwd=3)
  plot(out$time, out$I, xlab='Time', ylab='Infectious', type='l', col='steelblue4', lwd=3)
  plot(out$time, out$H, xlab='Time', ylab='Hospitalized', type='l', col='steelblue4', lwd=3)
  plot(out$time, out$F, xlab='Time', ylab='Funeral', type='l', col='steelblue4', lwd=3)
  plot(out$time, out$R, xlab='Time', ylab='Removed', type='l', col='steelblue4', lwd=3)
}

plot.legrand(out)

# Exercise 3: Modify the model to allow interventions at time T


# model of Legrand et al (2007) with interventions
legrand2 <- function(t, x, params){

  # states
  S <- x[1]
  E <- x[2]
  I <- x[3]
  H <- x[4]
  F <- x[5]
  R <- x[6]
  N <- S+E+I+H+F+R
  
  # calculated parameters
  gammaih = params$gammai*params$gammah/(params$gammah - params$gammai)
  gammadh = params$gammad*params$gammah/(params$gammah - params$gammad)
  
  # rates
  r1 <- (params$betaI*S*I*(t<T) + params$betaI*S*I*(1-params$z)*(t>=params$T) + params$betaH*S*H*(t<params$T) + params$betaF*S*F*(t<params$T))/N
  r2 <- params$alpha*E
  r3 <- params$gammah*params$theta1*I
  r4 <- gammadh*params$delta2*H
  r5 <- params$gammaf*F
  r6 <- params$gammai*(1-params$theta1)*(1-params$delta1)*I
  r7 <- params$delta1*(1-params$theta1)*params$gammad*I
  r8 <- gammaih*(1-params$delta2)*H
  
  # derivatives
  dS <- -r1
  dE <- r1 - r2
  dI <- r2 - r3 - r7 - r6
  dH <- r3 - r4 - r8
  dF <- r4 + r7 - r5
  dR <- r5 + r6 + r8
    
  #output
  out <- list(c(dS, dE, dI, dH, dF, dR))
}


# Here is and intervention at T = 7 
params2 <- list(betaI=2.532,
               betaH=0.012,
               betaF=0.462,
               alpha=7/12,
               gammah=7/4.2,
               gammai=7/10,
               gammaf=7/2,
               gammad = 7/8,
               theta1=0.65,
               delta1=0.47,
               delta2=0.42,
               T=7,
               z=0.88)
pop.size <- 470000
I0 <- 9
xstart <- c(S=pop.size-I0, E=0, I=I0 , H=0, F=0, R=0)
out2 <- as.data.frame(ode(xstart, times, legrand2, params2))

plot.legrand(out2)

get.size <- function(out) tail(out$R,1)

# Summarize paramteterization with the total epidemic size 
output0 <- rep(NA, 15)
intervention.times <- seq(1,15)
for(T in intervention.times){
  params2$T <- T
  out3 <- as.data.frame(ode(xstart, times, legrand2, params2))
  output0[T] <- get.size(out3)  
}
par(mfrow = c(1,1))
plot(intervention.times, output0, log='y', xlab='Intervention time', ylab='Epidemic size', las=1)

# Latin Hypercube Sampling 
require(lhs)            #add the lhs library
h <- 1000               #choose number of points
set.seed(6242015)
lhs<-maximinLHS(h,12)   #simulate

betaI.min <- 1
betaI.max <- 4
betaH.min <- 0.01
betaH.max <- 0.5
betaF.min <- 0.1
betaF.max <- 4
alpha.min <- 7/21
alpha.max <- 7/2
gammah.min <- 7/2
gammah.max <- 7/1
gammai.min <- 7/21
gammai.max <- 7/2
gammaf.min <- 7/7
gammaf.max <- 7/1
gammad.min <- 7/21
gammad.max <- 7/2
theta1.min <- 0
theta1.max <- 1
delta1.min <- 0
delta1.max <- 1
delta2.min <- 0
delta2.max <- 1
z.min <- 0
z.max <- 1


# create a parameter set by rescaling lhs
params.set <- cbind(
  betaI = lhs[,1]*(betaI.max-betaI.min)+betaI.min,
  betaH = lhs[,2]*(betaH.max-betaH.min)+betaH.min,
  betaF = lhs[,3]*(betaF.max-betaF.min)+betaF.min,
  alpha = lhs[,4]*(alpha.max-alpha.min)+alpha.min,
  gammah = lhs[,5]*(gammah.max-gammah.min)+gammah.min,
  gammai = lhs[,6]*(gammai.max-gammai.min)+gammai.min,
  gammaf = lhs[,7]*(gammaf.max-gammaf.min)+gammaf.min,
  gammad = lhs[,8]*(gammad.max-gammad.min)+gammad.min,  
  theta1 = lhs[,9]*(theta1.max-theta1.min)+theta1.min,  
  delta1 = lhs[,10]*(delta1.max-delta1.min)+delta1.min,  
  delta2 = lhs[,11]*(delta2.max-delta2.min)+delta2.min,    
  z = lhs[,12]*(z.max-z.min)+z.min)  

# How many levels of T to consider? 
levels <- 15
# Only sample a fraction of points, in this case 25%
h2 <-250
# create nested loop for different values of T and different parameter values
j <- 1  
data <- data.frame(matrix(rep(NA,levels*h2*14),nrow=levels*h2))
for(i in 1:h2){
  for (T in intervention.times){
    
    data[j,1:13] <- params <- as.list(c(params.set[i,], T=T))
    out <- as.data.frame(ode(xstart, times, legrand2, params))
    data[j,14] <- get.size(out)
    j <- j+1
    
  }
} # takes about a minute 
names(data) <- c(names(params),'outbreak.size')
save(data, file='data/sens_data.Rdata')

load('data/sens_data.Rdata')
plot(intervention.times, output0, type='l', lwd=3, ylim=c(10,2e6), log='y',
            xlab='Intervention times',
            ylab='Epidemic size')
points(data$T, data$outbreak.size, pch=19, cex=0.3, col='blue')

# summarise with box plots 
boxplot(data$outbreak.size~data$T, ylim=c(10,2e6), border='blue', log='y', pch='.')
par(new=TRUE)
plot(intervention.times, output0, type='l', lwd=3, ylim=c(10,2e6), log='y', xlab='Intervention times', ylab='Epidemic size', axes=FALSE)

# Exercise 4: write a function to see which variable is most highly correlated with infections prevented
c.matrix <- cor(data)
library(magrittr)
library(dplyr)
c.df <- as.data.frame(c.matrix)
c.df <- tibble::rownames_to_column(c.df, "parameter")
c.df %<>% 
  dplyr::select(parameter, outbreak.size)
# theta 1 is most negatively correlated with outbreak size, so positively correlated with cases prevented
# this is the fraction of infectious that go to the hospital (and thus then get the intervention)

# Exercise 5: Give the most sensitive variable a new lower or upper bound, rerun and determine what
# should be the parameter for future study 
# reset theta1
theta1.min <- 0.5
theta1.max <- 1

# rebind parameters 
params.set2 <- cbind(
  betaI = lhs[,1]*(betaI.max-betaI.min)+betaI.min,
  betaH = lhs[,2]*(betaH.max-betaH.min)+betaH.min,
  betaF = lhs[,3]*(betaF.max-betaF.min)+betaF.min,
  alpha = lhs[,4]*(alpha.max-alpha.min)+alpha.min,
  gammah = lhs[,5]*(gammah.max-gammah.min)+gammah.min,
  gammai = lhs[,6]*(gammai.max-gammai.min)+gammai.min,
  gammaf = lhs[,7]*(gammaf.max-gammaf.min)+gammaf.min,
  gammad = lhs[,8]*(gammad.max-gammad.min)+gammad.min,  
  theta1 = lhs[,9]*(theta1.max-theta1.min)+theta1.min,  
  delta1 = lhs[,10]*(delta1.max-delta1.min)+delta1.min,  
  delta2 = lhs[,11]*(delta2.max-delta2.min)+delta2.min,    
  z = lhs[,12]*(z.max-z.min)+z.min)  
# rerun simulation
for(i in 1:h2){
  for (T in intervention.times){
    
    data[j,1:13] <- params <- as.list(c(params.set2[i,], T=T))
    out <- as.data.frame(ode(xstart, times, legrand2, params))
    data[j,14] <- get.size(out)
    j <- j+1
    
  }
} # takes about a minute 
names(data) <- c(names(params),'outbreak.size')
save(data, file='data/sens_Tdata.Rdata')

load('data/sens_Tdata.Rdata')
plot(intervention.times, output0, type='l', lwd=3, ylim=c(10,2e6), log='y',
     xlab='Intervention times',
     ylab='Epidemic size')
points(data$T, data$outbreak.size, pch=19, cex=0.3, col='blue')

boxplot(data$outbreak.size~data$T, ylim=c(10,2e6), border='blue', log='y', pch='.')
par(new=TRUE)

c.matrix <- cor(data)
c.df <- as.data.frame(c.matrix)
c.df <- tibble::rownames_to_column(c.df, "parameter")
c.df %<>% 
  dplyr::select(parameter, outbreak.size)

# Based on rerunning with a much higher theta1 variable (0.5 - 1). 
# This changes the predicted epidemic size to lower than previous and makes the relationship 
# between theta1 and cases prevented even stronger. 
# I think gammaI is the next variable of interest (how do you increase recovery time?)
# This is particularly promising, as theta1 and gammaI are not strongly correlated

library(sensitivity)
bonferroni.alpha <- 0.05/12
prcc <- pcc(data[,1:12], data[,14], nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)
save(prcc, file='prcc.Rdata')

load('prcc.Rdata')
summary <- print(prcc)
par(mar=c(7,4,4,2)+0.1)
plot(summary$original, main='Partial rank correlation coefficients', ylim=c(-1,1),
     xlab='', ylab='Coefficient',
     axes=FALSE)
axis(2)
axis(1, at=seq(1:12), labels=row.names(summary), las=2)
mtext(text='Parameter', side=1, line=4.5)
box()
for(i in 1:12) lines(c(i,i),c(summary[i,4], summary[i,5]))
abline(h=0)


