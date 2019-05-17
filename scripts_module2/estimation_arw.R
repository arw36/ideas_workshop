### R code from vignette source 'estimation.rnw'
### Encoding: UTF-8
### Anna Willoughby
### May 16, 2019


load('data.RData')     #load the data and plot flu cases
plot(flu,type='b',log='y',main='Epidemic in a British boarding school', cex.main=0.85,
 xlab='Day', ylab='Active influenza cases')

# Exercise 1: 
ZNratio <- seq(0,1, by = 0.01)
plot(ZNratio, log(1-ZNratio)/-ZNratio)

model<-lm(log(flu[1:4])~day[1:4],data=flu);  #fit a linear model
summary(model)         #summary statistics for fit model
slope<-coef(model)[2]  #extract slope parameter
slope                 #print to screen

# Exercise 2 
# If time to infirmary cuts the infectious period, we need to adjust our value of gamma
# R0 will decreases 
# gamma = 1/1   # now 1 day instead of 2.5 days
R0 = 1.094913/1 + 1

# gamm = 1/2 # now 0.5 day, R0 decreases further
R0 = 1.094913/2 + 1


niamey[5,3]<-0  #replace a "NA"
niamey<-data.frame(biweek=rep(seq(1,16),3),site=c(rep(1,16),rep(2,16),rep(3,16)),
                   cases=c(niamey[,1],niamey[,2],niamey[,3])) #define "biweeks"
niamey1 <- filter(niamey, site ==1 )
niamey1$year <- niamey1$biweek/26 # convert biweek data to a year scale
plot(niamey1$year, niamey1$cases, type = 'b', xlab = 'Year', ylab = "Cases")
model_n <-lm(log(cases[1:8])~year[1:8],data=niamey1);  #fit a linear model
summary(model_n)         #summary statistics for fit model
slope_n<-coef(model_n)[2]  #extract slope parameter
slope_n     

gamma_n = (0.0384)^-1 # years 
# calculate R0 based on slope and gamma
R0 = slope_n/gamma_n + 1
R0 

# make a data frame to input values 
data <- data.frame(datapoints=seq(3, nrow(niamey1), by = 1))
get_SE <- function(x){
  model <-lm(log(cases[1:x])~year[1:x],data=niamey1);  #fit a linear model
  summary(model)         #summary statistics for fit model
  std_error <- coef(summary(model))[2, 2]
  return(std_error)
}

data %<>% mutate(SE = purrr::map(datapoints, get_SE))

gamma = (0.0384)^-1 
get_RO <- function(x){
  model <-lm(log(cases[1:x])~year[1:x],data=niamey1);  #fit a linear model
  summary(model)         #summary statistics for fit model
  slope <-coef(model)[[2]]
  RO = slope/gamma + 1
  return(RO)
}

data %<>% mutate(RO = purrr::map(datapoints, get_RO))
plot(data$RO, data$SE, type = 'b', xlab = 'RO', ylab = "Standard Error", col = datapoints)

load('data.RData')
niamey[5,3]<-0  #replace a "NA"
niamey<-data.frame(biweek=rep(seq(1,16),3),site=c(rep(1,16),rep(2,16),rep(3,16)),
                   cases=c(niamey[,1],niamey[,2],niamey[,3])) #define "biweeks"
plot(niamey$biweek,niamey$cases,type='p',col=niamey$site,xlab='Biweek',ylab='Cases')
lines(niamey$biweek[niamey$site==1],niamey$cases[niamey$site==1])
lines(niamey$biweek[niamey$site==2],niamey$cases[niamey$site==2],col=2)
lines(niamey$biweek[niamey$site==3],niamey$cases[niamey$site==3],col=3)

closed.sir.model <- function (t, x, params) {  #SIR model equations
  S <- x[1]
  I <- x[2]
  beta <- params
  dS <- -beta*S*I
  dI <- beta*S*I-(365/13)*I
  list(c(dS,dI))
}

sse.sir <- function(params0,data,site){  #function to calculate squared errors
  data<-data[data$site==site,]    #working dataset, based on site
  t <- data[,1]*14/365            #time in biweeks
  cases <- data[,3]               #number of cases
  beta <- exp(params0[1])            #parameter beta
  S0 <- exp(params0[2])           #initial susceptibles
  I0 <- exp(params0[3])           #initial infected        
  out <- as.data.frame(ode(c(S=S0,I=I0),times=t,closed.sir.model,beta,hmax=1/120))
  sse<-sum((out$I-cases)^2)       #sum of squared errors
}

library(deSolve)   #differential equation library
params0<-c(-3.2,7.3,-2.6)  #initial guess

fit1 <- optim(params0,sse.sir,data=niamey,site=1) #fit
exp(fit1$par)  #back-transform parameters
fit2 <- optim(params0,sse.sir,data=niamey,site=2) #fit
exp(fit2$par)  #back-transform parameters
fit3 <- optim(params0,sse.sir,data=niamey,site=3) #fit
exp(fit3$par)  #back-transform parameters

par(mfrow=c(3,1))   #set up plotting area for multiple panels
plot(cases~biweek,data=subset(niamey,site==1),type='b',col='blue', pch=21) #plot site 1
t <- subset(niamey,site==1)[,1]*14/365
mod.pred<-as.data.frame(ode(c(S=exp(fit1$par[2]),I=exp(fit1$par[3])),times=t,
                              closed.sir.model,exp(fit1$par[1]),hmax=1/120))
                              #obtain model predictions
lines(mod.pred$I~subset(niamey,site==1)[,1]) #and plot as a line

plot(cases~biweek,data=subset(niamey,site==2),type='b',col=site) #site 2
t <- subset(niamey,site==2)[,1]*14/365
mod.pred<-as.data.frame(ode(c(S=exp(fit2$par[2]),I=exp(fit2$par[3])),times=t,
                              closed.sir.model,exp(fit2$par[1]),hmax=1/120))
lines(mod.pred$I~subset(niamey,site==2)[,1])


plot(cases~biweek,data=subset(niamey,site==3),type='b',col=site) #site 3
t <- subset(niamey,site==3)[,1]*14/365
mod.pred<-as.data.frame(ode(c(S=exp(fit3$par[2]),I=exp(fit3$par[3])),times=t,
                              closed.sir.model,exp(fit3$par[1]),hmax=1/120))
lines(mod.pred$I~subset(niamey,site==3)[,1])

# Exercise 5: Estimate gamma and beta simultaneously 
# Exercise 6: What happens if one or both of the unknowns (I0 and S0) are fixed instead of gamma?
# Exercise 7: Modify your optimizer so that it returns the negative log-liklihood.What is it for the three districts of Niamey? 




