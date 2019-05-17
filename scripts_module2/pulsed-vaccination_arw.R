### R code from vignette source 'pulsed-vaccination.rnw'

require(deSolve)                          #deSolve library needed for this computing session

sir.model.open <- function (t, x, params) {    #here we begin a function with three arguments
  S <- x[1]                               #create local variable S, the first element of x
  I <- x[2]                               #create local variable I
  R <- x[3]                               #create local variable R
  with(                                   #we can simplify code using "with"
       as.list(params),                   #this argument to "with" lets us use the variable names
       {                                  #the system of rate equations
         dS <- mu*(S+I+R) - beta*S*I - mu*S
         dI <- beta*S*I - gamma*I - mu*I
         dR <- gamma*I - mu*R
         dx <- c(dS,dI,dR)                #combine results into a single vector dx
         list(dx)                         #return result as a list
       }
       )
}

R0 <- 10
N <-  1					                       #population size
mu <- 0.02                             #per capita birth/death rate
gamma <- 365/10  			                 #recovery rate (in years)
beta <- R0*(gamma+mu)/N 	             #transmission rate
xstart <- c(S=0.2, I=0.001, R=1-0.2-0.001)	 #initial conditions, must sum to one
Tmax <- 120                            #integrate for 200 years after transients
params <- c(beta=beta, gamma=gamma, mu=mu)  #parameter vector
tau <- 0.1                             #size of time step
times <- seq(0, Tmax, by=tau)          #function seq returns a sequence

out <- ode(xstart,times,sir.model.open,params, method='ode45', rtol=1e-7)  

op <- par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))                      #set graphical parameters
plot(I~time,data=out, type='l', lwd=2)                          #plot the I variable against time
par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)                      #re-set graphical parameters
plot(I~S,data=out,log='xy',yaxt='n',xlab='S', type='l', lwd=2)  #plot phase portrait
par(op)                                                         #re-set graphical parameters

xstart <- out[which(out[,1]==50),2:4] # start pulsed vaccination strategies once initial system reached equilibrium

pv <- 0.1                     # fraction of susceptibles vaccinated
Tv <- 4                       # number of years between pulses
vacc.events <- floor(Tmax/Tv) # number of pulses in Tmax years

data <- data.frame(S=out[which(out[,1]==50),2],
                   I=out[which(out[,1]==50),3],
                   R=out[which(out[,1]==50),4])

for(i in 1:vacc.events){
  out <- ode(xstart, seq(tau, Tv, by=tau), sir.model.open, params, method='ode45', rtol=1e-7)
  xstart <- out[dim(out)[1],2:4]        # reset initial condition
  xstart[1] <- (1-pv)*(tail(out,1)[2])  # vaccinate susceptibles
  xstart[3] <- xstart[3]+(pv)*(tail(out,1)[2])  # move to recovered class
  data <- rbind(data,out[,2:4])         # store result
}

data$time <- seq(50, Tmax+50, by=tau)
par(mar=c(5,4,4,4)+0.1)
plot(data$time[1:500], data$I[1:500], type='l', xlab='Time', ylab='', col='red', axes=FALSE)
axis(2, col.axis='red')
mtext(side=2, line=2.5, 'Infected', col='red')
box()
axis(1)
par(new=TRUE)
plot(data$time[1:500], data$S[1:500], type='l', xlab='', ylab='', axes=FALSE, col='black')
axis(4)
mtext(side=4, line=2.5, 'Susceptibles')


# Exercise 1. Why might this periodicity be of interest? 
#             I think periodicity is of interest because this allows strategic control methods. 
#             What might be the consequences of “peaks” and “troughs” for public health?
#             This could lead to targetted vaccination times, or the ordering/stockpiling schedules of certain treatments. 
# Exercise 2. By modifying the value of pv, see if you can locate the vaccination threshold. 
#             
pv.test <- 0.65 # varied this value from 0.1 - 0.8 in 0.2 increments. At 65% coverage, there is an initial infected but now epidemic

data.test <- data.frame(S=out[which(out[,1]==50),2],
                   I=out[which(out[,1]==50),3],
                   R=out[which(out[,1]==50),4])

for(i in 1:vacc.events){
  out <- ode(xstart, seq(tau, Tv, by=tau), sir.model.open, params, method='ode45', rtol=1e-7)
  xstart <- out[dim(out)[1],2:4]        # reset initial condition
  xstart[1] <- (1-pv.test)*(tail(out,1)[2])  # vaccinate susceptibles
  xstart[3] <- xstart[3]+(pv.test)*(tail(out,1)[2])  # move to recovered class
  data.test <- rbind(data.test,out[,2:4])         # store result
}

time <- seq(50, Tmax+50, by=tau) # idk why getting this error, appending time from data 
time <- time[time != "50"]
data.test$time <- time
par(mar=c(5,4,4,4)+0.1)
plot(data.test$time[1:500], data.test$I[1:500], type='l', xlab='Time', ylab='', col='red', axes=FALSE)
axis(2, col.axis='red')
mtext(side=2, line=2.5, 'Infected', col='red')
box()
axis(1)
par(new=TRUE)
plot(data.test$time[1:500], data.test$S[1:500], type='l', xlab='', ylab='', axes=FALSE, col='black')
axis(4)
mtext(side=4, line=2.5, 'Susceptibles')

# Does this agree with the analytically predicted threshold?
# Based on the analytically predicted threshold for an R0 of 10, 
# ((mu*T - pv)* (e^mu*T - 1) + mu*pv*T)/(mu*T(pv - 1 + e^mu*T)) < 1/R0 
threshold <- function(mu, Tv, pv, R0){
  ((mu*Tv - pv)* (exp(1)^(mu*Tv) - 1) + mu*pv*Tv)/(mu*Tv*(pv - 1 + exp(1)^(mu*Tv))) < 1/R0 
}

pv.values <- seq(0.01, 1, by = 0.01)
table(pv.values, threshold(0.02, 4, pv.values, 10))
# In fact, this calculation shows that the threshold is actually slightly lower, 0.55

# Exercise 3. One thing we might be interested in knowing is how our vaccination strategy aﬀects mean disease prevalance. 
#             Devise a strategy to study this question and produce a plot illustrate the result. 

# To study disease prevalence, we will calculate the number of mean I across the time period with varying vaccination prevalance (pv) and vaccination events (Tv)
# First we can calculate the mean I for the originally given conditions. 
mean(data$I)

# It seems reasonable to vary pv at 5% levels
pv.values <- seq(0.05, 1, by = 0.05)
# It seems reasonable to vary vaccination events from 1:50. 
Tv.values <- seq(1, 50, by = 1)

# cross these two vectors to determine all possible vaccination strategies
vacc.strategies <- tidyr::crossing(pv.values, Tv.values)

# Make a function to calculate mean I based on pv and Tv
calc.meanI <- function(pv, Tv){
  vacc.events <- floor(Tmax/Tv)
  data = data.frame(S=out[which(out[,1]==50),2],
                                    I=out[which(out[,1]==50),3],
                                    R=out[which(out[,1]==50),4])
for(i in 1:vacc.events){
  out <- ode(xstart, seq(tau, Tv, by=tau), sir.model.open, params, method='ode45', rtol=1e-7)
  xstart <- out[dim(out)[1],2:4]        # reset initial condition
  xstart[1] <- (1-pv)*(tail(out,1)[2])  # vaccinate susceptibles
  xstart[3] <- xstart[3]+(pv)*(tail(out,1)[2])  # move to recovered class
  data <- rbind(data,out[,2:4])         # store result
}
  mI <- mean(data$I)
  return(mI)
} 

# Run this function on all the possible vaccination parameters 
library(magrittr)
vacc.strategies %<>% mutate(mI = purrr::map2(pv.values, Tv.values, calc.meanI)) # takes about 3 minutes 
vacc.strategies$meanI <- unlist(vacc.strategies$mI)
# plot the mean infection prevalence as a heatmap based on pv and Tv
ggplot(vacc.strategies, aes(x = pv.values, y = Tv.values, fill = meanI)) + 
  geom_tile() +
  labs(x = "Vaccination prevalence", y = "Pulse vaccination frequency (years)", fill = "Infection prevalence (mean)") + 
  scale_fill_viridis_c(direction = -1, labels = function(x){sprintf("%.4f", x)})

 

