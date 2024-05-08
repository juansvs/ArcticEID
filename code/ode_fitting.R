#' Code taken from J. Arino's github, for fitting the solution to a 
#' ODE of an epidemiological model to actual data. 

library(lubridate)
library(tidyverse)
library(dplyr)
library(wbstats)
library(deSolve)
library(GA)
library(parallel)

# The RHS function for KMK
RHS_KMK_SIR <- function(t, x, p) {
  with(as.list(c(x,p)), {
    dS <- -beta*S*I
    dI <- beta*S*I-gamma*I
    dR <- gamma*I
    return(list(c(dS, dI, dR)))
  })
}

# Compute the error based on one parameter value. Will also need a few more
# arguments
error_incidence <- function(p_vary, 
                            params, 
                            incidence_data,
                            method = "rk4",
                            which_incidence = "Prevalence") {
  # Anything that changes during optimisation needs to be set here
  params$beta = as.numeric(p_vary["beta"])
  params$gamma = as.numeric(p_vary["gamma"])
  # Check value of R0. If it is less than 1, no need to compute ODE solution,
  # just return Inf
  R0 = params$beta / params$gamma
  if (R0<1) {
    return(Inf)
  }
  # I0 is the I that gives us the incidence we are matching in the data. Since
  # incidence=beta*S*I/N, we have
  # I(0) = incidence*N/(beta*S(0)) ~= incidence/beta
  S0 = params$pop
  I0 = S0*data_subset[[which_incidence]][1]
  IC = c(S = S0, I = I0, R = 0)
  # The times at which we compute the solution to compare with data
  times = incidence_data$monthnum
  sol = ode(IC, times, RHS_KMK_SIR, params, method = method)
  # Error checking
  if (sol[dim(sol)[1],"time"] < times[length(times)]) {
    return(Inf)
  }
  # Compute the error
  matchdb = sol[match(incidence_data$monthnum, sol[,"time"]),]
  # calculate the difference between observed and model predicted. This
  # difference is weighted by the number of samples
  diff_values =   matchdb[,"I"]/(matchdb[,"S"]+matchdb[,"I"]+matchdb[,"R"])-incidence_data[,which_incidence]
  diff_values_squared = diff_values^2
  error = sum(diff_values_squared)
  return(error)
}

# Read in the data
data = readr::read_csv("../data/caribou_seroprevalence.csv")
# Show the top of the table (always useful)
head(data)

# Type of incidence to use (column name in the data)
which_incidence = "Prevalence"

data_subset = group_by(data, Herd, monthnum) %>% filter(sum(Total)>8) %>% 
  ungroup(Herd) %>% summarise(pos = sum(Positives), Total = sum(Total), Prevalence = pos/Total)
  

ggplot(data,aes(monthnum,Prevalence))+geom_point(aes(color=Herd, size=Total))


params = list()
params$gamma = 0.3   # Let's see if we can fit with this simple value
# R0=beta/gamma, so beta=R0*gamma
params$beta = 3
params$pop = 0.2 # density of caribou

## Fit using a genetic algorithm
GA = ga(
  type = "real-valued",
  fitness = function(x) -error_incidence(p_vary = c(beta = x[1], gamma = x[2]),
                                         params = params,
                                         incidence_data = data_subset,
                                         which_incidence = "Prevalence",
                                         method = "rk4"),
  parallel = 4,
  lower = c(0.5, 0.05),
  upper = c(5, 2),
  pcrossover = 0.7,
  pmutation = 0.2,
  optim = TRUE,
  optimArgs = list(method = "CG"),
  suggestions = c(params$beta, params$gamma),
  popSize = 500,
  maxiter = 100
)

params$beta = GA@solution[1]
params$gamma = GA@solution[2]
S0 = 0.2
I0 = S0*data_subset[[which_incidence]][1]
IC = c(S = S0, I = I0, R = 0)
dates_num = data_subset$monthnum
times = seq(dates_num[1], dates_num[length(dates_num)], 0.1)
sol <- ode(IC, times, RHS_KMK_SIR, params)
sol_incidence =  sol[,"I"]/(sol[,"S"]+sol[,"I"]+sol[,"R"])
y_max = max(max(sol_incidence), max(data_subset[[which_incidence]]))

plot(sol[,"time"], sol_incidence, type = "l",
     xlab = "Month", ylab = "Prevalence", lwd = 2, col = "red",
     ylim = c(0, y_max))
points(data_subset$monthnum, data_subset[[which_incidence]])

#' When we aggregate the herds and years together, this produces beta
#' estimates of 2.99/individual/month, and recovery rates, 0.226/month.
#' This translates to beta = 35.9 /ind/y and gamma = 2.72 /y