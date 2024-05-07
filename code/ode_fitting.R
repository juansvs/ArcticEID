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
  I0 = data_subset[[which_incidence]][1]/(params$beta)
  IC = c(S = S0, I = I0, R = 0)
  # The times at which we compute the solution to compare with data
  times = as.numeric(as.Date(incidence_data$date))
  sol = ode(IC, times, RHS_KMK_SIR, params, method = method)
  # Error checking
  if (sol[dim(sol)[1],"time"] < times[length(times)]) {
    return(Inf)
  }
  # Values required to compute the error
  incidence_from_run = sol[,"I"]/(sol[,'S']+sol[,"I"]+sol[,"R"])
  # Compute the error
  diff_values = incidence_data[[which_incidence]] - incidence_from_run
  diff_values_squared = diff_values^2
  error = sum(diff_values_squared)
  return(error)
}

# Read in the data
data = readr::read_csv("../data/caribou_seroprevalence.csv")
# Show the top of the table (always useful)
head(data)
# It will be easier to work with days as days, so we convert dates
data$date = lubridate::parse_date_time(data$Month,'%B')

# Ah, yes, don't forget: let's flip the order of the table, as it is sorted
# by decreasing dates
data = data %>%  arrange(date)

# Type of incidence to use (column name in the data)
which_incidence = "Prevalence"

data_subset = filter(data, Total>8)

ggplot(data_subset,aes(format(date,"%m"),Prevalence))+geom_point(aes(color=Herd))


params = list()
params$gamma = 1/10   # Let's see if we can fit with this simple value
# R0=beta/gamma, so beta=R0*gamma
params$beta = 1
params$pop = 0.5 # density of caribou

## Fit using a genetic algorithm
GA = ga(
  type = "real-valued",
  fitness = function(x) -error_incidence(p_vary = c(beta = x[1], gamma = x[2]),
                                         params = params,
                                         incidence_data = data_subset,
                                         which_incidence = which_incidence,
                                         method = "rk4"),
  parallel = 4,
  lower = c(1, 0.01),
  upper = c(3, 0.1),
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
I0 = data_subset[[which_incidence]][1]/(params$beta)
IC = c(S = S0, I = I0, R = 0)
dates_num = as.numeric(as.Date(data_subset$date))
times = seq(dates_num[1], dates_num[length(dates_num)], 0.1)
sol <- ode(IC, times, RHS_KMK_SIR, params)
sol_incidence = as.numeric(params$beta) * sol[,"S"] * sol[,"I"]/params$pop
y_max = max(max(sol_incidence), max(data_subset[[which_incidence]]))

plot(sol[,"time"], sol_incidence, type = "l",
     xlab = "Date", ylab = "Incidence",
     xaxt = "n", lwd = 2, col = "red",
     ylim = c(0, y_max))
lines(as.numeric(as.Date(data_subset$date)), data_subset[[which_incidence]],
      type = "b")
dates_pretty = pretty(dates_num)
axis(1, at = dates_pretty, labels = as.Date(dates_pretty, origin="1970-01-01"))
