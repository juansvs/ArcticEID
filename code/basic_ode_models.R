library(deSolve)

# Basic SIR model: spatially implicit, no group structure, continuous time

# Parameters. Rates are per year
pars <- c(
  gamma_H = 1.1, # fecundity
  gamma_I = 0.2,
  nu = 1, # loss of immunity in one year
  beta = 20, # transmission. wild guess at this point
  mu_N = 1/12, # based on 15 year lifespan (12-20 yr)
  rho = 365/5, # # let's assume you either die or recover quickly, in five days, and that 80% of infected die
  phi = 0.8, # probability of death
  delta = 1/2, # a carcass can last for 2 years? (doesn't matter for this model)
  K = 0.4 # carrying capacity
)

# State variables, densities of individuals per km^-2
state <- c(
  S = 0.3,
  I = 0.01,
  D = 0,
  R = 0
)

# ODE system as function
ode1 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    Nt <- S+I+R
    dS <- gamma_H*(Nt)*(1-Nt/K)+nu*R-beta*S*(I+D)-mu_N*S
    dI <- beta*S*(I+D)-rho*I-mu_N*I
    dD <- phi*rho*I-delta*D
    dR <- (1-phi)*rho*I-nu*R-mu_N*R
    
    list(c(dS,dI,dD,dR))
  })
}

# get the solution up to 20 years, with weekly resolution
times <- seq(0,20,by = 1/52)

# solve numerically
out <- ode(y = state, times = times, func = ode1, parms = pars)

# plot
matplot.deSolve(out, 
                main="", ylab = expression(paste("Density (",km^-2,")")), xlab = "Time (years)", lty = 1, 
                col = hcl.colors(4),cex.lab=1.5)

                