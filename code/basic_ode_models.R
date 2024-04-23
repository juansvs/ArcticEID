library(deSolve)

# Basic SIR model: spatially implicit, no group structure, continuous time

# Parameters. Rates are per year
pars <- c(
  phi = 0.5, # fecundity of healthy individuals
  sigma = 1, # proportion of fecundity of disease individuals, relative to healthy ones
  lambda = 1, # loss of immunity in one year
  beta = 20, # transmission. wild guess at this point. Units are km^2/ind/yr
  mu_N = 1/12, # based on 15 year lifespan (12-20 yr)
  mu_I = 365/5, # disease-induced mortality rate
  theta = 0.8, # proportion of individuals that develop a pathology and die
  rho = 365/15, # let's assume you either die or recover quickly, in five days, and that 80% of infected die
  delta = 1/2 # a carcass can last for 2 years? (doesn't matter for this model)
)

# State variables, densities of individuals per km^-2
state <- c(
  S = 0.3,
  I = 0.01,
  R = 0,
  D = 0
)

# ODE system as function
ode1 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- phi*(S+R)+sigma*phi*I-mu_N*S-beta*S*(I+D)+lambda*R
    dI <- -mu_N*I+beta*S*(I+D)-theta*mu_I*I-(1-theta)*rho*I
    dR <- -mu_N*R+(1-theta)*rho*I-lambda*R
    dD <- theta*mu_I*I-delta*D
    
    list(c(dS,dI,dR,dD))
  })
}

odeMI <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- phi*(S)+sigma*phi*I-mu_N*S-beta*S*(I+D)+lambda*R
    dI <- -mu_N*I+beta*S*(I+D)-theta*mu_I*I-(1-theta)*rho*I
    dR <- phi*R-mu_N*R+(1-theta)*rho*I-lambda*R
    dD <- theta*mu_I*I-delta*D
    
    list(c(dS,dI,dR,dD))
  })
}

odeSSMI <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSj <- phi*(Sa+Ia)+sigma*phi*Ia-gama*Sj-mu_J*Sj-beta*Sj*(Ij+Ia+D)+lambda*Rj
    dIj <- -mu_J*Ij-gama*Ij+beta*Sj*(Ij+Ia+D)-theta*mu_I*Ij-(1-theta)*rho*Ij
    dRj <- phi*Ra-gama*Rj-mu_J*Rj+(1-theta)*rho*Ij-lambda*Rj
    dSa <- gama*Sj-mu_N*Sa-beta*Sa*(Ij+Ia+D)+lambda*Ra
    dIa <- gama*Ij-mu_N*Ia+beta*Sa*(Ij+Ia+D)-theta*mu_I*Ia-(1-theta)*rho*Ia
    dRa <- gama*Rj-mu_N*Ra+(1-theta)*rho*Ia-lambda*Ra
    dD <- theta*mu_I*(Ij+Ia)-delta*D
    
    list(c(dSj,dIj,dRj,dSa,dIa,dRa,dD))
  })
}
                