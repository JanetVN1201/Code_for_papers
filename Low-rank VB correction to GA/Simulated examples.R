library(tictoc)
library(mvtnorm)
###Example in section 3.2

x0 <- seq(-10,5,by = 0.001)

f1 <- function(x){
  a = aa
  ll = 0
  for (i in 1:N){
    ll = ll + dpois(a[i], lambda = exp(x[1] + x[2]*cov[i] + x[i+2] - mean(x)), log = T)
  }
  ll = ll + dt(x[2], df = 5, log = TRUE) + dunif(x[1], min = -3, max = 3, log = T) + dmvnorm(x[1:N+2], mean = rep(0, N), sigma = 0.25*diag(N), log = T)
  return(ll)
  }


mu = -1
beta1 = -0.5
N = 100
x_i = rnorm(N, mean = 0, sd = 0.5)
cov = rnorm(N)
aa = rpois(N, lambda  = exp(mu + beta1*cov + x_i))
par(mar = c(4,2,2,1))

barplot(table(as.factor(aa)), main = "", ylim = c(0, N))

tic("Gaussian")
if (T){
  tol = 1000
  mode1 = c(-1,-0.5, rep(0, N))
  while (tol>10){
    hess2 <- numDeriv::hessian(func=f1, x = mode1)
    grad2 <- numDeriv::grad(func=f1, x = mode1)
    mode2 = solve(-hess2, grad2-hess2%*%mode1)
    tol = sum(abs(mode1 - mode2))
    mode1 = mode2
  }
}
toc()

m = mode1

cov1 = solve(-hess2)

n_corr = 2

u_A = diag(N)
A = matrix(data = cbind(rep(1, N), cov, u_A),
                            ncol = 2+N)
cov_eta = A%*%cov1%*%t(A)

minf = function(delta) {
  #Loglik
  eta = A%*%(m + cov1[,1:n_corr]%*%delta)
  ll = 0
  for (i in 1:N){
    loglik = exp(eta[i] + cov_eta[i,i]/2) - 
      aa[i]*eta[i]
    ll = ll + loglik}
  
  #KLD
  xx = rbind(rep(seq(-3,  3, length = 20), each = 100),rep(seq( - 3,  3, length = 20), n = 100))
  area = (max(xx[1,]) - min(xx[1,]))*(max(xx[2,]) - min(xx[2,]))
  prior = rep(NA, ncol(xx))
  post = prior
  f = prior
  f2 = prior
  for (i in 1:ncol(xx)){
    prior[i] = dt(xx[2,i], df = 5)*dunif(xx[1,i], min = -3, max = 3)
    post[i] = dmvnorm(xx[,i], mean = m[1:n_corr] + cov1[1:n_corr,1:n_corr]%*%delta, sigma = cov1[1:n_corr,1:n_corr], log = F)
    f[i] = post[i]*(log(post[i]+0.001) - log(prior[i]+0.001))
  }
  kl = 2*area*sum(f)/ncol(xx)
  return (ll + kl)
}


tic("VBC")
d_est = optim(rep(0,n_corr), minf)
toc()
delta = d_est$par
improved_mean = m + cov1[,1:2]%*%d_est$par
improved_mean_eta = A%*%(m + (cov1[,1:2]%*%d_est$par))
delta_eta = A%*%(cov1[,1:2]%*%d_est$par)


#MCMC

model = "model {
  for(i in 1:N){
    y[i] ~ dpois(lambda[i])
    log(lambda[i]) = beta0 + beta1*x[i] + ustar[i]
}
  beta0 ~ dunif(-3,3)
  beta1 ~ dt(0,1,5)
  for (j in 1:N){
  u[j] ~ dnorm(0,4)
  ustar[j] <- u[j] - mean(u[])
  }
}"

model = "model {
  for(i in 1:N){
    y[i] ~ dpois(lambda[i])
    log(lambda[i]) = beta0 + beta1*x[i] + u[i]
}
  beta0 ~ dunif(-3,3)
  beta1 ~ dt(0,1,5)
  for (j in 1:N){
  u[j] ~ dnorm(0,4)
  }
}"

library(runjags)
tic("MCMC")
res.1 = run.jags(model = model,
                 monitor = c("beta0","beta1", "u[1]", "u[8]","u[15]"),
                 data = list('y' = aa, 'x' = cov, 'N' = N),
                 n.chains = 1,
                 inits = list(beta0 = 1, beta1 = -0.5, u = rep(0.5,N)),
                 burnin = 10^3,
                 sample = 10^5,
                 method = "parallel")

trace <- combine.mcmc(res.1)
toc()
mean.mcmc <- c(mean(trace[, "beta0"]), mean(trace[, "beta1"]), mean(trace[,"u[1]"]), mean(trace[,"u[8]"]), mean(trace[,"u[15]"]))
print(mean.mcmc)

####Stan using BRMS
library(rstan)
library(brms)
library(tidybayes)
library(dplyr)


tic("STAN")
sv <- stanvar(scode = "sum(z_1[1])~normal(0, N_1*0.001);",
              block = "model", position = "end")

mod_sim <- brm(y ~ 1 + xx + (1|group),
               family = poisson,
               data = data.frame(y = aa, xx = cov, group = 1:N),
               refresh = 0, silent = 2, cores = 4,
               prior = c(
                 prior(normal(0,1), class = Intercept),
                 prior(normal(0.5,0.001), class = sd)
               ),
               stanvars = sv,
               warmup = 10000, iter = 100000)
toc()

samplesss <- data.frame(mod_sim %>%
  spread_draws(r_group[group,]))
fixed_e <- data.frame(mod_sim %>%
                        spread_draws(c(b_Intercept,b_xx)))

stan_sam <- matrix(data = rep(NA, 360000*N/num_pg), ncol = N/num_pg, nrow = 360000)
for (i in 1:(N/num_pg)){
  stan_sam[,i] <- samplesss[samplesss$group == i, 2]
}

###########Example in Section 4.3
##MCMC - define
cat("model
    {

    for ( i in 1:N ) {
    u[i] ~ dnorm(0, tau)
    lambda[i] = exp(beta0 + beta1*x[i] + u[i] )
    y[i] ~ dpois( lambda[i])

    }

    ### Define the priors
    beta0 ~ dnorm( 0, 1 )
    beta1 ~ dnorm( 0, 1 )
    ltau ~ dgamma(1,1)
     tau = exp(ltau)

    }", file="sim_jags_VB.txt")
inits <- list(beta0 = 0, beta1 = 0, ltau = 0.5)
params <- c("beta0", "beta1", "tau")

###Data

if(TRUE){
  set.seed(0)
  N = 100
  b0 = -1
  u = rnorm(N, mean = 0, sd = sqrt(1/100))
  x = rnorm(N, mean = 0, sd = 1)
  b1 = -0.5
  eta = b0 + b1*x + u 
  Y = rpois(n = N, lambda = exp(eta))
  prec = 1
  
  par(mfrow = c(1,1))
  par(mar = c(4,2,2,1))
  barplot(table(as.factor(Y)), main = "", ylim = c(0, N))
  
  #Gaussian strategy
  INLA_result = inla(formula = Y ~ 1 + x + f(ID, model = "iid"), data=data.frame(Y = Y, x = x, ID = 1:N), 
                     family = "poisson", 
                     control.fixed = list(prec = prec,   
                                          prec.intercept = prec,
                                          correlation.matrix=TRUE),
                     control.compute=list(config = TRUE),
                     control.inla = list(strategy = "gaussian",
                                         control.vb = list(enable = FALSE)))
  
  #Laplace strategy
  INLA_result_L = inla(formula = Y ~ 1 + x + f(ID, model = "iid"), data=data.frame(Y = Y, x = x, ID = 1:N), 
                       family = "poisson", 
                       control.compute=list(config = TRUE),
                       control.fixed = list(prec = prec,   
                                            prec.intercept = prec,
                                            correlation.matrix=TRUE), 
                       control.inla = list(strategy = "laplace"))
  
  
  #INLA-VBC
  INLA_result_VB = inla(formula = Y ~ 1 + x + f(ID, model = "iid"), data=data.frame(Y = Y, x = x, ID = 1:N), 
                        family = "poisson", 
                        control.compute=list(config = TRUE),
                        control.fixed = list(prec = prec,   
                                             prec.intercept = prec,
                                             correlation.matrix=TRUE), 
                        control.inla = list(strategy = "gaussian",
                                            control.vb = list(enable = TRUE)))
  
  #MCMC
  library(rjags)
  dat1 = list(y = Y, x = x, N = N)
  jags.m <- jags.model( file = "sim_jags_VB.txt", data=dat1, inits=inits, n.chains=1, n.adapt=500 )
  samps <- coda.samples( jags.m, params, n.iter=10^5, start = (10^2+1) )
  MCtime = system.time(coda.samples( jags.m, params, n.iter=10^5, start = (10^2+1) ))
  
}


