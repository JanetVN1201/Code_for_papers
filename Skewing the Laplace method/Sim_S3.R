#Section 3 examples
library(rjags)
library(INLA)

##Poisson - Ex 3.3.1
set.seed(101)
if (T){
n <- 300   # number of observations  
ID <- 1:n

alpha_true <- -1.5
tau_true <- 1
u_true <- rnorm(n, 0, 1/sqrt(tau_true))
lambda <- exp(alpha_true + u_true)
y <- rpois(n, lambda)
hist(y)
# JAGS model
model_string <- "
model {
  # Likelihood
  for (i in 1:N) {
    y[i] ~ dpois(lambda[i])
    log(lambda[i]) <- alpha + u[i]
  }

  # Random effects
  for (j in 1:N) {
    u[j] ~ dnorm(0, 1)
  }

  # Priors
  alpha ~ dnorm(0, 0.001)
}
"

# Prepare data for JAGS
data_jags <- list(
  y = y,
  N = n
)

# Initial values
inits <- function() {
  list(alpha = 0)
}

# Parameters to monitor
params <- c("alpha")

# Run model
model <- jags.model(textConnection(model_string),
                    data = data_jags,
                    inits = inits,
                    n.chains = 3,
                    quiet = TRUE)

update(model, 10000) # burn-in

samples <- coda.samples(model,
                        variable.names = params,
                        n.iter = 40000)

#INLA + VB
res_INLA_M <- inla(y ~ 1 + f(ID, model = "iid", hyper = list(prec = list(initial = 0, fixed = T))),
                   data = data.frame(y = y,
                                     ID = ID),
                   family = "poisson")
res_INLA_MV <- inla(y ~ 1 + f(ID, model = "iid", hyper = list(prec = list(initial = 0, fixed = T))),
                   data = data.frame(y = y,
                                     ID = ID),
                   family = "poisson",
                   control.inla = list(control.vb = list(strategy = "variance")))

#Plot
par(mfrow = c(1,1))
hist(data.frame(samples[[1]])[,1], n = 50, prob = T)
lines(res_INLA_M$marginals.fixed$`(Intercept`[,1],
res_INLA_M$marginals.fixed$`(Intercept`[,2], col = "red")
lines(res_INLA_MV$marginals.fixed$`(Intercept`[,1],
      res_INLA_MV$marginals.fixed$`(Intercept`[,2], col = "blue", lty = 2)
}

##Student t - Ex 3.3.2
set.seed(111)
if (T){
n <- 10
x <- rnorm(n)
alpha_true <- 0
beta_true <- 1
tau_true <- 1
nu_true <- 4

# Generate Student-t distributed data
y <- alpha_true + beta_true * x + 
  rt(n, df = nu_true)/sqrt(nu_true/(nu_true - 2))

# JAGS model string
model_string <- "
model {
  for (i in 1:N) {
    y[i] ~ dt(mu[i], 2, 4)   # Student-t
    mu[i] <- alpha + beta * x[i]
  }
  
  # Priors
  alpha ~ dnorm(0, 0.001)
  beta  ~ dnorm(0, 0.001)
}
"

# Data list for JAGS
data_jags <- list(
  y = y,
  x = x,
  N = n
)

# Initial values
inits <- function() {
  list(alpha = 0, beta = 0)
}

# Parameters to monitor
params <- c("alpha", "beta")

# Run JAGS
model <- jags.model(textConnection(model_string),
                    data = data_jags,
                    inits = inits,
                    n.chains = 3,
                    quiet = TRUE)

update(model, 10000) # burn-in

samples <- coda.samples(model,
                        variable.names = params,
                        n.iter = 90000)

#INLA + VB
res_INLA_M <- inla(y ~ 1 + x,
                   data = data.frame(y = y,
                                     x = x),
                   family = "T",
                   control.family = list(hyper = list(theta1 = list(initial = 0, fixed = T), theta2 = list(initial = log(2), fixed = T))))
res_INLA_MV <- inla(y ~ 1 + x,
                    data = data.frame(y = y,
                                      x = x),
                    family = "T",
                    control.inla = list(control.vb = list(strategy = "variance")),
                    control.family = list(hyper = list(theta1 = list(initial = 0, fixed = T), theta2 = list(initial = log(2), fixed = T))))

                                 

par(mfrow = c(1,1))
hist(data.frame(samples[[1]])[,1], n = 50, prob = T)
lines(res_INLA_M$marginals.fixed$`(Intercept`[,1],
      res_INLA_M$marginals.fixed$`(Intercept`[,2], col = "red")
lines(res_INLA_MV$marginals.fixed$`(Intercept`[,1],
      res_INLA_MV$marginals.fixed$`(Intercept`[,2], col = "blue", lty = 2)

hist(data.frame(samples[[1]])[,2], n = 50, prob = T)
lines(res_INLA_M$marginals.fixed$x[,1],
      res_INLA_M$marginals.fixed$x[,2], col = "red")
lines(res_INLA_MV$marginals.fixed$x[,1],
      res_INLA_MV$marginals.fixed$x[,2], col = "blue", lty = 2)

print(sd(data.frame(samples[[1]])[,1]))
print(sd(data.frame(samples[[1]])[,2]))
print(res_INLA_M$summary.fixed$sd)
print(res_INLA_MV$summary.fixed$sd)
}

##G Pareto - Ex 3.3.3
library(Rcpp)
sourceCpp("~/Downloads/GPareto.cpp")

##Spec-Sens - Ex 3.3.4
set.seed(123)
if (T){
library(statmod)
library(numDeriv)


n_sample <- 50
inverse_link <- function(eta) return(1/(1 + exp(-eta)))

p0 <- 0.8
p1 <- 0.985

pdf <- function(eta) return(p0*inverse_link(eta) + (1-p1)*(1-inverse_link(eta)))


eta_i <- rnorm(n_sample, mean = 0, sd = 1000)
binomial_probs <- pdf(eta_i)
y_i <- rbinom(n_sample, size = 1, prob = binomial_probs)


log_likelihood <- function(y, eta){
  
  probs <-  pdf(eta)
  to_return <- sum(y*log(probs) + (1-y) * log(1 - probs))
  return(to_return) 
}


# Plot a
eta_to_plot <- seq(-1.5, 3, by = 0.01)
log_likelihood_to_plot <- unlist(lapply(eta_to_plot, log_likelihood, y = y_i))
likelihood_to_plot <- exp(log_likelihood_to_plot)

# Plot b 
eta_to_plot1 <- seq(-10, 10, by = 0.01)
log_likelihood_to_plot1 <- unlist(lapply(eta_to_plot1, log_likelihood, y = y_i))

## Gaussian approximation
log_likelihood_fn <- splinefun(eta_to_plot1, log_likelihood_to_plot1)
find_mode <- optim(par = 0,
                   log_likelihood_fn,
                   method = "BFGS",
                   control = list(fnscale = -1),
                   hessian = TRUE)
mode_like <- find_mode$par
std_dev_like <- sqrt(-1/drop(hessian(log_likelihood_fn, x = find_mode$par)))

gauss_approx_log_like_to_plot <- dnorm(x = seq(-5, 5, 0.01), mean = mode_like, sd = std_dev_like, log = T)

# Plot c
second_deriv <- function(eta) return(hessian(log_likelihood_fn, eta))
second_deriv_vals <- unlist(lapply(eta_to_plot1, second_deriv))

# Plot d


## True Posterior
log_prior_vals <- dnorm(eta_to_plot, mean = 0, sd = 1000, log = T)
log_likelihood_vals <- log_likelihood_to_plot
log_posterior <- log_likelihood_vals + log_prior_vals
true_posterior <- exp(log_posterior)/sum(exp(log_posterior)*0.01)

## Gaussian Approximation
log_true_posterior_fn <- splinefun(eta_to_plot, log(true_posterior))
true_posterior_optim <- optim(par = 0,
                              log_true_posterior_fn,
                              method = "BFGS",
                              control = list(fnscale = -1),
                              hessian = TRUE)

gauss_approx_mean <- true_posterior_optim$par
gauss_approx_sd <- sqrt(-1/true_posterior_optim$hessian)

## VB - Variational Bayes

kl_normal <- function(mu1, sigma1, mu2, sigma2) {
  term1 <- log(sigma2 / sigma1)
  term2 <- (sigma1^2 + (mu1 - mu2)^2) / (2 * sigma2^2)
  kl <- term1 + term2 - 0.5
  return(kl)
}


vb_var <- function(delta){
  
  out <- gauss.quad.prob(n = 21,
                         dist = "normal",
                         mu = gauss_approx_mean,
                         sigma = drop(gauss_approx_sd)*exp(delta))
  
  exp_negative_log_like <- unlist(lapply(out$nodes, log_likelihood, y = y_i))
  exp_negative_log_like <- -sum(exp_negative_log_like*out$weights)
  
  KLD <- drop(kl_normal(mu1 = gauss_approx_mean,
                        sigma1 = gauss_approx_sd*exp(delta),
                        mu2 = 0,
                        sigma2 = 1000))
  
  return(exp_negative_log_like + KLD)
}

vb_optimise <- optim(par = 0, 
                     vb_var,
                     method = "BFGS",
                     hessian = TRUE)


# Plot
par(mfrow = c(2,2))
# Plot a
plot(eta_to_plot, likelihood_to_plot, type = "l")
# Plot b
plot(eta_to_plot1, log_likelihood_to_plot1, type = "l")
lines(seq(-5, 5, 0.01), gauss_approx_log_like_to_plot + find_mode$val, col = "red", lty = 2, lwd = 2)
# Plot c
plot(eta_to_plot1, second_deriv_vals, type = "l")
abline(h = -1/std_dev_like^2, col = "blue", lwd = 2, lty = 2)
abline(h = 0, col = "red", lwd = 2)
# Plot d
plot(eta_to_plot, true_posterior, type =  "l", col = "green", lwd = 1, lty = 3)
lines(eta_to_plot, dnorm(eta_to_plot, mean = gauss_approx_mean, sd = gauss_approx_sd), col = "red", lwd =1)
lines(eta_to_plot, dnorm(eta_to_plot, mean = gauss_approx_mean, sd = gauss_approx_sd*exp(vb_optimise$par)), col = "blue", lty = 2, lwd = 1)

#Values
print(gauss_approx_sd)
print(gauss_approx_sd*exp(vb_optimise$par))
eta_mean = sum(eta_to_plot*true_posterior*0.01)
print(sqrt(sum((eta_to_plot - eta_mean)^2*true_posterior*0.01)))
}