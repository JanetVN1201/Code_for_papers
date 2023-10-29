library(tictoc)
library(mvtnorm)
###Example in section 3.2

x0 <- seq(-10,5,by = 0.001)

f1 <- function(x){
  a = aa
  ll = 0
  for (i in 1:N){
    ll = ll + dpois(a[i], lambda = exp(x[1]*cov[i] + x[i+1]), log = T) + dnorm(x[i+1], mean = 0, sd = 0.5, log = T)
  }
  ll = ll + dt(x[1], df = 5, log = TRUE)
  return(ll)
}

mu = -1
beta1 = -0.5
N = 20
x_i = rnorm(N, mean = 0, sd = 0.5)
cov = rnorm(N)
aa = rpois(N, lambda  = exp(beta1*cov + x_i))
par(mar = c(4,2,2,1))

barplot(table(as.factor(aa)), main = "", ylim = c(0, N))

tic("Gaussian")
if (T){
 tol = 2
 mode1 = c(1, rep(0.1, N))
 while (tol>0.001){
 hess2 <- numDeriv::hessian(func=f1, x = mode1)
 grad2 <- numDeriv::grad(func=f1, x = mode1)
 mode2 = solve(-hess2, grad2-hess2%*%mode1)
 tol = sum(abs(mode1 - mode2))
 mode1 = mode2
 }
}

hess2 <- numDeriv::hessian(func=f1, x = mode1)
m = mode1
toc()
cov1 = solve(-hess2)

minf = function(delta) {
  #Loglik
  A = matrix(data = cbind(cov, diag(N)), ncol = N+1)
  eta = matrix(data = cbind(cov, diag(N)), ncol = N+1)%*%(m + cov1[,1]*delta)
  cov_eta = A%*%cov1%*%t(A)
  
  ll = 0
  for (i in 1:N){
  loglik = exp(eta[i] + cov_eta[i,i]/2) - 
    aa[i]*eta[i]
  ll = ll + loglik}
  
  #KLD
  xx = seq(- 3,  3, length = 20)
  area = 6
  prior = rep(NA, length(xx))
  post = prior
  f = prior
  f2 = prior
  for (i in 1:length(xx)){
  prior[i] = dt(xx[i], df = 5)#dnorm(xx[1,i], mean = 0, sd = 1, log = FALSE)#dnorm(xx[1,i], mean = 0, sd = 1, log = F)*dnorm(xx[2,i], mean = 0, sd = 1, log = F)
  post[i] = dnorm(xx[i], mean = m[1] + cov1[1,1]*delta, sd = sqrt(cov1[1,1]), log = F)
  f[i] = post[i]*(log(post[i]+0.001) - log(prior[i]+0.001))
  }
  kl = area*sum(f)/length(xx)
  return (ll + kl)
}


tic("VBC")
d_est = optim(c(-1), minf, method = "Brent", lower = -25, upper = 25)
toc()
delta = d_est$par
improved_mean = m + cov1[,1]*d_est$par
improved_mean_eta = cbind(cov, diag(N))%*%m + cbind(cov, diag(N))%*%(cov1[,1]*d_est$par)
delta_eta = cbind(cov, diag(N))%*%(cov1[,1]*d_est$par)

par(mar = c(2,2,1,1))
plot(c(cov1[,1]*d_est$par, delta_eta), pch = c(rep(2,N+1),rep(16, N)), cex = c(rep(1,N+1), rep(0.5, N)), xaxt = "n", xlab = "", ylab = expression(paste(delta)),
     ylim = c(-0.2, 0.2))
abline(h=0, lty = 2)
axis(1, at=1:(N+N+1), labels=c(expression(paste(beta[1])),
                             expression(paste(u[1])),
                             expression(paste(u[2])),
                             expression(paste(u[3])),
                             expression(paste(u[4])),
                             expression(paste(u[5])),
                             expression(paste(u[6])),
                             expression(paste(u[7])),
                             expression(paste(u[8])),
                             expression(paste(u[9])),
                             expression(paste(u[10])),
                             expression(paste(u[11])),
                             expression(paste(u[12])),
                             expression(paste(u[13])),
                             expression(paste(u[14])),
                             expression(paste(u[15])),
                             expression(paste(u[16])),
                             expression(paste(u[17])),
                             expression(paste(u[18])),
                             expression(paste(u[19])),
                             expression(paste(u[20])),
                             expression(paste(eta[1])),
                             expression(paste(eta[2])),
                             expression(paste(eta[3])),
                             expression(paste(eta[4])),
                             expression(paste(eta[5])),
                             expression(paste(eta[6])),
                             expression(paste(eta[7])),
                             expression(paste(eta[8])),
                             expression(paste(eta[9])),
                             expression(paste(eta[10])),
                             expression(paste(eta[11])),
                             expression(paste(eta[12])),
                             expression(paste(eta[13])),
                             expression(paste(eta[14])),
                             expression(paste(eta[15])),
                             expression(paste(eta[16])),
                             expression(paste(eta[17])),
                             expression(paste(eta[18])),
                             expression(paste(eta[19])),
                             expression(paste(eta[20]))))



#MCMC
  model = "model {
  for(i in 1:n){
    y[i] ~ dpois(lambda[i])
    log(lambda[i]) = beta1*x[i] + u[i]
  }

  beta1 ~ dt(0,1,5)
  for (j in 1:n){
  u[j] ~ dnorm(0,4)
  }
}"

library(runjags)
  tic("MCMC")
res.1 = run.jags(model = model,
                 monitor = c("beta1", "u"),
                 data = list(y = aa, n = N, x = cov),
                 n.chains = 6,
                 inits = list(beta1 = -0.5, u = rep(0.5,N)),
                 burnin = 10^2,
                 sample = 10^5,
                 method = "parallel")

trace <- combine.mcmc(res.1)
toc()
mean.mcmc <- c(mean(trace[, "beta1"]), mean(trace[,"u[1]"]), mean(trace[,"u[2]"]),
               mean(trace[,"u[3]"]), mean(trace[,"u[4]"]), mean(trace[,"u[5]"]), mean(trace[,"u[6]"]),
               mean(trace[,"u[7]"]), mean(trace[,"u[8]"]), mean(trace[,"u[9]"]), mean(trace[,"u[10]"]),
               mean(trace[,"u[11]"]), mean(trace[,"u[12]"]), mean(trace[,"u[13]"]), mean(trace[,"u[14]"]),
               mean(trace[,"u[15]"]), mean(trace[,"u[16]"]), mean(trace[,"u[17]"]), mean(trace[,"u[18]"]),
               mean(trace[,"u[19]"]), mean(trace[,"u[20]"]))
print(mean.mcmc)

dev.off()
par(mar = c(4,2,1,1))
hist(trace[, "beta1"], col = "azure1", prob = T, breaks = 50, xlab = expression(paste(beta[0])), main = "", ylab = "")
lines(x0, dnorm(x0, mean = mode1[1], sd = sqrt(cov1[1,1]), log = F), lty = 2, col = "black", lwd = 2)
lines(x0,dnorm(x0, mean = improved_mean[1], sd = sqrt(cov1[1,1]), log = F), lty = 1, col = "blue", lwd = 2)
lines(x0, dt(x0, 0, df = 5), lty = 4, lwd = 2, col = "red")

hist(trace[, "beta1"], col = "azure1", prob = T, breaks = 50, xlab = expression(paste(beta[1])), main = "", ylab = "")
lines(x0, dnorm(x0, mean = mode1[2], sd = sqrt(cov1[2,2]), log = F), lty = 2, col = "black", lwd = 2)
lines(x0,dnorm(x0, mean = improved_mean[2], sd = sqrt(cov1[2,2]), log = F), lty = 1, col = "blue", lwd = 2)
lines(x0, dunif(x0, min = -3, max = 3), lty = 4, lwd = 2, col = "red")

for (i in 1:N) {
hist(trace[, paste0("u[",i,"]")], col = "azure1", prob = T, breaks = 50, xlab = expression(paste("u[",i,"]")), main = "", ylab = "")
lines(x0, dnorm(x0, mean = mode1[i+1], sd = sqrt(cov1[i+1,i+1]), log = F), lty = 2, col = "black", lwd = 2)
lines(x0,dnorm(x0, mean = improved_mean[i+1], sd = sqrt(cov1[i+1,i+1]), log = F), lty = 1, col = "blue", lwd = 2)
lines(x0, dnorm(x0, mean = 0, sd = 1), lty = 4, lwd = 2, col = "red")
}

res= rbind(t(m), t(improved_mean), t(mean.mcmc))
colnames(res) = c("beta1")
rownames(res) = c("gaussian", "vb.corrected", "mcmc")
print(res)



#############With overall intercept

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

res_GG <- inla(formula = aa ~ 1 + cov + f(group, model = "iid", constr = T),
               family = "poisson",
               data = data.frame(aa = aa, cov = cov, group = group),
               control.inla = list(strategy = "gaussian"),
               control.compute = list(config = T))


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
#mode1 = c(res_GG$summary.fixed$mode, res_GG$summary.random$group$mode)
#hess2 <- numDeriv::hessian(func=f1, x = mode1)
m = mode1

cov1 = solve(-hess2)

# m = c(res_GG$summary.fixed$mode, res_GG$summary.random$group$mode)
# cov1a = matrix(res_GG$misc$configs$config[[1]]$Qinv[,N+1:2], byrow = F, nrow = N+2)
# cov1a = rbind(cov1a[N+1:2,], cov1a[1:N,])
# cov1=cov1a
# cov1[2,1] <- cov1[1,2]
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

par(mar = c(2,2,1,1))
plot(c(cov1[,1:2]%*%d_est$par, delta_eta), pch = c(rep(2,N+2),rep(16, N)), cex = c(rep(1,N+2), rep(0.5, N)), xaxt = "n", xlab = "", ylab = expression(paste(delta)),
     ylim = c(-0.2, 0.2))
abline(h=0, lty = 2)
# axis(1, at=1:(N+N+2), labels=c(expression(paste(beta[0])),
#                                expression(paste(beta[1])),
#                                expression(paste(u[1])),
#                                expression(paste(u[2])),
#                                expression(paste(u[3])),
#                                expression(paste(u[4])),
#                                expression(paste(u[5])),
#                                expression(paste(u[6])),
#                                expression(paste(u[7])),
#                                expression(paste(u[8])),
#                                expression(paste(u[9])),
#                                expression(paste(u[10])),
#                                expression(paste(u[11])),
#                                expression(paste(u[12])),
#                                expression(paste(u[13])),
#                                expression(paste(u[14])),
#                                expression(paste(u[15])),
#                                expression(paste(u[16])),
#                                expression(paste(u[17])),
#                                expression(paste(u[18])),
#                                expression(paste(u[19])),
#                                expression(paste(u[20])),
#                                expression(paste(eta[1])),
#                                expression(paste(eta[2])),
#                                expression(paste(eta[3])),
#                                expression(paste(eta[4])),
#                                expression(paste(eta[5])),
#                                expression(paste(eta[6])),
#                                expression(paste(eta[7])),
#                                expression(paste(eta[8])),
#                                expression(paste(eta[9])),
#                                expression(paste(eta[10])),
#                                expression(paste(eta[11])),
#                                expression(paste(eta[12])),
#                                expression(paste(eta[13])),
#                                expression(paste(eta[14])),
#                                expression(paste(eta[15])),
#                                expression(paste(eta[16])),
#                                expression(paste(eta[17])),
#                                expression(paste(eta[18])),
#                                expression(paste(eta[19])),
#                                expression(paste(eta[20]))))



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

# dev.off()
# par(mar = c(4,2,1,1))
# hist(trace[, "beta0"], col = "azure1", prob = T, breaks = 50, xlab = expression(paste(beta[0])), main = "", ylab = "")
# lines(x0, dnorm(x0, mean = mode1[1], sd = sqrt(cov1[1,1]), log = F), lty = 2, col = "black", lwd = 2)
# lines(x0,dnorm(x0, mean = improved_mean[1]+0.05, sd = sqrt(cov1[1,1]), log = F), lty = 1, col = "blue", lwd = 2)
# lines(x0, dt(x0, 0, df = 5), lty = 4, lwd = 2, col = "red")
# 
# hist(trace[, "beta1"], col = "azure1", prob = T, breaks = 50, xlab = expression(paste(beta[1])), main = "", ylab = "")
# lines(x0, dnorm(x0, mean = mode1[2], sd = sqrt(cov1[2,2]), log = F), lty = 2, col = "black", lwd = 2)
# lines(x0,dnorm(x0, mean = improved_mean[2]-0.05, sd = sqrt(cov1[2,2]), log = F), lty = 1, col = "blue", lwd = 2)
# lines(x0, dunif(x0, min = -3, max = 3), lty = 4, lwd = 2, col = "red")
# 
# for (i in 1:N) {
#   hist(trace[, paste0("u[",i,"]")], col = "azure1", prob = T, breaks = 50, xlab = expression(paste("u[",i,"]")), main = "", ylab = "")
#   lines(x0, dnorm(x0, mean = mode1[i+2], sd = sqrt(cov1[i+2,i+2]), log = F), lty = 2, col = "black", lwd = 2)
#   lines(x0,dnorm(x0, mean = improved_mean[i+2], sd = sqrt(cov1[i+2,i+2]), log = F), lty = 1, col = "blue", lwd = 2)
#   lines(x0, dnorm(x0, mean = 0, sd = 1), lty = 4, lwd = 2, col = "red")
# }

res= rbind(t(m[1:2]), t(improved_mean[1:2]), t(mean.mcmc[1:2]))
colnames(res) = c("beta0","beta1")
rownames(res) = c("gaussian", "vb.corrected", "mcmc")
print(res)


####Stan using BRMS
library(rstan)
library(brms)
library(tidybayes)
library(dplyr)


tic("STAN")
# make_stancode(y ~ 1 + xx + (1|group),
#               family = poisson,
#               data = data.frame(y = aa, xx = cov, group = group),
#               refresh = 0, silent = 2, cores = 4,
#               prior = c(
#                 prior(normal(0,1), class = Intercept),
#                 prior(normal(0.5,0.001), class = sd)
#               ))

# m_code <- 'functions {
# }
# data {
#   int<lower=1> N;  // total number of observations
#   array[N] int Y;  // response variable
#   int<lower=1> K;  // number of population-level effects
#   matrix[N, K] X;  // population-level design matrix
#   int<lower=1> Kc;  // number of population-level effects after centering
#   // data for group-level effects of ID 1
#   int<lower=1> N_1;  // number of grouping levels
#   int<lower=1> M_1;  // number of coefficients per level
#   array[N] int<lower=1> J_1;  // grouping indicator per observation
#   // group-level predictor values
#   vector[N] Z_1_1;
#   int prior_only;  // should the likelihood be ignored?
# }
# transformed data {
#   matrix[N, Kc] Xc;  // centered version of X without an intercept
#   vector[Kc] means_X;  // column means of X before centering
#   for (i in 2:K) {
#     means_X[i - 1] = mean(X[, i]);
#     Xc[, i - 1] = X[, i] - means_X[i - 1];
#   }
# }
# parameters {
#   vector[Kc] b;  // regression coefficients
#   real Intercept;  // temporary intercept for centered predictors
#   vector<lower=0>[M_1] sd_1;  // group-level standard deviations
#   array[M_1] vector[N_1] z_1;  // standardized group-level effects
# }
# transformed parameters {
#   vector[N_1] r_1_1;  // actual group-level effects
#   real lprior = 0;  // prior contributions to the log posterior
#   r_1_1 = (sd_1[1] * (z_1[1]));
#   lprior += normal_lpdf(Intercept | 0, 1);
#   lprior += normal_lpdf(sd_1 | 0.5, 0.001)
#   - 1 * normal_lccdf(0 | 0.5, 0.001);
# }
# model {
#   // likelihood including constants
#   if (!prior_only) {
#     // initialize linear predictor term
#     vector[N] mu = rep_vector(0.0, N);
#     mu += Intercept;
#     for (n in 1:N) {
#       // add more terms to the linear predictor
#       mu[n] += r_1_1[J_1[n]] * Z_1_1[n];
#     }
#     target += poisson_log_glm_lpmf(Y | Xc, mu, b);
#   }
#   // priors including constants
#   target += lprior;
#   target += std_normal_lpdf(z_1[1]);
# }
# generated quantities {
#   // actual population-level intercept
#   real b_Intercept = Intercept - dot_product(means_X, b);
# }'

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
               warmup = 1000, iter = 10000)
toc()

samplesss <- data.frame(mod_sim %>%
  spread_draws(r_group[group,]))
fixed_e <- data.frame(mod_sim %>%
                        spread_draws(c(b_Intercept,b_xx)))

stan_sam <- matrix(data = rep(NA, 36000*N/num_pg), ncol = N/num_pg, nrow = 36000)
for (i in 1:(N/num_pg)){
  stan_sam[,i] <- samplesss[samplesss$group == i, 2]
}

dev.off()
par(mar = c(4,2,1,1))
hist(trace[, "beta0"], col = "azure1", prob = T, breaks = 50, xlab = expression(paste(beta[0])), main = "", ylab = "", ylim = c(0,10), xlim = c(-1.3,-0.8))
hist(fixed_e$b_Intercept, col = rgb(255,100,100, max = 255,alpha = 50), prob = T, breaks = 50, add = T)
lines(x0, dnorm(x0, mean = mode1[1], sd = sqrt(cov1[1,1]), log = F), lty = 2, col = "black", lwd = 3)
lines(x0,dnorm(x0, mean = improved_mean[1], sd = sqrt(cov1[1,1]), log = F), lty = 1, col = "blue", lwd = 4)
lines(x0, dunif(x0, min = -3, max = 3), lty = 4, lwd = 4, col = "red")

hist(trace[, "beta1"], col = "azure1", prob = T, breaks = 50, xlab = expression(paste(beta[1])), main = "", ylab = "", ylim = c(0,10), xlim = c(-0.8,-0.3))
hist(fixed_e$b_xx, col = rgb(255,100,100, max = 255,alpha = 50), prob = T, breaks = 50, add = T)
lines(x0, dnorm(x0, mean = mode1[2], sd = sqrt(cov1[2,2]), log = F), lty = 2, col = "black", lwd = 4)
lines(x0,dnorm(x0, mean = improved_mean[2], sd = sqrt(cov1[2,2]), log = F), lty = 1, col = "blue", lwd = 4)
lines(x0, dt(x0, 0, df = 5), lty = 4, lwd = 4, col = "red")

 for (i in c(1,8,15)) {
   hist(trace[, paste0("u[",i,"]")], col = "azure1", prob = T, breaks = 50, xlab = expression(paste("u[",i,"]")), main = "", ylab = "")
   hist(stan_sam[,i], col = rgb(255,100,100, max = 255,alpha = 50), prob = T, breaks = 50, add = TRUE)
   lines(x0, dnorm(x0, mean = mode1[i+2], sd = sqrt(cov1[i+2,i+2]), log = F), lty = 2, col = "black", lwd = 2)
   lines(x0,dnorm(x0, mean = improved_mean[i+2], sd = sqrt(cov1[i+2,i+2]), log = F), lty = 1, col = "blue", lwd = 2)
   lines(x0, dnorm(x0, mean = 0, sd = 1), lty = 4, lwd = 2, col = "red")
 }

# plot(mod_sim)
# summary(mod_sim)
# mod_sim$fit
# plot(mod_sim, variable = "^r", regex = T)

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


