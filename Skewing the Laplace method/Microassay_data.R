library(INLA)
library(rstan)
library(tictoc)

#https://onlinelibrary.wiley.com/doi/full/10.1002/bimj.201600207

load("/Users/vanniej/Downloads/bimj1793-sup-0002-codeanddata/data/base/Datengrundlage_unbalanciert_10vs80.RData")

data10$class <- as.numeric(class10)-1
data10$ID <- extID10

summary(data10)
summary(data10$class)

n = nrow(data10)

#INLA
res1_bimj <- inla(class ~ 1 + Biomarker1 + Biomarker2 + Biomarker3,
                  data = data10,
                  family = "binomial",
                  control.compute = list(config=T, return.marginals.predictor=T),
                  control.predictor = list(compute = T),
                  control.inla = list(control.vb = list(strategy = "variance")),
                  control.fixed=list(prec.intercept = 0.01, prec = 0.01))
summary(res1_bimj)

##MCMC

stanmodelcode <- "
data {
  int<lower=0> N;
  vector[N] x1;
vector[N] x2;
vector[N] x3;
  array[N] int<lower=0, upper=1> y;
}
parameters {
  real alpha;
  real beta1;
real beta2;
real beta3;
}
model {
  y ~ bernoulli_logit(alpha + beta1 * x1 + beta2 * x2 + beta3 * x3);
}"

data_stan <- data10[c(-12)]

N = nrow(data_stan)
colnames(data_stan) <- c(paste0("x",1:10),"y")
data_stan <- as.list(data_stan)
data_stan$N <- N
tic()
res_stan <- stan(model_code = stanmodelcode,
                 model_name = "Imbalanced data_logit",
                 data = data_stan,
                 iter = 10000,
                 chains = 3)
toc()
samps <- rstan:::extract(res_stan)
mcmc.skew <- rep(NA,11)

skew <- function(x, na.rm = FALSE){
  if(na.rm) x <- x[!is.na(x)]
  n <- length(x)
  sum((x - mean(x))^3)/(n - 2)/var(x)^(3/2)
}

mcmc.skew[1] <- skew(samps$alpha)
mcmc.skew[2] <- skew(samps$beta1)
mcmc.skew[3] <- skew(samps$beta2)
mcmc.skew[4] <- skew(samps$beta3)


##SKEW corr
library(parallel)
library(Matrix)
library(statmod)
library(runjags)
library(rjags)
library(INLA)
library(sn)
library(numDeriv)
library(mvtnorm)
library(gtools)
library(e1071)
library(parallel)
library(LaplacesDemon)
library(SynchWave)

No_of_fixed_effects = 4
covariate <- res1_bimj$misc$configs$A

sum.of.sn <- function(mu, var, skew) {
  m <- length(mu)

  stopifnot(length(var) == m)
  if (missing(skew)) skew <- rep(0, m)
  stopifnot(length(skew) == m)

  sum.mu <- sum(mu)
  sum.var <- sum(var)
  ## don't know how many points are required here...
  n <- as.integer(2L^10L)

  ## need to get one element that is exactly zero
  x <- seq.int(-n, n-1L, by = 2L) / n
  x <- x * 6.0 * sqrt(sum.var)  ##  need to choose the range...
  idx.zero <- n %/% 2L+1L
  stopifnot(x[idx.zero] == 0.0)

  sum.fft <- rep(1, n)
  for(k in 1:m) {
    if (var[k] == 0) {
      p <- rep(0, n)
      p[idx.zero] <- 1
    } else {
      if (skew[k] == 0) {
        p <- dnorm(x, sd = sqrt(var[k]))
      } else {
        par <- INLA:::inla.sn.reparam(moments = c(0, var[k], skew[k]))
        p <- dsn(x, xi = par$xi, omega = par$omega, alpha = par$alpha)
      }
      p <- p / sum(p)
    }

    sum.fft <- sum.fft * (fft(p))
  }

  sum.f <- Re(fft(sum.fft, inverse = TRUE))
  sum.f <- sum.f / (sum(sum.f) * (x[2] - x[1]))
  if(even(m)) sum.f <- ifftshift(sum.f)
  return (list(x = x + sum.mu, y = sum.f, mean = sum.mu, var = sum.var))
}


# Function for square root of eigen value decomposition
mat.power <- function(A, power = 1) {
  e <- eigen(A)
  return (
    matrix(e$vectors %*% diag(e$values^power) %*% t(e$vectors),
           nrow(A), ncol(A)))
}


# Function for re-ordering the varaince co-variance matrix
re_ordered_mat <- function(S, to_correct){
  dim_of_mat <- dim(S)[1]
  new_order <- c(to_correct, setdiff(1:dim_of_mat, to_correct))
  new_mat <- S[new_order, new_order]
  return(new_mat)
}

#I2 numerical integration
# fun <- function(x1, x2, x3, x4, Sigma, skew) {
#   ## vectorize in the first argument. x2 can be a vector that is kept fixed
#   val <- numeric(length(x1))
#   Q <- solve(Sigma)
#   par1 <- INLA:::inla.sn.reparam(moments = c(0, Sigma[1, 1], skew[1]))
#   par2 <- INLA:::inla.sn.reparam(moments = c(0, Sigma[2, 2], skew[2]))
#   par3 <- INLA:::inla.sn.reparam(moments = c(0, Sigma[3, 3], skew[3]))
#   par4 <- INLA:::inla.sn.reparam(moments = c(0, Sigma[4, 4], skew[4]))
#   for(idx in 1:length(x1)) {
#     xx1 <- x1[idx]
#     for (jdx in 1:length(x2)) {
#       xx2 <- x2[jdx]
#       for (kdx in 1:length(x3)) {
#         xx3 <- x3[kdx]
#         x12 <- matrix(cbind(xx1, xx2, xx3, x4), ncol = 4)
#     d <- dmvnorm(x12, sigma = Sigma, log = FALSE)
#     z12 <- cbind(qsn(pnorm(x12[, 1] / sqrt(Sigma[1, 1])), dp = unlist(par1), solver = "RFB"),
#                  qsn(pnorm(x12[, 2] / sqrt(Sigma[2, 2])), dp = unlist(par2), solver = "RFB"),
#                  qsn(pnorm(x12[, 3] / sqrt(Sigma[3, 3])), dp = unlist(par3), solver = "RFB"),
#                  qsn(pnorm(x12[, 4] / sqrt(Sigma[4, 4])), dp = unlist(par4), solver = "RFB"))
#     vals <- (unlist(lapply(1:nrow(x12),
#                            function(idx) {
#                              xx12 <- x12[idx, ]
#                              zz12 <- z12[idx, ]
#                              return (-0.5 * t(xx12) %*% Q %*% xx12 + 0.5 * t(zz12) %*% Q %*% zz12)
#                            })))
#     val[idx] <- sum(d * vals)
#   }
#   return (val)
#     }
#   }
# }
#
# x1 <- seq(-5, 5, by = 0.025)
# x2 <- seq(-5, 5, by = 0.025)
# x3 <- seq(-5, 5, by = 0.025)
# x4 <- seq(-5, 5, by = 0.025)
#
# rho <- 0.75
# Sigma = matrix(rho, ncol = 4, nrow = 4)
# diag(Sigma) <- rep(1, 4)
# Q <- solve(Sigma)
# skew <- c(0.3, -0.4, 0.6, -0.2)
#
# I2 <- sum(fun(x1, x2, x3, x4, Sigma, skew)) * diff(x1)[1] * diff(x2)[1]* diff(x3)[1] * diff(x4)[1]
# I2
d_SGC <- function(x, Sigma1, par1){

  p0 = length(x)
  d0 = 1
  d1 = 1
  vec0 = c(0)

  for (i in 1:p0){
    d0 = d0 * dsn(x[i], dp = unlist(par1[[i]]))
    d1 = d1 * dnorm(qnorm(psn(x[i], dp = unlist(par1[[i]]))))
    vec0 = c(vec0, qnorm(psn(x[i], dp = unlist(par1[[i]]))))
  }

  vec0 = vec0[-1]

  return(d0*dmvnorm(x = as.vector(vec0), sigma = cov2cor(matrix(Sigma1, ncol = 4)))/d1)


}

kld_f = function(mu, Sigma, s){
  # mu = mu0
  # Sigma = Sigma0
  # s = s0

  #Make grid of values and par sets
  p = length(mu)
  NN = 5
  x0 = list()
  step = 1
  par = list()
  for (ii in 1:p){
    x1 = seq(mu[ii]-1*sqrt(Sigma[ii,ii]), mu[ii]+1*sqrt(Sigma[ii,ii]), length = NN)
    x0 = c(x0, list(x1))
    par = c(par, list(INLA:::inla.sn.reparam(moments = c(mu[ii], Sigma[ii,ii], s[ii]))))
    step = step*diff(x1)[1]
  }

  grid0 = as.matrix(drop(expand.grid(x0)), ncol = p, nrow = NN^p)

  sum = 0

  for (jj in 1:(NN^p)){
    snjj = as.numeric(d_SGC(x = c(grid0[jj,]), Sigma1 = Sigma, par1 = par))
    sum = sum + (snjj*log(snjj/as.numeric(dmvnorm(grid0[jj,], sigma = Sigma, mean = mu))))
  }

  return(sum*step)

}

# Function for optimising the skewness
VB_skew <- function(SKEW, to_correct) {

  #  print(SKEW)
  s <- rep(0, No_of_fixed_effects)
  s[1:len_to_correct] <- SKEW


  kld = kld_f(mu = mu_ind, Sigma = S, s = s)

  negative_log_likelihood_sum <- 0

  for(k in 1:n){

    mu1 <- mu_ind

    coefficient <- c(covariate[k,])
    coefficient <- coefficient[c(to_correct, setdiff(1:No_of_fixed_effects, to_correct))]
    coefficient <- as.vector(coefficient %*% L)

    mu1 <- coefficient * mu1
    mu2 <- sum(mu1[-c(1:len_to_correct)])  # clubbing all the normal distributions
    mu1 <- c(mu1[c(1:len_to_correct)], mu2)

    S1 <- coefficient^2
    S2 <- sum(S1[-c(1:len_to_correct)]) # clubbing all the normal distributions
    S1 <- c(S1[c(1:len_to_correct)], S2)


    #s <- c(SKEW[1], rep(0,9), SKEW[2:3])
    s <- rep(0, No_of_fixed_effects)
    s[1: len_to_correct] <- SKEW
    s <- sign(coefficient) * s
    s <- c(s[1: len_to_correct], 0)

    dist <- sum.of.sn(mu = mu1, var = S1, skew = s)

    probability <- exp(dist$x)/(1 + exp(dist$x))

    if(any(probability > 1 - 1e-7)){
      to.replace <-  which(probability > 1-1e-7)
      probability[to.replace] = 1 - 1e-7
    }


    expected_log_like <- sum(-dbinom(y[k], size = Ntrials, prob = probability,
                                     log = TRUE) * dist$y * (dist$x[2] - dist$x[1]))

    negative_log_likelihood_sum <- negative_log_likelihood_sum + expected_log_like

  }

  return(negative_log_likelihood_sum + kld)
}


### Covariance matrix from INLA
S0 <- matrix(solve(forceSymmetric(res1_bimj$misc$configs$config[[1]]$Q)), ncol = No_of_fixed_effects)
S0

### Covariance matrix from INLA
CORR_MAT <- cov2cor(S0)
CORR_MAT

# To declare which components to correct
to_be_corrected <- 1


# Declaring an empty vector to store the results
beta_skew_est <- numeric(length(to_be_corrected))

y = data10$class
Ntrials = 1

# For loop for the main optimisation function
tic()
for(i in seq_along(to_be_corrected)){

  print(i)

  S <- S0

  to_correct <- order(abs(CORR_MAT[,to_be_corrected[i]]), decreasing = TRUE)
  to_correct <- to_correct[1:ifelse(length(to_be_corrected) > 2, 2, length(to_be_corrected))]
  len_to_correct <- 1 #length(to_correct)
  len_not_to_correct <- No_of_fixed_effects - len_to_correct

  S <- re_ordered_mat(S, to_correct)

  S11 <- S[1:len_to_correct, 1:len_to_correct] # components to correct
  S22 <- S[len_to_correct + 1:len_not_to_correct, len_to_correct + 1:len_not_to_correct]
  # components not to correct
  S21 <- S[len_to_correct + 1:len_not_to_correct, 1:len_to_correct]
  # Covariace part of to_correct and not_to_correct

  L11 <- t(chol(S11))
  L11_inv <- solve(L11)

  L22 <- S22 - S21 %*% solve(S11) %*% t(S21)
  L22_eigen <- mat.power(S22 - S21 %*% solve(S11) %*% t(S21), -1/2)
  L21 <- -L22_eigen %*% S21 %*% solve(S11)
  L12 <- matrix(0, nrow = len_to_correct, ncol = len_not_to_correct)

  LInv <- rbind(cbind(L11_inv, L12), cbind(L21, L22_eigen))
  L <- solve(LInv)

  #	mu <- apply(trace, 2, mean) #res$misc$configs$config[[1]]$improved.mean
  mu <- res1_bimj$misc$configs$config[[1]]$improved.mean
  mu <- mu[c(1:4, 1:n)]
  mu_new <- mu[c(to_correct, setdiff(1:No_of_fixed_effects, to_correct))]
  mu_ind <- LInv %*% mu_new

  len_to_correct <- length(to_correct)
  r <- optim(rep(0.5, length(to_correct)),
             VB_skew,
             NULL,
             method = "BFGS",
             control = list(ndeps = rep(1/1000, length(to_correct))),
             to_correct = to_correct)

  beta_skew_est[i] <- r$par[1]
}
toc()
#Plots

mean_of_corrected <- res1_bimj$misc$configs$config[[1]]$improved.mean
varaince_of_corrected <- diag(S0[(1:4),(1:4)])

par(mar = c(4,4,0.1,0.1))

hist(samps$alpha, breaks = 50,main = "",  prob = T, ylim = c(0,0.5),
     xlab = expression(beta[0]), ylab = expression(paste("p(", beta[0],"|data)")))
lines(inla.smarginal(res1_bimj$marginals.fixed$`(Intercept)`), lwd = 2)
j = 1

if (j %in% to_be_corrected) {
skn_param <- INLA:::inla.sn.reparam(moments = c(mean_of_corrected[j], varaince_of_corrected[j], beta_skew_est[j]))
#

x_coordinate <- inla.smarginal(res1_bimj$marginals.fixed[[j]])$x
y_coordinate <- dsn(x_coordinate, xi = skn_param$xi, omega = skn_param$omega, alpha = skn_param$alpha)

lines(x_coordinate, y_coordinate, type = "l", lwd = 2, col = "red", lty = 2)
}

hist(samps$beta1, breaks = 50, main = "",  prob = T, ylim = c(0,0.6),
     xlab = expression(beta[1]), ylab = expression(paste("p(", beta[1],"|data)")))
lines(inla.smarginal(res1_bimj$marginals.fixed$Biomarker1), lwd = 2)

j = 2

if (j %in% to_be_corrected) {
skn_param <- INLA:::inla.sn.reparam(moments = c(mean_of_corrected[j], varaince_of_corrected[j], beta_skew_est[j]))
#
x_coordinate <- inla.smarginal(res1_bimj$marginals.fixed[[j]])$x
y_coordinate <- dsn(x_coordinate, xi = skn_param$xi, omega = skn_param$omega, alpha = skn_param$alpha)

lines(x_coordinate, y_coordinate, type = "l", lwd = 2, col = "red", lty = 2)
}

hist(samps$beta1, breaks = 50, main = "",  prob = T, ylim = c(0,0.3), xlim = c(0,2.5),
     xlab = expression(beta[1]), ylab = expression(paste("p(", beta[1],"|data)")))
lines(inla.smarginal(res1_bimj$marginals.fixed$Biomarker1), lwd = 2)

j = 2

if (j %in% to_be_corrected) {
  skn_param <- INLA:::inla.sn.reparam(moments = c(mean_of_corrected[j], varaince_of_corrected[j], beta_skew_est[j]))
  #
  x_coordinate <- inla.smarginal(res1_bimj$marginals.fixed[[j]])$x
  y_coordinate <- dsn(x_coordinate, xi = skn_param$xi, omega = skn_param$omega, alpha = skn_param$alpha)

  lines(x_coordinate, y_coordinate, type = "l", lwd = 2, col = "red", lty = 2)
}

hist(samps$beta2, breaks = 50, main = "",  prob = T, ylim = c(0,0.4),
     xlab = expression(beta[2]), ylab = expression(paste("p(", beta[2],"|data)")))
lines(inla.smarginal(res1_bimj$marginals.fixed$Biomarker2), lwd = 2)

j = 3
if (j %in% to_be_corrected) {
skn_param <- INLA:::inla.sn.reparam(moments = c(mean_of_corrected[j], varaince_of_corrected[j], beta_skew_est[j]))
#
x_coordinate <- inla.smarginal(res1_bimj$marginals.fixed[[j]])$x
y_coordinate <- dsn(x_coordinate, xi = skn_param$xi, omega = skn_param$omega, alpha = skn_param$alpha)

lines(x_coordinate, y_coordinate, type = "l", lwd = 2, col = "red", lty = 2)
}

hist(samps$beta3, breaks = 50, main = "",  prob = T, ylim = c(0,0.4),
     xlab = expression(beta[3]), ylab = expression(paste("p(", beta[3],"|data)")))
lines(inla.smarginal(res1_bimj$marginals.fixed$Biomarker3), lwd = 2)

j = 4
if (j %in% to_be_corrected) {
  skn_param <- INLA:::inla.sn.reparam(moments = c(mean_of_corrected[j], varaince_of_corrected[j], beta_skew_est[j]))
  #
  x_coordinate <- inla.smarginal(res1_bimj$marginals.fixed[[j]])$x
  y_coordinate <- dsn(x_coordinate, xi = skn_param$xi, omega = skn_param$omega, alpha = skn_param$alpha)

  lines(x_coordinate, y_coordinate, type = "l", lwd = 2, col = "red", lty = 2)
}

hist(samps$beta3, breaks = 50, main = "", prob = T, ylim = c(0,0.25), xlim = c(-0.9,1.5),
     xlab = expression(beta[3]), ylab = expression(paste("p(", beta[3],"|data)")))
lines(inla.smarginal(res1_bimj$marginals.fixed$Biomarker3), lwd = 2)

j = 4
if (j %in% to_be_corrected) {
skn_param <- INLA:::inla.sn.reparam(moments = c(mean_of_corrected[j], varaince_of_corrected[j], beta_skew_est[j]))
#
x_coordinate <- inla.smarginal(res1_bimj$marginals.fixed[[j]])$x
y_coordinate <- dsn(x_coordinate, xi = skn_param$xi, omega = skn_param$omega, alpha = skn_param$alpha)

lines(x_coordinate, y_coordinate, type = "l", lwd = 2, col = "red", lty = 2)
}





print("MCMC skewness")
mcmc.skew
print("INLA skewness")
beta_skew_est



