data1 <- read.csv("~/Google Drive/My Drive/KAUST/Research/VB/Skew and Var/diabetes_binary_health_indicators_BRFSS2015.csv")
data1 <- data1[sample(nrow(data1), size = 5000),]
library(rstan)
library(INLA)
time_I1 <- Sys.time()
res1 <- inla(Diabetes_binary ~ 1 + HighBP + HighChol +
                CholCheck + BMI + Smoker + Stroke + 
               HeartDiseaseorAttack + PhysActivity +
               Fruits + Veggies + HvyAlcoholConsump +
               AnyHealthcare + NoDocbcCost +
               GenHlth + MentHlth + PhysHlth + DiffWalk +
               Sex + Age + Education + Income,
             family = "binomial",
             data = data1,
             control.predictor = list(compute = T),
             control.compute = list(config=T, return.marginals.predictor=T),
             control.inla = list(control.vb = list(strategy = "variance")), 
             control.fixed=list(prec.intercept = 0.01, prec = 0.01)
             )
time_I2 <- Sys.time()
time_INLAG <- time_I2-time_I1

summary(res1)

library(rstan)
library(e1071)

stanmodelcode <- "
data {
  int<lower=0> N;
  vector[N] x1;
vector[N] x2;
vector[N] x3;
vector[N] x4;
vector[N] x5;
vector[N] x6;
vector[N] x7;
vector[N] x8;
vector[N] x9;
vector[N] x10;
  vector[N] x11;
vector[N] x12;
vector[N] x13;
vector[N] x14;
vector[N] x15;
vector[N] x16;
vector[N] x17;
vector[N] x18;
vector[N] x19;
vector[N] x20;
vector[N] x21;
  array[N] int<lower=0, upper=1> y;
}
parameters {
  real alpha;
real beta1;
real beta2;
real beta3;
real beta4;
real beta5;
real beta6;
real beta7;
real beta8;
real beta9;
real beta10;
real beta11;
real beta12;
real beta13;
real beta14;
real beta15;
real beta16;
real beta17;
real beta18;
real beta19;
real beta20;
real beta21;
}
model {
  y ~ bernoulli_logit(alpha + beta1 * x1 + beta2 * x2 + beta3 * x3 + beta4 * x4 + beta5 * x5 +
  beta6*x6 + beta7*x7 + beta8*x8 + beta9*x9 + beta10*x10 + beta11 * x11 + beta12 * x12 + beta13 * x13 + beta14 * x14 + beta15 * x15 +
  beta16*x16 + beta17*x17 + beta18*x18 + beta19*x19 + beta20*x20 +beta21*x21);
}"

data_stan <- data1[]

N = nrow(data_stan)

colnames(data_stan) <- c("y", paste0("x",1:21))
data_stan <- as.list(data_stan)
data_stan$N <- N

if (T){
time_M1 <- Sys.time()

  res_stan <- stan(model_code = stanmodelcode,
                   model_name = "Imbalanced data_logit",
                   data = data_stan,
                   iter = 10^5,
                   chains = 1)
time_M2 <- Sys.time()
time_M4 <- time_M2-time_M1

samps <- rstan:::extract(res_stan)

mcmc.skew <- apply(data.frame(samps[1:22]), 2, skewness)
}

#Skew corr
n = nrow(data1)
No_of_fixed_effects = length(res1$summary.fixed$mean)
covariate <- res1$misc$configs$A

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

# Function for optimising the skewness
VB_skew <- function(SKEW, to_correct) {
  
  #  print(SKEW)
  s <- rep(0, No_of_fixed_effects)
  s[1:len_to_correct] <- SKEW
  
  
  sk <- sum((0.1696405/factorial(2) * s^2) - (0.3984908/factorial(4) * s^4) 
            + (49.3041004/factorial(6) * s^6))
  
  negative_log_likelihood_sum <- 0
  
  #   browser()
  
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
  
  return(negative_log_likelihood_sum + sk)
}

## Automatic Version for optimisation

### Covariance matrix from INLA
S0 <- solve(forceSymmetric(res1$misc$configs$config[[1]]$Q))
S0 <- as.matrix(S0)

### Covariance matrix from INLA
CORR_MAT <- cov2cor(S0)

# To declare which components to correct
to_be_corrected <- c(1:4)
library(sn)
library(gtools)

# Declaring an empty vector to store the results
beta_skew_est <- numeric(length(to_be_corrected))

y = data1$Diabetes_binary
Ntrials = 1
library(SynchWave)
# For loop for the main optimisation function
timeS1 <- Sys.time()
for(i in seq_along(to_be_corrected)){
  
  print(i)
  
  S <- S0
  
  to_correct <- order(abs(CORR_MAT[,to_be_corrected[i]]), decreasing = TRUE)
  to_correct <- to_correct[1:ifelse(length(to_be_corrected) > 4, 4, length(to_be_corrected))]
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
  mu <- res1$misc$configs$config[[1]]$improved.mean
  mu <- mu[c(1:No_of_fixed_effects, 1:n)]
  mu_new <- mu[c(to_correct, setdiff(1:No_of_fixed_effects, to_correct))]
  mu_ind <- LInv %*% mu_new
  
  len_to_correct <- length(to_correct)
  r <- optim(rep(1e-2, length(to_correct)), 
             VB_skew, 
             NULL,
             method = "BFGS",
             control = list(ndeps = rep(1/1000, length(to_correct))),
             to_correct = to_correct
  
  beta_skew_est[i] <- r$par[1]
}
timeS2 <- Sys.time()
time_S <- timeS2-timeS1

