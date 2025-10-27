library(rjags)
library(runjags)
library(INLA)
library(sn)
library(mvtnorm)
library(e1071)
library(pracma)

# Simulate data-------------------------------------------------------------------------------------------------
set.seed(10)
n <- 100
beta1 <- -2 #-2
beta2 <- -3 #-3
beta3 <- -3
beta4 <- 1
No_of_fixed_effects <- n + 4
rho <- sqrt(3)/2

covariate1 <- rep(1, n) 
covariate2 <- rnorm(n, mean = 0, sd = 0.5)
covariate3 <- rnorm(n, mean = 0.2, sd = 0.25)
covariate4 <- rnorm(n, mean = 0.1, sd = 0.5)

s <- 1
prec <- 1/s^2
S0 <- s^2 * toeplitz(rho^seq(0, n-1))
Q0 <- solve(S0)

u <-  rmvnorm(1, sigma = S0)[1,] 
eta <- beta1 * covariate1 + beta2 * covariate2 + beta3 * covariate3 + beta4 * covariate4 + u
mu <- exp(eta) / (1 + exp(eta))

Ntrials <- 2
y <-  rbinom (n = n, prob = mu, size = Ntrials)
table(y)

idx = 1:n
prior.ui <- list(prec = list(prior = "loggamma", param = c(1, 1), initial = log(prec), fixed = T),
                 rho = list(prior = "loggamma", param = c(1, 1), initial = log((1+rho)/(1-rho)), fixed = T))


# INLA VB-M+V----------------------------------------------------------------------------------------------------
formula = y ~ -1 + x1 + x2 + x3 + x4 + f(idx,model="ar1",hyper = prior.ui)
timeI1 <- Sys.time()

res = inla(formula = formula,
       family = "binomial",
       data = data.frame(y=y, x1 = covariate1, x2 = covariate2, x3 = covariate3, x4 = covariate4 ,idx = idx, rep(Ntrials, n)),
       Ntrials = rep(Ntrials, n),
       control.compute = list(config=T, return.marginals.predictor=T),
       control.predictor = list(compute = T),
       control.inla = list(control.vb = list(strategy = "variance")), 
       control.fixed=list(prec.intercept = 0.01, prec = 0.01))
timeI2 <- Sys.time()
time_INLA <- timeI2 - timeI1

covariate <- res$misc$configs$A
covariate <- covariate[,c(n + (1:4), 1:n)]

# MCMC -----------------------------------------------------------------------------------------------------------
model = "model {
	for(i in 1:n){
	y[i] ~ dbinom(p[i], Ntrials)
	logit(p[i]) =  beta1 * u1[i] + beta2 * u2[i] + beta3 * u3[i] + beta4 * u4[i] + u[i]
	}

	beta1 ~ dnorm(0, 0.01)
	beta2 ~ dnorm(0, 0.01)
	beta3 ~ dnorm(0, 0.01)
	beta4 ~ dnorm(0, 0.01)
	u[1] ~ dnorm(0, PREC)
	for(i in 2:n){
		u[i] ~ dnorm(RHO * u[i-1], PREC / (1 - RHO * RHO))
	}
}" 

time_MCMC1 <- Sys.time()
trace = combine.mcmc(run.jags(model = model,
			 monitor = c( "beta1", "beta2", "beta3",  "beta4"),
			 data = list(y = y, u1 = covariate1, u2 = covariate2, u3 = covariate3, u4 = covariate4, n = n, tau = 0.25, Ntrials = Ntrials, RHO = rho, PREC = prec),
			 n.chains = 2,
			 thin = 20,
			 inits = list(beta0 = 0, beta1 = 0, ".RNG.name"="base::Super-Duper", ".RNG.seed"=7),
			 burnin = 10^3,
			 sample = 10^4,
			 method = "parallel"))
time_MCMC2 <- Sys.time()
time_MCMC <- time_MCMC2 - time_MCMC1

beta_skew <- apply(trace, 2, skewness)

# Helper functions for SGC-VB-------------------------------------------------------------------------------------------

# Function for summing up skew normal
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
    if(m %% 2 == 0) sum.f <- ifftshift(sum.f)
    return (list(x = x + sum.mu, y = sum.f, mean = sum.mu, var = sum.var))
}

# Function for square root of eigen value decomposition
mat.power <- function(A, power = 1) {
    e <- eigen(A)
    return (
        matrix(e$vectors %*% diag(e$values^power) %*% t(e$vectors),
               nrow(A), ncol(A)))
}

# Function for re-ordering the variance co-variance matrix
re_ordered_mat <- function(S, to_correct){
	dim_of_mat <- dim(S)[1]
	new_order <- c(to_correct, setdiff(1:dim_of_mat, to_correct))
	new_mat <- S[new_order, new_order]
	return(new_mat)
}

# Main function in SGC-VB
# Function for optimizing the skewness
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

# Extract first two posterior moments-----------------------------------------------------------------------
### Covariance matrix from INLA
S0 <- solve(forceSymmetric(res$misc$configs$config[[1]]$Q))
S0 <- S0[c(n + 1:4, 1:n), c(n + 1:4, 1:n)]
S0 <- as.matrix(S0)
### Correlation matrix from INLA
CORR_MAT <- cov2cor(S0)

# Mean 
mu <- res$misc$configs$config[[1]]$improved.mean
mu <- mu[c(n + 1:4, 1:n)]

# To declare which components to correct--------------------------------------------------------------------
to_be_corrected <- c(1:4)

# Declaring an empty vector to store the results------------------------------------------------------------
beta_skew_est <- numeric(length(to_be_corrected))

# Main optimization function--------------------------------------------------------------------------------
timeS1 <- Sys.time()
for(i in seq_along(to_be_corrected)) {

	print(i)
	
	S <- S0
	
	to_correct <- order(abs(CORR_MAT[,to_be_corrected[i]]), decreasing = TRUE)
	to_correct <- to_correct[1:ifelse(length(to_be_corrected) > 4, 4, length(to_be_corrected))]
	len_to_correct <- 1 #correct one at a time
	len_not_to_correct <- No_of_fixed_effects - len_to_correct
	
	S <- re_ordered_mat(S, to_correct)
	
	S11 <- S[1:len_to_correct, 1:len_to_correct] # components to correct
	S22 <- S[len_to_correct + 1:len_not_to_correct, len_to_correct + 1:len_not_to_correct] 
		# components not to correct
	S21 <- S[len_to_correct + 1:len_not_to_correct, 1:len_to_correct] 
		# Covariance part of to_correct and not_to_correct

	L11 <- t(chol(S11))
	L11_inv <- solve(L11)
	
	L22 <- S22 - S21 %*% solve(S11) %*% t(S21)
	L22_eigen <- mat.power(S22 - S21 %*% solve(S11) %*% t(S21), -1/2)
	L21 <- -L22_eigen %*% S21 %*% solve(S11)
	L12 <- matrix(0, nrow = len_to_correct, ncol = len_not_to_correct)
	
	LInv <- rbind(cbind(L11_inv, L12), cbind(L21, L22_eigen))
	L <- solve(LInv)
	
	mu_new <- mu[c(to_correct, setdiff(1:No_of_fixed_effects, to_correct))]
	mu_ind <- LInv %*% mu_new
	
	len_to_correct <- length(to_correct)
	
	r <- optim(rep(1e-2, length(to_correct)), 
		VB_skew, 
		NULL,
		method = "BFGS",
		control = list(ndeps = rep(1/1000, length(to_correct))),
		to_correct = to_correct)
	
	beta_skew_est[i] <- r$par[1]
}
timeS2 <- Sys.time()
time_SC <- timeS2- timeS1

# Comparing the results between MCMC, INLA VB-M+V and SGC-VB-----------------------------------------------
results <- data.frame(beta_skew_est, beta_skew[to_be_corrected])
colnames(results) <- c("SGC", "MCMC")

mean_of_corrected <- mu[to_be_corrected]
variance_of_corrected <- diag(S0[to_be_corrected,to_be_corrected])

beta_mean <- apply(trace, 2, mean)
beta_var <- (apply(trace, 2, sd))^2

results$MCMC_mean <- beta_mean[to_be_corrected]
results$INLA_mean <- mean_of_corrected
results$MCMC_var <- beta_var[to_be_corrected]
results$INLA_var <- variance_of_corrected
print("n=100, iter = 10e3")
print(round(results, 2))
time_INLA
time_INLA+time_SC
time_MCMC
