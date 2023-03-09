library(tictoc)

###Example in section 3.2

x0 <- seq(-10,5,by = 0.001)

f1 <- function(x){
  a = aa
  ll = 0
  for (i in 1:N){
    ll = ll + dpois(a[i], lambda = exp(x[1] + x[2]*cov[i]), log = T) 
  }
  ll = ll + dt(x[1], df = 5, log = TRUE) + dunif(x[2], min = -3, max = 3, log = TRUE) 
  return(ll)
}

mu = -1
beta1 = -0.5
N = 20
cov = rnorm(N)
aa = rpois(N, lambda  = exp(mu + beta1*cov))
par(mar = c(4,2,2,1))

barplot(table(as.factor(aa)), main = "", ylim = c(0, N))

tic("Gaussian")
if (T){
 tol = 2
 mode1 = c(4,1)
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
  A = matrix(data = cbind(rep(1, N), cov), ncol = 2)
  eta = matrix(data = cbind(rep(1, N), cov), ncol = 2)%*%(m + delta)
  cov_eta = A%*%cov1%*%t(A)
  
  ll = 0
  for (i in 1:N){
  loglik = exp(eta[i] + cov_eta[i,i]/2) - 
    aa[i]*eta[i]
  ll = ll + loglik}
  
  #KLD
  xx = rbind(rep(seq(-3,  1, length = 20), each = 20),rep(seq( - 3,  3, length = 20), n = 20))
  area = (max(xx[1,]) - min(xx[1,]))*(max(xx[2,]) - min(xx[2,]))
  prior = rep(NA, ncol(xx))
  post = prior
  f = prior
  f2 = prior
  for (i in 1:ncol(xx)){
  prior[i] = dt(xx[1,i], df = 5)*dunif(xx[2,i], min = -3, max = 3)#dnorm(xx[1,i], mean = 0, sd = 1, log = FALSE)#dnorm(xx[1,i], mean = 0, sd = 1, log = F)*dnorm(xx[2,i], mean = 0, sd = 1, log = F)
  post[i] = dmvnorm(xx[,i], mean = m + delta, sigma = cov1, log = F)
  f[i] = post[i]*(log(post[i]+0.001) - log(prior[i]+0.001))
  }
  kl = area*sum(f)/ncol(xx)
  return (ll + kl)
}


tic("VBC")
d_est = optim(c(0,0), minf)
toc()
d_est$par
improved_mean = m + d_est$par
improved_mean_eta = cbind(rep(1, N), cov)%*%m + cbind(rep(1, N), cov)%*%d_est$par
delta_eta = cbind(rep(1, N), cov)%*%d_est$par

par(mar = c(2,2,1,1))
plot(c(delta, delta_eta), pch = c(rep(2,2),rep(16, N)), cex = c(rep(1,2), rep(0.5, N)), xaxt = "n", xlab = "", ylab = expression(paste(delta)),
     ylim = c(-0.2, 0))
abline(h=0, lty = 2)
axis(1, at=1:(N+2), labels=c(expression(paste(beta[0])), 
                             expression(paste(beta[1])),
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
    log(lambda[i]) = beta0 + beta1*x[i]
  }

  beta0 ~ dt(0,1,5)
  beta1 ~ dunif(-3,3)
}"

library(runjags)
  tic("MCMC")
res.1 = run.jags(model = model,
                 monitor = c("beta0", "beta1"),
                 data = list(y = aa, n = N, x = cov),
                 n.chains = 6,
                 inits = list(beta0 = -1, beta1 = -0.5),
                 burnin = 10^2,
                 sample = 10^5,
                 method = "parallel")

trace <- combine.mcmc(res.1)
toc()
mean.mcmc <- c(mean(trace[, "beta0"]), mean(trace[, "beta1"]))
print(mean.mcmc)

dev.off()
par(mar = c(4,2,1,1))
hist(trace[, "beta0"], col = "azure1", prob = T, breaks = 50, xlab = expression(paste(beta[0])), main = "", ylab = "")
lines(x0, dnorm(x0, mean = mode1[1], sd = sqrt(cov1[1,1]), log = F), lty = 2, col = "black", lwd = 2)
lines(x0,dnorm(x0, mean = improved_mean[1], sd = sqrt(cov1[1,1]), log = F), lty = 1, col = "blue", lwd = 2)
lines(x0, dt(x0, 0, df = 5), lty = 4, lwd = 2, col = "red")

hist(trace[, "beta1"], col = "azure1", prob = T, breaks = 50, xlab = expression(paste(beta[1])), main = "", ylab = "")
lines(x0, dnorm(x0, mean = mode1[2], sd = sqrt(cov1[2,2]), log = F), lty = 2, col = "black", lwd = 2)
lines(x0,dnorm(x0, mean = improved_mean[2], sd = sqrt(cov1[2,2]), log = F), lty = 1, col = "blue", lwd = 2)
lines(x0, dunif(x0, min = -3, max = 3), lty = 4, lwd = 2, col = "red")

res= rbind(t(m), t(improved_mean), t(mean.mcmc))
colnames(res) = c("beta0", "beta1")
rownames(res) = c("gaussian", "vb.corrected", "mcmc")
print(res)

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


