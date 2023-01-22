## utility function for quantile-regression.
##     - continous Poisson and quantile regression for Poisson

inla.incGamma = function(x, lambda=0, log=FALSE) {
  ##library(gsl); return (gamma_inc(x, lambda))
  lres = lgamma(x) + pgamma(lambda, x, lower = FALSE, log.p=TRUE)
  return (if (log) lres else exp(lres))
}

inla.qcontpoisson = function(p, lambda, print.level = 0) {
  fun.min = function(log.x, lambda, prob) {
    p.est = inla.pcontpoisson(exp(log.x), lambda)
    eps = 1E-12
    p.est = max(eps, min(p.est, 1 - eps))
    return ((inla.link.logit(p.est) - inla.link.logit(prob))^2)
  }
  initial.value = log(qpois(p, lambda) + 0.5)
  return (exp(nlm(fun.min, p = initial.value, print.level = print.level,
                  lambda = lambda, prob = p)$estimate))
}

inla.pcontpoisson = function(x, lambda, log=FALSE, deriv=0) {
  ## > F := (x, lambda) -> GAMMA(x, lambda) / GAMMA(x);
  ##                                     GAMMA(x, lambda)
  ##                F := (x, lambda) -> ----------------
  ##                                         GAMMA(x)
  ##
  ## > simplify(diff(F(x, lambda), lambda));           
  ##                             (x - 1)
  ##                       lambda        exp(-lambda)
  ##                     - --------------------------
  ##                                GAMMA(x)
  ##
  ## > simplify(diff(F(x, lambda), lambda$2));         
  ##                   (x - 2)
  ##             lambda        exp(-lambda) (-x + 1 + lambda)
  ##             --------------------------------------------
  ##                               GAMMA(x)
  ##
  if (deriv == 0) {
    ## can use gsl::gamma_inc_Q() instead
    lres = inla.incGamma(x, lambda, log=TRUE) - inla.incGamma(x, log=TRUE)
    return (if (log) lres else exp(lres))
  } else if (deriv == 1) {
    stopifnot(!log)
    return (-exp((x-1)*log(lambda) -lambda -inla.incGamma(x, log=TRUE)))
  } else if (deriv == 2) {
    stopifnot(!log)
    return ((lambda-x+1) * exp((x-2)*log(lambda) -lambda -inla.incGamma(x, log=TRUE)))
  } else {
    stop("deriv != 0, 1, 2")
  }
}

inla.pcontpoisson.eta = function(x, eta, deriv = 0, log=FALSE) {
  ## the cdf of the contpoisson parameterised by the linear predictor and the log-link
  lambda = exp(eta)
  if (deriv == 0) {
    return (inla.pcontpoisson(x, lambda, log=log))
  } else if (deriv == 1) {
    stopifnot(!log)
    return (inla.pcontpoisson(x, lambda, deriv=1) * lambda)
  } else if (deriv == 2) {
    stopifnot(!log)
    return (lambda * (inla.pcontpoisson(x, lambda, deriv=1) +
                        inla.pcontpoisson(x, lambda, deriv=2) * lambda))
  } else {
    stop("deriv != 0, 1, 2")
  }
}

inla.contpoisson.solve.lambda = function(quantile, alpha, iter.max = 1000, max.step = 3,
                                         tol = sqrt(.Machine$double.eps), verbose=FALSE)
{
  ## solve quantile=inla.pcontpoisson(lambda, alpha),  for lambda
  stopifnot(length(quantile) == 1 && length(alpha) == 1)
  return(exp(inla.contpoisson.solve.eta(quantile, alpha, iter.max, max.step, tol, verbose)))
}

inla.contpoisson.solve.eta = function(quantile, alpha, iter.max = 1000, max.step = 3,
                                      tol = sqrt(.Machine$double.eps), verbose=FALSE) 
{
  ## solve quantile=inla.pcontpoisson(lambda=exp(eta), alpha),  for eta
  stopifnot(length(quantile) == 1 && length(alpha) == 1)
  eta.0 = log((sqrt(quantile) - qnorm(alpha,  sd = 0.5))^2)
  for(i in 1:iter.max) {
    f = inla.pcontpoisson(quantile, lambda = exp(eta.0)) - alpha
    fd = inla.pcontpoisson.eta(quantile, eta.0, deriv=1)
    d = -min(max.step, max(-max.step, f/fd))
    eta = eta.0 = eta.0 + d
    if (verbose)
      print(round(c(iter=i, eta = eta, f=f, f.deriv = fd, err = d), digits = 6))
    if (abs(d) < tol)
      return(eta)
  }
  stop(paste("FAILURE", quantile, alpha))
}



## Continuous poisson functions
library(INLA)
library(bayesQR)
source("cont-poisson.R")



foo.fun1 = function(lambda.grid){
  
  nt = length(lambda.grid)
  colo = viridis::viridis(nt, option = "A", direction = 1, begin = .3, end = .6)
  lty1 = 1:length(lambda.grid)
  
  par(bty = "n")  
  curve(inla.pcontpoisson(x, lambda.grid[1]), from = -1, 10, xlab = "y", 
        ylab = expression(F[lambda](y)), lwd = 2, col = colo[1], ylim = c(0,1), lty = lty1[1])
  curve(ppois(x, lambda.grid[1]), add = T, col = colo[1], lwd = 2, lty = lty1[1])
  
  for(k in 2:nt){
    curve(inla.pcontpoisson(x, lambda.grid[k]), from = -1, 10, xlab = "y", 
          ylab = expression(F[lambda](y)), lwd = 2, col = colo[k], add = T, lty = lty1[k])
    curve(ppois(x, lambda.grid[k]), add = T, col = colo[k], lwd = 2, lty = 3, lty1[k])
    
  }
  legend("topleft", legend = lambda.grid, col = colo, lwd = 2, bty = "n", title = expression(lambda), lty = lty1)
  
}

Qcp  = Vectorize(inla.qcontpoisson)
foo.fun2 = function(tau.grid, xmax, lambda.grid){
  
  nt = length(tau.grid)
  colo = viridis::viridis(nt+1, option = "A", direction = -1, begin = .5, end = .9)[1:nt]
  qp = qpois(.5, lambda.grid)
  pp = ppois(qp, lambda.grid)
  lty1 = 1:length(lambda.grid)
  
  
  par(bty = "n")  
  plot(x=lambda.grid, Qcp(tau.grid[1], lambda.grid), type="l", 
       ylim= c(-1, xmax), col = colo[1], lwd = 2, 
       xlab = expression(lambda), ylab = expression(q[alpha](lambda)), lty = 1)
  lines(x=lambda.grid, qpois(tau.grid[1], lambda.grid), col = colo[1], lty = 1, lwd = 2)
  
  for(k in 2:nt){
    lines(x=lambda.grid, Qcp(tau.grid[k], lambda.grid), type="l", ylim= c(-1, xmax), col = colo[k], lwd = 2, lty = lty1[k])
    lines(x=lambda.grid, qpois(tau.grid[k], lambda.grid), col = colo[k], lwd = 2, lty = lty1[k])
    
  }
  legend("topleft", legend = tau.grid, col = colo[1:nt], lwd = 2, bty = "n", title = expression(alpha), lty = lty1)
}

pdf(file = "poisson.pdf", width = 8, height = 10)


par(mfrow= c(1,1), mar = c(4,4,0,0))
lam.grid = c( 1, 2.5, 5, 7, 10)
foo.fun1(lam.grid) 
lam.grid = seq(0.001, 6, length.out = 1000)
foo.fun2(c(0.05, .25, .5, .75, 0.95), 8, lam.grid) 
par(mfrow=c(1,1))
dev.off()


####Quantile crossing investigation


gen.crossing = function(seed, n, a=1 , b=0.5 ,m = 10, add.z = F ){
  set.seed(seed)
  x = rnorm(n)
  eta = a + b*x
  lambda = exp(eta)
  y = numeric(n)
  
  y = rpois(n, lambda = lambda)
  df = data.frame(x = x, y = y)
  
  
  
  col = viridis::viridis(m, direction = -1, option = "C")
  lty1 = 1:m
  res.list = list()
  alpha.seq = seq(.2, .9, length.out = m)
  
  
  par(mfrow=c(1,1))
  plot.quantiles = function(x, obj, col, add.plot = T){
    log.out = obj$summary.fixed[1,1] + obj$summary.fixed[2,1]*x
    return(exp(log.out))
    
    if(add.plot){
      lines(exp(log.out), col = col)
    }
  }
  
  
  plot.quantiles.machado = function(x, obj){
    log.out = obj$coefficients[1] + obj$coefficients[2]*x
    return(exp(log.out))
    
  }
  
  
  par(bty= "n")
  plot(df$x, y, pch = 20, col = rgb(0,0,0,.5), main = "",  xlab = "x", ylab = "y")#, xlim = c(2,3), ylim = c(20, 45))
  for(alpha in 1:length(alpha.seq)){
    
    
    temp = inla(y ~ 1 + x, data = df, family = "poisson", 
                control.family = list(control.link = 
                                        list(model = "quantile", quantile = alpha.seq[alpha])),
                control.predictor=list(compute=TRUE),
                verbose = TRUE)
    res.list[[alpha]] = temp
    curve(plot.quantiles(x, temp, col = "red"), add = T, col = col[alpha], lty = lty1[alpha], lwd =2)
    
    
  }
  
  
  
  par(bty= "n")
  plot(df$x, y, pch = 20, col = rgb(0,0,0,.5), main = "", xlab = "x", ylab = "y")
  machado.list = list()
  for(alpha in 1:length(alpha.seq)){
    temp = Qtools::rq.counts(y ~ 1 + x, data = df, tau = alpha.seq[alpha])
    machado.list[[alpha]] = temp
    curve(plot.quantiles.machado(x, temp), add = T, col = col[alpha], lty = lty1[alpha], lwd = 2)
    
  }
  
  
  
  if(add.z){res.list2 = list()
  df2 = df
  df2$z = rnorm(n, sd = 4)
  par(bty= "n")
  plot(df$x, df$y, pch = 20, col = rgb(0,0,0,.5), main = "The Wrong Model")
  for(alpha in 1:length(alpha.seq)){
    
    temp = inla(y ~ 1 + x + z, data = df2, family = "poisson", 
                control.family = list(control.link = 
                                        list(model = "quantile", quantile = alpha.seq[alpha])),
                control.predictor=list(compute=TRUE))
    res.list2[[alpha]] = temp
    curve(plot.quantiles(x, temp, col = "red"), add = T, col = col[alpha])
    
  }
  
  return(list(res.list, machado.list, res.list2, df))}
  else   return(list(res.list, machado.list, df))
  par(mfrow=c(1,1))
}

xx = gen.crossing(10, 70)
