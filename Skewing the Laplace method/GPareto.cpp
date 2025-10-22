#include <numeric>
#include <math.h>
#include <algorithm>
#include <iostream>
//#include <Rcpp.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace R;
using namespace RcppArmadillo;
using namespace std;


/*** R
library(runjags)
library(INLA)
library(evd)

rgp = function(n, sigma, eta, alpha, xi = 0.001)
	{
	if (missing(sigma)) {
	stopifnot(!missing(eta) && !missing(alpha))
	sigma = exp(eta) * xi / ((1.0 - alpha)^(-xi) -1.0)
	}
	#set.seed(1)
	return (sigma / xi * (runif(n)^(-xi) - 1.0))
	}

n = 10
set.seed(6) 
x = runif(n) - 3
x = c(x, -1)
n = length(x)
eta = 1 + x
alpha = 0.95
xi = 0.75
y = rgp(n, eta = eta, alpha = alpha, xi=xi)
#y = c(y, 0.56)
#n = n + 1

xi.high <- 0.9
to.theta <- inla.models()$likelihood$gp$hyper$theta$to.theta
xi.intern <- to.theta(xi, interval = c(0, xi.high))

r = inla(y ~ 1 + x,
	inla.mode = "classic",
	data = data.frame(y, x), 
	family = "gp",  
	control.family = list(control.link = list(quantile = alpha), hyper = list(xi = list(initial = xi.intern, param = c(7, 0, xi.high), fixed = TRUE))),  
	control.fixed = list(prec = 0.001, prec.intercept = 0.001), 
	control.predictor = list(compute=TRUE), 
	control.compute = list(config = TRUE),
	verbose=FALSE)

tik1 <- Sys.time()
r.vb.mean = inla(y ~ 1 + x,
	inla.mode = "experimental",
	data = data.frame(y, x), 
	family = "gp",  
	control.family = list(control.link = list(quantile = alpha), hyper = list(xi = list(initial = xi.intern, param = c(7, 0, xi.high), fixed = TRUE))),  
	control.fixed = list(prec = 0.001, prec.intercept = 0.001), 
	control.predictor = list(compute=TRUE), 
	control.compute = list(config = TRUE),
	verbose=FALSE, 
	control.inla = list(control.vb = list(strategy = "mean", f.enable.limit = c(n, n))))
tok1 <- Sys.time()
print(tok1 - tik1)	

tik2 <- Sys.time()
r.vb.var = inla(y ~ 1 + x,
	inla.mode = "experimental",
	data = data.frame(y, x), 
	family = "gp",  
	control.family = list(control.link = list(quantile = alpha), hyper = list(xi = list(initial = xi.intern, param = c(7, 0, xi.high), fixed = TRUE))),  
	control.fixed = list(prec = 0.001, prec.intercept = 0.001), 
	control.predictor = list(compute=TRUE), 
	control.compute = list(config = TRUE),
	verbose=FALSE, 
	control.inla = list(control.vb = list(strategy = "variance", f.enable.limit = c(n, n))))
tok2 <- Sys.time()
print(tok2 - tik2)		
*/

// [[Rcpp::export]]
NumericVector dGenPareto(NumericVector x, NumericVector eta, double alpha, double xi, String logarithm){
	
	int number = x.size();
	NumericVector sigma(number);
	NumericVector z(number);
	NumericVector result(number);
		
	for (int i = 0; i < number; i += 1)
	{
		sigma[i] = exp(eta[i])* xi  / (pow((1.0 - alpha), -xi) - 1.0);
		//sigma[i] = pow((1.0 - alpha), -xi);
	}
	
	for (int i = 0; i < number; i += 1)
	{
		z[i] = x[i]/sigma[i];
	}
	
	for (int i = 0; i < number; i += 1)
	{
		result[i] = -log(sigma[i]) - ((1/xi + 1)*log(1 + xi*z[i]));
	}
	
	if (logarithm == "TRUE")
	{
		return(result);
	}else{
		for (int i = 0; i < number; i += 1)
		{
		result[i] = exp(result[i]);
		}
		return(result);
	}

	return(sigma);
}


// [[Rcpp::export]]
double lposterior(NumericVector theta, NumericVector x, NumericVector y, double alpha, double xi){
	
	int number = x.size();
	double intercept = theta[0];
	double beta = theta[1];
	
	NumericVector eta(number);
	for (int i = 0; i < number; i += 1)
	{
		eta[i] = intercept + beta*x[i];
	}
	
	
	NumericVector loglik = dGenPareto(y, eta, alpha, xi = xi, "TRUE");
	double logposterior = std::accumulate(loglik.begin(), loglik.end(), 0.0) + R::dnorm(intercept, 0, sqrt(1000), TRUE) + R::dnorm(beta, 0, sqrt(1000), TRUE);

	
	return(logposterior);
}


// [[Rcpp::export]]
List MetroPolis(long int niter, NumericVector theta, NumericVector x, NumericVector y, double alpha, double xi){
	
	NumericMatrix theta_trace(niter, 2);
	long int naccept = 0;
	double s = 1; 
	
	
	
	for (int iter = 0; iter < niter; iter += 1)
	{	
		NumericVector theta_prime;
		
		theta_prime = clone(theta);

		double idx = rbinom(1, 1, 0.5)[0];
		
		theta_prime[idx] = theta_prime[idx] + rnorm(1, 0, s)[0];
		
		double accept_prob = exp(min(0.0, lposterior(theta_prime, x, y, alpha, xi) - lposterior(theta = theta, x, y, alpha, xi)));
		
		double r = runif(1,0,1)[0];
		if (r < accept_prob)
		{
			theta = theta_prime;
			naccept = naccept + 1;
		}
		
		
		theta_trace(iter, 0) = theta[0];
		theta_trace(iter, 1) = theta[1];
		
	}
	
	
	List ret;
	ret["naccept"] = naccept;
	ret["theta_trace"] = theta_trace;
	
	
	
	return(ret);
}





/*** R
dGenPareto(y, eta, alpha, xi = 0.75, "TRUE")
dgpd(y, loc = 0, scale = 0.75*exp(eta)/(0.05^(-0.75) - 1), shape = 0.75, log = T)
dGenPareto(y, eta, alpha, xi = 0.75, "FALSE")
theta <- c(0, 0)

lposterior(theta = theta, x = x, y = y, alpha = alpha, xi = xi)

niter = 10000000
args = list(x = x, y = y)
tik3 <- Sys.time()
aa <- MetroPolis(niter, theta, x, y, alpha, xi)
tok3 <- Sys.time()
print(tok3 - tik3)
aa$naccept <- aa$naccept/niter;

aa$theta_trace <- aa$theta_trace[-(1:10000),]
niter <- nrow(aa$theta_trace)
print(niter)
print(aa$naccept)


#dev.new()
hist(aa$theta_trace[, 1], breaks = 150, prob = TRUE, main = "Intercept")
#lines(inla.smarginal(r$marginals.fixed$`(Intercept)`), col = "blue", lwd = 2)
lines(inla.smarginal(r.vb.mean$marginals.fixed$`(Intercept)`), col = "red", lwd = 2)
lines(inla.smarginal(r.vb.var$marginals.fixed$`(Intercept)`), col = "green", lwd = 2)

#dev.new()
hist(aa$theta_trace[, 2], n = 150, prob = TRUE, main = "Beta")
#lines(inla.smarginal(r$marginals.fixed$x), col = "blue", lwd = 2)
lines(inla.smarginal(r.vb.mean$marginals.fixed$x), col = "red", lwd = 2)
lines(inla.smarginal(r.vb.var$marginals.fixed$x), col = "green", lwd = 2)

eta1 <- seq(-5.5, -4.2, 0.05)
feta <- numeric(length(eta1))
#for(i in seq_along(eta)) feta[i] <- sum(dt( sqrt((df - 2)/df)*(y - eta[i]), df = df, log = T))
for(i in seq_along(eta)) feta[i] <- sum(dgpd(y, loc = 0, scale = 0.75*exp(eta1[i])/(0.05^(-0.75) - 1), shape = 0.75, log = T))
detat <- splinefun(eta1, feta)

r <- optim(0, detat, NULL, hessian = TRUE, method = "BFGS", control = list(fnscale = -1))

layout(matrix(c(1,1,2,2,3,3,4,4,4,5,5,5), nrow = 2, ncol = 6, byrow = TRUE))
plot(eta1, exp(feta), type = "l")
plot(eta1, feta, type = "l")

plot(eta1, detat(eta1, deriv = 1), type = "l")


xx <- seq(-25, 25, 0.01)
fxx <- numeric(length(xx))
for (i in seq_along(xx)) fxx[i] <- sum(dgpd(y, loc = 0, scale = 0.75*exp(xx[i])/(0.05^(-0.75) - 1), shape = 0.75, log = T))
detaf <- splinefun(xx, fxx)
r <- optim(-3, detaf, NULL, hessian = TRUE, method = "BFGS", control = list(fnscale = -1))
#r$hessian <- (sum(dgpd(y, loc = 0, scale = 0.75*exp(r$par + 1e-3)/(0.05^(-0.75) - 1), shape = 0.75, log = F))
#             - 2*sum(dgpd(y, loc = 0, scale = 0.75*exp(r$par)/(0.05^(-0.75) - 1), shape = 0.75, log = F))
#             + sum(dgpd(y, loc = 0, scale = 0.75*exp(r$par - 1e-3)/(0.05^(-0.75) - 1), shape = 0.75, log = F)))/1e-6

# library(pracma)
# first_moment <- trapz(xx, xx*fxx/sum(fxx*0.01))
# second_moment <- trapz(xx, xx*xx*fxx/sum(fxx*0.01))


layout(matrix(c(1,1,2,2,3,3,4,4,4,5,5,5), nrow = 2, ncol = 6, byrow = TRUE))
idx <- c(2051:2751)
plot(xx[idx], exp(fxx[idx]), type = "l", lwd =2, main = "(a)", ylab = " ", xlab = expression(eta))
#abline(v = r$par)
#lines(xx, dnorm(xx,  mean = r$par, sd = sqrt(3.368376)), type = "l", col = "red")

idx <- c(1:length(xx))
plot(xx[idx], (fxx[idx])- max((fxx[idx])), type = "l", lwd =2, main = "(b)", ylab = " ", xlab = expression(eta))
lines(xx[idx], dnorm(xx[idx], mean = r$par, sd = sqrt(-1/r$hessian), log = T) - max(dnorm(xx[idx], mean = r$par, sd = sqrt(-1/r$hessian), log = T)), col = "red", lty = 2, lwd = 2)

plot(xx, detaf(xx, deriv = 2), type = "l", lwd =2, main = "(c)", ylab = " ", xlab = expression(eta))
abline(h = 0, col = "red", lwd =2)
abline(h = 1/r$hessian, col = "blue", lwd = 2, lty = 2)

hist(aa$theta_trace[, 1], breaks = 150, prob = TRUE, main = "(d)", ylim = c(0, 0.2), xlim = c(-12, 6), xlab = expression(beta[0]), ylab = " ")
#lines(inla.smarginal(r$marginals.fixed$`(Intercept)`), col = "blue", lwd = 2)
lines(inla.smarginal(r.vb.mean$marginals.fixed$`(Intercept)`), col = "red", lwd = 2)
lines(inla.smarginal(r.vb.var$marginals.fixed$`(Intercept)`), col = "blue", lwd = 2, lty = 2)

#dev.new()
hist(aa$theta_trace[, 2], n = 150, prob = TRUE, main = "(e)", xlim = c(-5,3), ylim = c(0, 0.45), xlab = expression(beta[1]), ylab = " ")
#lines(inla.smarginal(r$marginals.fixed$x), col = "blue", lwd = 2)
lines(inla.smarginal(r.vb.mean$marginals.fixed$x), col = "red", lwd = 2)
lines(inla.smarginal(r.vb.var$marginals.fixed$x), col = "blue", lwd = 2, lty = 2)

r.vb.mean$summary.fixed[, "sd"]
r.vb.var$summary.fixed[, "sd"]
r.mcmc <- cbind(mean = apply(aa$theta_trace, 2, mean), sd = apply(aa$theta_trace, 2, sd))
*/











