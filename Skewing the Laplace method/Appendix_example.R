library(sn)
library(mvtnorm)

NN = 2000
p=2
mu = rep(0,p)
Sigma0 = matrix(c(1,0.5,0.5,1), byrow = T, ncol = p)
s0 = c(-0.8,-0.8)
Corr0 = cov2cor(Sigma0)
chol0 = chol(Corr0)

params1 = INLA:::inla.sn.reparam(moments = c(mu[1], Sigma0[1,1], s0[1]))
params2 = INLA:::inla.sn.reparam(moments = c(mu[2], Sigma0[2,2], s0[2]))

library(mvtnorm)
UU = cbind(rep(NA, NN), rep(NA, NN))
YY = UU
XX = rmvnorm(n = NN, mean = mu, sigma = Corr0)
for (i in 1:NN){
YY[i,1] = qsn(pnorm(XX[i,1], mean = 0, sd = 1), dp = unlist(params1))
YY[i,2] = qsn(pnorm(XX[i,2], mean = 0, sd = 1), dp = unlist(params2))
}

summary(YY)

x_c = cut(YY[,1], 30)
y_c = cut(YY[,2], 30)
z <- table(x_c, y_c)

# mode0_ind = which(z == max(z), arr.ind = TRUE)
# mode0 = c(levels(x_c)[mode0_ind[1,1]], levels(y_c)[mode0_ind[1,2]])

library(plot3D)
hist3D(z = z, border = "black")
#image2D(z = z, border = "black")
scatter2D(x = YY[,1], y = YY[,2], col = "lightgrey", pch = 16, cex = 0.5)

Sigma0_est = matrix(c(sd(YY[,1])^2, cov(YY[,1],YY[,2]), cov(YY[,1],YY[,2]),sd(YY[,2])^2),
byrow = T, ncol = 2)

mu_est = c(mean(YY[,1]), mean(YY[,2]))

delta = rep(0.5,2)

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
  
  return(d0*dmvnorm(x = as.vector(vec0), sigma = cov2cor(Sigma1))/d1)
  
  
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

opt_f <- function(delta){
 
  kk = nrow(y)
  E_ll = 0
  
  par0 = list(INLA:::inla.sn.reparam(moments = c(0, sigma0[1,1], delta[1])),
                                     INLA:::inla.sn.reparam(moments = c(0, sigma0[2,2], delta[2])))
  
    
  for (i in 1:kk){
    E_ll = E_ll + log(d_SGC(y[i,], sigma0, par0))
  }
  
  KLD0 = kld_f(mu0, sigma0, delta)
  
  return(-E_ll + KLD0) 
  
}

mu0 = mu_est
sigma0 = Sigma0_est

print(mu0)
print(sigma0)

# mu0 = mu
# sigma0 = Sigma0
y = YY

results <- optim(delta, opt_f, lower = c(-0.9,-0.9), upper = c(0.9, 0.9), method = "L-BFGS-B")

results$par

#Gaussian approx
x_g <- seq(-6, 6, by = 0.1)
y_g = x_g
z_g = matrix(rep(NA, length(x_g)*length(y_g)), nrow = length(y_g), ncol = length(x_g))
z_sgc = z_g
z_sgc_T = z_g

par0 = list(INLA:::inla.sn.reparam(moments = c(mu0[1], sigma0[1,1], results$par[1])),
            INLA:::inla.sn.reparam(moments = c(mu0[2], sigma0[2,2], results$par[2])))


for (i in 1:length(x_g)){
  for (j in 1:length(y_g)){
    z_g[i,j] = dmvnorm(cbind(x_g[i], y_g[j]), mean = mu0, sigma = sigma0)
    z_sgc[i,j] = d_SGC(c(x_g[i], y_g[j]), sigma0, par0)
    z_sgc_T[i,j] = d_SGC(c(x_g[i], y_g[j]), Sigma0, list(params1, params2))
  }
}

par(mar = c(0,0,0,0))
scatter2D(x = YY[,1], y = YY[,2], col = "lightgrey", pch = 16, cex = 0.5, xlim = c(-3,3), ylim = c(-3,3))
contour(x_g, y_g, z_g, nlevels = 10, add = T, lty = 2, lwd = 2, drawlabels = F)
contour(x_g, y_g, z_sgc, nlevels = 10, col = "blue", add = T, lwd = 2, drawlabels = F)
contour(x_g, y_g, z_sgc_T, nlevels = 10, col = "red", add = T, drawlabels = F, lty = 3, lwd = 2)
legend(legend = c("Laplace method", "SGC-VB", "True"), lwd = rep(2,3), lty = c(2, 1, 3), col = c("black","blue","red"),
       x = 1.5, y = -1.8)
