library(sn)
library(mvtnorm)

skews <- c(-0.8, -0.5, -0.2, 0, 0.2, 0.5, 0.8)

p = 2
mu = rep(0,p)
Sigma0 = matrix(c(1,0.5,0.5,1), byrow = T, ncol = p)

Corr0 = cov2cor(Sigma0)
chol0 = chol(Corr0)

for (ii in 1:length(skews)){
  for (jj in 1:length(skews)){
    s0 = c(skews[ii],skews[jj])
  
    
    params1 = INLA:::inla.sn.reparam(moments = c(mu[1], Sigma0[1,1], s0[1]))
    params2 = INLA:::inla.sn.reparam(moments = c(mu[2], Sigma0[2,2], s0[2]))
    
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
    
    x_g <- seq(-3, 3, by = 0.1)
    y_g = x_g
    z_sgc_T = matrix(rep(NA, length(x_g)*length(y_g)), nrow = length(y_g), ncol = length(x_g))
    
    
    par0 = list(params1,
                params2)
    
    
    for (i in 1:length(x_g)){
      for (j in 1:length(y_g)){
        z_sgc_T[i,j] = d_SGC(c(x_g[i], y_g[j]), Sigma0, list(params1, params2))
      }
    }
    
    par(mar = c(2,2,1,1))
    plot(x_g, y_g, type = "n")
    contour(x_g, y_g, z_sgc_T, nlevels = 10, col = "black", add = T, drawlabels = F, lty = 1, lwd = 1)
    points(x = 0, y = 0, pch = 16, cex = 1, lwd = 2)
    
  }
}
