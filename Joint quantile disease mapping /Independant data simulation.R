# The code of the simulated independent data 

library(pander)
library(SpatialEpi)
library(BayesX)
library(spdep)
library(INLA)
library(mvtnorm)
map <- pennLC$spatial.polygon

nb <- poly2nb(map)
nb2INLA("map.adj", nb)
g1 <- inla.read.graph(filename = "map.adj")
adj <- poly2nb(map)
R <- nb2mat(adj , style = "B", zero.policy = TRUE)
num_n <- colSums(R)
Q <- -R
diag(Q) <- num_n + 1
tau_x = 1
Q = tau_x*Q
n=67
n.rep = 10


# The simulated data

sigma1 = solve(Q)

x1 = c()
for (i in 1:n.rep){
  x1 = c(x1,rmvnorm(1, sigma= sigma1))
}

# The cases of disease 1
n1 <- n*n.rep
eta1 <- 1 + x1
alpha1 <- 0.5
shape1 <- exp(eta1)+1
lambda1 <- qgamma(alpha1, shape = shape1, lower.tail = FALSE, rate = 1)
y1 <- rpois(n1, lambda1)
y1




sigma1 = solve(Q)

x2 = c()
for (i in 1:n.rep){
  x2 = c(x2,rmvnorm(1, sigma= sigma1))
}

# The cases of disease 2

n2 <- n*n.rep
eta2 <- 1 + x2
alpha2 <- 0.5
shape2 <- exp(eta2)+1
lambda2 <- qgamma(alpha2, shape = shape2, lower.tail = FALSE, rate = 1)
y2 <- rpois(n2, lambda2)
y2




# # Seperate model of disease 1 

u = c(1:67)
alpha = alpha1
rep_1= c(rep(1:n.rep, each = n))

rsm = inla(formula = y1 ~ 1+
             f(u,model = "besagproper", scale.model = T,
               graph = g1, replicate = rep_1)   ,
           family = "poisson",
           data =  data.frame(u , y1),
           control.predictor = list(compute = T),
           control.family = list(control.link = list(model = "quantile",
                                                     quantile = alpha)),
           control.compute = list(dic = T, waic = T, cpo = T))

# Separate model of disease  2


u = c(1:67)
alpha = alpha1
rep_2= c(rep(1:n.rep, each = n))
rsm2 = inla(formula = y2 ~ 1+
              f(u,model = "besagproper", scale.model = T,
                graph = g1, replicate = rep_2)   ,
            family = "poisson",
            data =  data.frame(u , y2),
            control.predictor = list(compute = T),
            control.family = list(control.link = list(model = "quantile",
                                                      quantile = alpha)),
            control.compute = list(dic = T, waic = T, cpo = T))




# Joint quantile model

b = length(y1)
y1a = c(y1 , rep(NA,b))
y2a = c(rep(NA,b) , y2)
yy = cbind(y1a,y2a)
mu = c(rep(1, b) , rep(2,b))
Besag = c(rep(1:n, n.rep) , rep(NA ,b))
Besag.c = c(rep(NA , b), rep(1:n, n.rep) )
rep_i1 = c(rep(1:n.rep, each = n), rep(NA, b))
rep_i2 = c(rep(NA, b), rep(1:n.rep, each = n))

u = Besag
v = Besag.c


m = as.factor(mu)
d = data.frame(yy,m , Besag , Besag.c, rep_i1, rep_i2, u, v)  

formula = yy ~ -1 + m + f(u,model = "besagproper",
                          graph = g1, replicate = rep_i1) + 
  f(v,model = "besagproper",
    graph = g1, replicate = rep_i2) +
  f(Besag, model = "besagproper",graph=g1, replicate = rep_i1)+
  f(Besag.c, copy="Besag",
    hyper = list(beta = list(fixed = FALSE, params = c(0,0.001))), replicate = rep_i2)

r1 <- inla(formula,
           family = c("poisson","poisson"),
           control.family = list(list(control.link = list(model = "quantile",quantile = alpha1)),
                                 
                                 list(control.link = list(model = "quantile",
                                                          quantile = alpha2))), 
           data = d,
           verbose = F,
           control.predictor = list(compute = T),
           control.compute = list(dic = T, waic = T, cpo = T),
           inla.mode = "experimental")
summary(r1)



dt = data.frame(
  DIC = c(rsm$dic$dic, rsm2$dic$dic, rsm$dic$dic + rsm2$dic$dic  ,r1$dic$dic),
  WAIC = c(rsm$waic$waic, rsm2$waic$waic, rsm$waic$waic+rsm2$waic$waic  ,r1$waic$waic))

rownames(dt) <- c("Separate 1", "Separate 2", "Sum of Separates "  ,"Joint qunatile")
pander(dt)