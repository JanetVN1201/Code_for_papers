
#####Tokyo example in Section 5.1
data(Tokyo)
formula = y ~ -1 + f(time, model="rw2", cyclic=TRUE,
                     initial = -4,
                     constr = FALSE,
                     scale.model = TRUE,
                     fixed = TRUE,
                     vb.correct = TRUE) 
#GA
r = inla(formula,
         family="binomial",
         Ntrials=n,
         data=Tokyo,
         control.inla = list(strategy = "gaussian",
                             control.vb = list(enable = FALSE)))
#INLA-VBC
rr = inla(formula,
          family="binomial",
          Ntrials=n,
          data=Tokyo,
          control.inla = list(strategy = "gaussian", 
                              control.vb = list(enable = TRUE)), 
          verbose = TRUE)
#INLA
rrr = inla(formula,
           family="binomial",
           Ntrials=n, data=Tokyo,
           control.inla = list(strategy = "laplace"))

####### JAGS
library(runjags)
model = "model {
  for(i in 1:N){
    y[i] ~ dbinom(p[i],trials[i])
    logit(p[i]) = u[i]
  }
  
u[1] ~ dnorm(-1.6, tau)
u[2] ~ dnorm(-1.5, tau)
for (t in 3:(N-1)){
u[t] ~ dnorm(2*u[t-1] - u[t-2], tau)
}
u[N] ~ dnorm(-1.6, tau)
tau ~ dgamma(25,0.10)
}"

library(tictoc)
tic("MCMC")
res.1 = run.jags(model = model,
                 monitor = c("u", "tau"),
                 data = list('y' = Tokyo$y, 'trials' = Tokyo$n, 'N' = nrow(Tokyo)),
                 n.chains = 6,
                 inits = list(u = rep(0.5,nrow(Tokyo))),
                 burnin = 10^3,
                 sample = 10^4,
                 method = "parallel")

trace <- combine.mcmc(res.1)
toc()


##Leuk
data(Leuk)
data1 <- Leuk

nwseg <- inla.sp2segment(nwEngland)
bnd1 <- inla.nonconvex.hull(nwseg$loc, 0.03, 0.1, resol = 50)

os=365 #Yearly

data1$time = data1$time/os

#Scale the covariates
data1$age = scale(data1$age)
data1$tpi = scale(data1$tpi)
data1$wbc = scale(data1$wbc)


c <- inla.coxph(inla.surv(time, cens) ~ 0 + a0 + sex + age + wbc + tpi,
                data = data.frame(a0 = 1, data1),
                control.hazard = list(n.intervals = 100, model = "rw2"))
coo = cbind(c$data$xcoord,c$data$ycoord)


mesh <- inla.mesh.2d(coo, max.edge = c(0.05,0.1), boundary = bnd1, offset = c(0.1,0.1))
plot(mesh, asp = 1, col = "grey", lwd = 0.1, main = "")
points(coo, col='red', pch = 20, cex = 0.5 ) 
points(nwseg$loc, col = "black", pch = 19, cex = 0.1)


A <- inla.spde.make.A(mesh, coo)
spde <- inla.spde2.pcmatern(mesh = mesh,
                            prior.range = c(0.05, 0.5), # P(range < 0.05) = 0.01
                            prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01
stk <- inla.stack(
  data = c(list(E = c$E), c$data[c('y..coxph')]), 
  A = list(A, 1), 
  effect = list(
    list(spatial = 1:spde$n.spde), 
    c$data[c('baseline.hazard', 'a0', 
                    'age', 'sex', 'wbc', 'tpi')])) 

data2 <- c(inla.stack.data(stk), c$data.list)

f1 <- update(c$formula, '. ~ . + f(spatial, model = spde)')
#GA
r <- inla(f1, 
           family = "poisson",
           E = data2$E,
           data = data2, 
           control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
          control.inla = list(strategy = "gaussian",
                               control.vb = list(enable = FALSE)),
          inla.mode = "experimental") 
summary(r)

#INLA-VBC
r1 <- inla(f1, 
           family = "poisson",
           E = data2$E,
           data = data2, 
           control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
           control.inla = list(strategy = "gaussian",
                               control.vb = list(enable = TRUE)),
           inla.mode = "experimental") 
summary(r1)

#INLA
rr <- inla(f1, 
           family = "poisson",
           E = data2$E,
           data = data2, 
           control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
           control.inla = list(strategy = "laplace"), inla.mode = "experimental") 
summary(rr)



