####Example 1

##MCMC - define
cat("model
    {

    for ( i in 1:N ) {
    u[i] ~ dnorm(0, u_sd)
    lambda[i] = exp(beta0 + beta1*x[i] + beta2*x1[i] + beta3*x2[i] + u[i])
    y[i] ~ dpois( lambda[i])

    }

    ### Define the priors
    beta0 ~ dnorm( 0, 1 )
    beta1 ~ dnorm( 0, 1 )
     beta2 ~ dnorm( 0, 1 )
      beta3 ~ dnorm( 0, 1 )
      u_sd ~ dgamma(1,1)

    }", file="sim_jags_VB.txt")
cat("model
    {

    for ( i in 1:N ) {
    lambda[i] = exp(beta0 + beta1*x[i] )
    y[i] ~ dpois( lambda[i])

    }

    ### Define the priors
    beta0 ~ dnorm( 0, 1 )
    beta1 ~ dnorm( 0, 1 )
     

    }", file="sim_jags_VB.txt")
inits <- list(beta0 = 0, beta1 = 0)
params <- c("beta0", "beta1")

##VB - define
Model <- function(parm, Data)
{
  ### Parameters
  beta <- parm[Data$pos.beta]
  ### Log-Prior
  beta.prior <- sum(dnormv(beta, 0, 1, log=TRUE))
  ### Log-Likelihood
  mu <- exp(tcrossprod(Data$X, t(beta)))
  LL <- sum(dpois(Data$y, mu, log=TRUE))
  ### Log-Posterior
  LP <- LL + beta.prior 
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=mu[1],
                   yhat=rnorm(length(mu), mu), parm=parm)
  return(Modelout)
}


###Data

if(TRUE){
set.seed(0)
N = 10000
b0 = -1
x = rnorm(N, mean = 0, sd = 1)
b1 = -0.5
eta = b0 + b1*x 
Y = rpois(n = N, lambda = exp(eta))
prec = 1

par(mfrow = c(1,1))
par(mar = c(4,2,2,1))
barplot(table(as.factor(Y)), main = "", ylim = c(0, N))

#Gaussian
  
INLA_result = inla(formula = Y ~ 1 + x, data=data.frame(Y = Y, x = x), 
                   family = "poisson", 
                   control.fixed = list(prec = prec,   
                                       prec.intercept = prec,
                                       correlation.matrix=TRUE),
                   control.compute=list(config = TRUE),
                   control.inla = list(strategy = "gaussian",
                                       control.vb = list(enable = FALSE)))
Gtime = system.time(inla(formula = Y ~ 1 + x  , data=data.frame(Y = Y, x = x), 
                               family = "poisson", 
                               control.fixed = list(prec = prec,   
                                                    prec.intercept = prec,
                                                    correlation.matrix=TRUE),
                               control.compute=list(config = TRUE),
                         control.inla = list(strategy = "gaussian",
                                             control.vb = list(enable = FALSE))))

#Laplace
INLA_result_L = inla(formula = Y ~ 1 + x  , data=data.frame(Y = Y, x = x), 
                     family = "poisson", 
                     control.compute=list(config = TRUE),
                     control.fixed = list(prec = prec,   
                                          prec.intercept = prec,
                                          correlation.matrix=TRUE), 
                     control.inla = list(strategy = "laplace"))

Ltime = system.time(inla(formula = Y ~ 1 + x  , data=data.frame(Y = Y, x = x), 
                         family = "poisson", 
                         control.compute=list(config = TRUE),
                         control.fixed = list(prec = prec,   
                                              prec.intercept = prec,
                                              correlation.matrix=TRUE), 
                         control.inla = list(strategy = "laplace")))

#Gaussian with VB
INLA_result_VB = inla(formula = Y ~ 1 + x , data=data.frame(Y = Y, x = x), 
                   family = "poisson", 
                   control.compute=list(config = TRUE),
                   control.fixed = list(prec = prec,   
                                        prec.intercept = prec,
                                        correlation.matrix=TRUE), 
                   control.inla = list(strategy = "gaussian",
                                       control.vb = list(enable = TRUE)))
VBCtime = system.time(inla(formula = Y ~ 1 + x , data=data.frame(Y = Y, x = x), 
                          family = "poisson", 
                          control.compute=list(config = TRUE),
                          control.fixed = list(prec = prec,   
                                               prec.intercept = prec,
                                               correlation.matrix=TRUE), 
                          control.inla = list(strategy = "gaussian",
                                              control.vb = list(enable = TRUE))))

# summary(INLA_result_VB)
# lines(INLA_result_VB$marginals.fixed$`(Intercept)`, type="l", lty = 3, col = "blue")
# plot(INLA_result_VB$marginals.fixed$x, type="l")
# par(mar = c(2,2,2,2))
# barplot(table(Y), col = "lightgrey", ylim = c(0,60), main = "", xlab = "")

#MCMC
library(rjags)
dat1 = list(y = Y, x = x, N = N)
jags.m <- jags.model( file = "sim_jags_VB.txt", data=dat1, inits=inits, n.chains=1, n.adapt=500 )
samps <- coda.samples( jags.m, params, n.iter=10^5, start = 10001 )
MCtime = system.time(coda.samples( jags.m, params, n.iter=10^5, start = 10001 ))

#VB
# library(LaplacesDemon)
# ##############################  Data from above  ###############################
# y <- Y
# X <- cbind(1, as.matrix(x), as.matrix(x1), as.matrix(x2))
# J <- ncol(X)
# 
# #########################  Data List Preparation for VB  #########################
# mon.names <- "mu[1]"
# parm.names <- as.parm.names(list(beta=rep(0,J)))
# pos.beta <- grep("beta", parm.names)
# PGF <- function(Data) {
#   beta <- rnorm(Data$J)
#   
#   return(c(beta))
# }
# MyData <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
#                parm.names=parm.names, pos.beta=pos.beta, y=y)
# Initial.Values <- rep(0,J)

#Fit <- VariationalBayes(Model, Initial.Values, Data=MyData, Covar=NULL,
#     Iterations=5000, Method="Salimans2", Stop.Tolerance=1e-3, CPUs=1)
#VBtime = system.time(VariationalBayes(Model, Initial.Values, Data=MyData, Covar=NULL,
#                                      Iterations=1000, Method="Salimans2", Stop.Tolerance=1e-3, CPUs=1))

#Results
results <- data.frame(Name = c("Gaussian beta0", "Gaussian beta1", "MCMC beta0", "MCMC beta1", "Laplace beta0", "Laplace beta1",  "VBC beta0", "VBC beta1"),
                      Post_mean = c(INLA_result$summary.fixed$mean,list(summary(samps))[[1]][1]$statistics[1:2],INLA_result_L$summary.fixed$mean,INLA_result_VB$summary.fixed$mean))
results

times <- data.frame(Name = c("Gaussian", "MCMC", "Laplace", "VBC"),
                    time = c(Gtime[3], MCtime[3], Ltime[3], VBCtime[3]))
times
}


######
dev.off()
par(mar = c(4,2,1,1))
hist(samps[[1]][,1], breaks = seq(-1.2, -0.8, length = 100), prob = T, col = "azure1",ylim = c(0,25), ylab = "", xlab = expression(paste(beta[0])), main = "", xlim = c(-1.15, -0.9))
lines(INLA_result$marginals.fixed$`(Intercept)`, xlim = c(-2, 0), ylim = c(0,3), col = "black", lty = 2, ylab = "", xlab = expression(paste(beta[0])), lwd = 2)
#lines(INLA_result_L$marginals.fixed$`(Intercept)`, type="l", lty = 2, col = "red", lwd = 2)
lines(INLA_result_VB$marginals.fixed$`(Intercept)`, type="l", lty = 1, col = "blue", lwd = 2)
abline(v = b0, lwd = 3)

hist(samps[[1]][,2], breaks = seq(-0.7, -0.3, length = 100), prob = T, col = "azure1", , ylim = c(0,25), ylab = "", xlab = expression(paste(beta[1])), main = "", xlim = c(-0.6, -0.4))
lines(INLA_result$marginals.fixed$x, col = "black", xlim = c(-2, 0), ylim = c(0, 5), lty = 2, lwd = 2, ylab = "", xlab = expression(paste(beta[1])))
#lines(INLA_result_L$marginals.fixed$x, type="l", lty = 2, col = "red", lwd = 2)
lines(INLA_result_VB$marginals.fixed$x, type="l", lty = 1, col = "blue", lwd = 2)
abline(v = b1, lwd = 3)





#MLE with VB
mle <- glm(formula = Y ~ x , 
           data = data.frame(Y = Y, x = x),
           family = "poisson")
summary(mle)

G <- list(n = n, y = y, x = x)
#Create minimization function to correct MLE (i.e. uniform prior)
n = N
G$n = n
G$N = length(mle$coefficients) + n
G$y = Y

#creat matrix for premultiplication
xx = rbind(cbind(rep(1,G$n), x), diag(G$N - G$n))
Spost = xx %*% vcov(mle) %*% t(xx) + 1E-10*diag(G$N)
Spost[(G$n+1):G$N, (G$n+1):G$N]
Qpost = solve(Spost)
Qpost = Qpost + t(Qpost)
diag(Qpost) <- diag(Qpost)/2
G$post <- list(S = Spost, Q = Qpost, mu = c(log(mle$fitted.values),mle$coefficients), var = diag(Spost))

Qprior = G$post$Q
diag(Qprior) <- diag(Qprior) - c(exp(G$post$mu[1:G$n]),rep(0, G$N-G$n)) - 1E-6
Qprior[(G$n+1):G$N, (G$n+1):G$N]
G$prior <- list(S = solve(Qprior), Q = Qprior, mu = rep(0, G$N), var = diag(solve(Qprior)))

vb.fun <-function(delta)
{mEll <-function(y, mu, sd) 
{xx <-seq(mu-7*sd, mu+7*sd, by = 0.05*sd)
dd <-dnorm(xx, mean = mu, sd = sd)
ll <--dpois(y, lambda =exp(xx), log = TRUE)
ss <-sum(dd*ll)/sum(dd)
return(ss)}

kld <- function(prior, post) {
  val <- c(1/2 * (trace(prior$Q %*% solve(post$Q)) +
                    c(post$mu %*% prior$Q %*% post$mu) -dim(post$Q)[1] - ldet(prior$Q) + ldet(post$Q)))
  return(val) }


d.mu <-numeric(G$N)

for(j in seq_along(delta))
{d.mu <- d.mu+delta[j]*G$post$S[G$n+j, ]}

mmu <- G$post$mu+d.mu
value <- 0 

for(i in 1:G$n) 
{value <- value+ mEll(G$y[i], mmu[i],sqrt(G$post$var[i]))}
value <- value + kld(G$prior, list(Q = G$post$Q, mu = mmu))

G$sol <<- list(mu = mmu[(G$n+1):G$N], sd = sqrt(G$post$var[(G$n+1):G$N]))
return(value)
}


delta <- rep(0.1, G$N-G$n)
res <- optim(delta, vb.fun, NULL, method = "BFGS")

Sol <- cbind(mle.vbcor = c(G$sol$mu, G$sol$sd),
             mle = c(mle$coefficients,
                     sqrt(diag(vcov(mle)))), 
             r.la = c(INLA_result_L$summary.fixed[, "mean"],
                      INLA_result_L$summary.fixed[, "sd"]),
             r.ga = c(INLA_result$summary.fixed[, "mean"],
                      INLA_result$summary.fixed[, "sd"]),
             ga.vbcor = c(INLA_result_VB$summary.fixed[, "mean"],
                          INLA_result_VB$summary.fixed[, "sd"]))

m <- dim(Sol)[1] %/% 2
rownames(Sol) <- c(paste0("mean", 1:m), paste0("sd", 1:m)) 
print(round(dig = 4, Sol))


#####Tokyo example
data(Tokyo)
formula = y ~ -1 + f(time, model="rw2", cyclic=TRUE,
                     ## fix the precision here so you can see the effect better
                     initial = -4,
                     constr = FALSE,
                     scale.model = TRUE,
                     fixed = TRUE,
                     vb.correct = TRUE) 

r = inla(formula,
         family="binomial",
         Ntrials=n,
         data=Tokyo,
         control.inla = list(strategy = "gaussian"))

rr = inla(formula,
          family="binomial",
          Ntrials=n,
          data=Tokyo,
          control.inla = list(strategy = "gaussian", 
                              control.vb = list(enable = TRUE)), 
          verbose = TRUE)

rrr = inla(formula,
           family="binomial",
           Ntrials=n, data=Tokyo,
           control.inla = list(strategy = "laplace"))

summary(r)
summary(rr)
summary(rrr)

abs_errorG = sum(abs(r$summary.random$time$mean-rrr$summary.random$time$mean))
abs_errorVB = sum(abs(rr$summary.random$time$mean-rrr$summary.random$time$mean))

par(mar = c(2,2,1,1))
plot(r$summary.random$time$mean[1:60], col = "black", pch = 19, cex = 0.3, ylim = c(-2.0,-1.4))
lines(rr$summary.random$time$mean[1:60], lwd = 2, col = "blue")
lines(rrr$summary.random$time$mean[1:60], lwd = 4, lty = 2, col = "red")

plot(r$summary.random$time$mean, col = "black", pch = 19, cex = 0.3)
lines(rr$summary.random$time$mean, lwd = 2, col = "blue")
lines(rrr$summary.random$time$mean, lwd = 4, lty = 2, col = "red")

plot(r$marginals.random$time[[339]], col = "black", pch = 19, cex = 0.3, xlim = c(-4, -1))
abline(v = r$summary.random$time$mean[339], lwd = 2)
lines(rr$marginals.random$time[[339]], lwd = 2, col = "blue")
abline(v = rr$summary.random$time$mean[339], col = "blue", lwd = 2)
lines(rrr$marginals.random$time[[339]], lwd = 2, lty = 2, col = "red")
abline(v = rrr$summary.random$time$mean[339], col = "red", lty = 2, lwd = 4)



######Real example
setwd("~/")
PRSA_d1 <- read.csv("Downloads/PRSA_Data_20130301-20170228/PRSA_Data_Aotizhongxin_20130301-20170228.csv")
PRSA_d2 <- read.csv("Downloads/PRSA_Data_20130301-20170228/PRSA_Data_Changping_20130301-20170228.csv")
PRSA_d3 <- read.csv("Downloads/PRSA_Data_20130301-20170228/PRSA_Data_Dingling_20130301-20170228.csv")
PRSA_d4 <- read.csv("Downloads/PRSA_Data_20130301-20170228/PRSA_Data_Dongsi_20130301-20170228.csv")
PRSA_d5 <- read.csv("Downloads/PRSA_Data_20130301-20170228/PRSA_Data_Guanyuan_20130301-20170228.csv")
PRSA_d6 <- read.csv("Downloads/PRSA_Data_20130301-20170228/PRSA_Data_Gucheng_20130301-20170228.csv")
PRSA_d7 <- read.csv("Downloads/PRSA_Data_20130301-20170228/PRSA_Data_Huairou_20130301-20170228.csv")
PRSA_d8 <- read.csv("Downloads/PRSA_Data_20130301-20170228/PRSA_Data_Nongzhanguan_20130301-20170228.csv")
PRSA_d9 <- read.csv("Downloads/PRSA_Data_20130301-20170228/PRSA_Data_Shunyi_20130301-20170228.csv")
PRSA_d10 <- read.csv("Downloads/PRSA_Data_20130301-20170228/PRSA_Data_Tiantan_20130301-20170228.csv")
PRSA_d11 <- read.csv("Downloads/PRSA_Data_20130301-20170228/PRSA_Data_Wanliu_20130301-20170228.csv")
PRSA_d12 <- read.csv("Downloads/PRSA_Data_20130301-20170228/PRSA_Data_Wanshouxigong_20130301-20170228.csv")

PM_data = rbind(PRSA_d1, PRSA_d2, PRSA_d3, PRSA_d4, PRSA_d5,
                PRSA_d6, PRSA_d7, PRSA_d8, PRSA_d9, PRSA_d10,
                PRSA_d11, PRSA_d12)

saveRDS(PM_data, file = "Downloads/PM_data.RDS")
PM_data <- readRDS("~/Data/PM_data.RDS")
PM_data$f_month <- PM_data$month
PM_data$f_station <- PM_data$station
PM_data$f_day <- PM_data$day
PM_data$SO2 <- scale(PM_data$SO2)
PM_data$NO2 <- scale(PM_data$NO2)
PM_data$CO <- scale(PM_data$CO)
PM_data$O3 <- scale(PM_data$O3)
PM_data$lPM2.5 <- log(PM_data$PM2.5+0.0001)
PM_data$Ind <- rep(0, nrow(PM_data))
PM_data$Ind[PM_data$lPM2.5 > log(25)] <- 1

PM_data = PM_data[which(PM_data$f_station == "Guanyuan"),]

hist(PM_data$lPM2.5)

res <- inla(lPM2.5 ~ 1 + SO2 + f(No, model = "ar1"),
            family = "gaussian",
            data = PM_data,
            control.inla = list(strategy = "gaussian",
                                control.vb = list(enable = FALSE)),
            control.compute = list(config = TRUE),
            inla.mode = "experimental",
            control.predictor = list(compute = TRUE))

resVB <- inla(lPM2.5 ~ 1 + SO2 + CO + NO2 + f(No, model = "ar1") + f(f_station, model = "iid"),
            family = "gaussian",
            data = PM_data,
            control.inla = list(strategy = "gaussian",
                                control.vb = list(enable = TRUE)),
            inla.mode = "experimental",
            control.predictor = list(compute = TRUE))

resL <- inla(lPM2.5 ~ 1 + SO2 + f(inla.group(TEMP), model = "rw2", scale.model = TRUE),
            family = "gaussian",
            data = PM_data,
            control.inla = list(strategy = "laplace"),
            inla.mode = "experimental")

saveRDS(res, file = "~/Data/PMGau.RDS")
saveRDS(resVB, file = "~/Data/PMVB.RDS")
saveRDS(resL, file = "~/Data/PML.RDS")


resI <- inla(Ind ~ 1 + SO2 + NO2 + CO + O3 + f(f_month, model = "rw2", scale.model = TRUE) +
              f(f_station, model = "iid"),
            family = "binomial",
            data = PM_data,
            control.predictor = list(compute = TRUE),
            control.inla = list(strategy = "gaussian",
                                control.vb = list(enable = TRUE)),
            inla.mode = "experimental")

res1I <- inla(Ind ~ 1 + SO2 + NO2 + CO + O3 + f(f_month, model = "rw2", scale.model = TRUE) +
              f(f_station, model = "iid"),
            family = "binomial",
            data = PM_data,
            control.predictor = list(compute = TRUE),
            control.inla = list(strategy = "gaussian",
                                control.vb = list(enable = FALSE)),
            inla.mode = "experimental")

res2I <- inla(Ind ~ 1 + SO2 + NO2 + CO + O3 + f(f_month, model = "rw2", scale.model = TRUE) +
              f(f_station, model = "iid"),
            family = "binomial",
            data = PM_data,
            control.predictor = list(compute = TRUE),
            control.inla = list(strategy = "laplace"),
            inla.mode = "experimental")

summary(abs(PM_data$lPM2.5-res$summary.fitted.values$mean))
hist(PM_data$lPM2.5)
hist(res$summary.fitted.values$mean, add = TRUE, col="blue", transparency = 0.5)
plot(res$summary.random$f_month$ID, res$summary.random$f_month$mean, type = "l")


##Leuk
data(Leuk)
data1 <- Leuk

nwseg <- inla.sp2segment(nwEngland)
bnd1 <- inla.nonconvex.hull(nwseg$loc, 0.03, 0.1, resol = 50)
#bnd2 <- inla.nonconvex.hull(nwseg$loc, 0.25)
os=365 #Yearly
#lines(bnd2)
data1$time = data1$time/os
data1$age = scale(data1$age)
data1$tpi = scale(data1$tpi)
data1$wbc = scale(data1$wbc)

c <- inla.coxph(inla.surv(time, cens) ~ 0 + a0 + sex + age + wbc + tpi,
                data = data.frame(a0 = 1, data1),
                control.hazard = list(n.intervals = 100, model = "rw2"))
coo = cbind(c$data$xcoord,c$data$ycoord)

#r=10
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


#Scaled covariates- good idea
r <- inla(f1, 
           family = "poisson",
           E = data2$E,
           data = data2, 
           control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
          control.inla = list(strategy = "gaussian",
                               control.vb = list(enable = FALSE)),
          inla.mode = "experimental") 
summary(r)

r1 <- inla(f1, 
           family = "poisson",
           E = data2$E,
           data = data2, 
           control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
           control.inla = list(strategy = "gaussian",
                               control.vb = list(enable = TRUE)),
           inla.mode = "experimental") 
summary(r1)
rr <- inla(f1, 
           family = "poisson",
           E = data2$E,
           data = data2, 
           control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
           control.inla = list(strategy = "laplace"), inla.mode = "experimental") 
summary(rr)
par(mfrow = c(1,1),mar = c(4,3,0.5,0.1))
plot(r$summary.random$baseline.hazard$ID*os, r$summary.random$baseline.hazard$mean, pch = 16, cex = 0.5, ylim = c(-3,3), xlab = "Time", ylab = "Baseline hazard")
points(r$summary.random$baseline.hazard$ID*os, r$summary.random$baseline.hazard$`0.025quant` , pch = 16, cex = 0.3, ylim = c(-2,2), xlab = "Time", ylab = "")
points(r$summary.random$baseline.hazard$ID*os, r$summary.random$baseline.hazard$`0.975quant`, pch = 16, cex = 0.3, ylim = c(-2,2), xlab = "Time", ylab = "")

lines(r$summary.random$baseline.hazard$ID*os, r1$summary.random$baseline.hazard$mean, type = "l", lty = 1, col = "blue", lwd = 2)
lines(r$summary.random$baseline.hazard$ID*os, r1$summary.random$baseline.hazard$`0.025quant`, type = "l", lty = 1, col = "blue", lwd = 1)
lines(r$summary.random$baseline.hazard$ID*os, r1$summary.random$baseline.hazard$`0.975quant`, type = "l", lty = 1, col = "blue", lwd = 1)

lines(r$summary.random$baseline.hazard$ID*os, rr$summary.random$baseline.hazard$mean, type = "l", lty = 2, col = "red", lwd = 3)
lines(r$summary.random$baseline.hazard$ID*os, rr$summary.random$baseline.hazard$`0.025quant`, type = "l", lty = 2, col = "red", lwd = 2)
lines(r$summary.random$baseline.hazard$ID*os, rr$summary.random$baseline.hazard$`0.975quant`, type = "l", lty = 2, col = "red", lwd = 2)

plot(r$summary.random$baseline.hazard$ID[1:20]*os, r$summary.random$baseline.hazard$mean[1:20], pch = 16, cex = 0.5, ylim = c(-0,3), xlab = "Time", ylab = "")
points(r$summary.random$baseline.hazard$ID[1:20]*os, r$summary.random$baseline.hazard$`0.025quant`[1:20] , pch = 16, cex = 0.3, ylim = c(-2,2), xlab = "Time", ylab = "")
points(r$summary.random$baseline.hazard$ID[1:20]*os, r$summary.random$baseline.hazard$`0.975quant`[1:20], pch = 16, cex = 0.3, ylim = c(-2,2), xlab = "Time", ylab = "")

lines(r$summary.random$baseline.hazard$ID[1:20]*os, r1$summary.random$baseline.hazard$mean[1:20], type = "l", lty = 1, col = "blue", lwd = 2)
lines(r$summary.random$baseline.hazard$ID[1:20]*os, r1$summary.random$baseline.hazard$`0.025quant`[1:20], type = "l", lty = 1, col = "blue", lwd = 1)
lines(r$summary.random$baseline.hazard$ID[1:20]*os, r1$summary.random$baseline.hazard$`0.975quant`[1:20], type = "l", lty = 1, col = "blue", lwd = 1)

lines(r$summary.random$baseline.hazard$ID[1:20]*os, rr$summary.random$baseline.hazard$mean[1:20], type = "l", lty = 2, col = "red", lwd = 3)
lines(r$summary.random$baseline.hazard$ID[1:20]*os, rr$summary.random$baseline.hazard$`0.025quant`[1:20], type = "l", lty = 2, col = "red", lwd = 2)
lines(r$summary.random$baseline.hazard$ID[1:20]*os, rr$summary.random$baseline.hazard$`0.975quant`[1:20], type = "l", lty = 2, col = "red", lwd = 2)


###Results
setwd("~/Documents/GitHub/Variational-Bayes/VB_mean paper")
PM_G <- readRDS("~/Documents/GitHub/Variational-Bayes/VB_mean paper/PMGaus.RDS")
PM_VB <- readRDS("~/Documents/GitHub/Variational-Bayes/VB_mean paper/PMVB.RDS")
PM_L <- readRDS("~/Documents/GitHub/Variational-Bayes/VB_mean paper/PMLaplace.RDS")



##Leuk results
Leuk_G <- readRDS("~/Documents/GitHub/Variational-Bayes/VB_mean paper/leukGaus.RDS")
Leuk_VB <- readRDS("~/Documents/GitHub/Variational-Bayes/VB_mean paper/leukVB.RDS")
Leuk_L <- readRDS("~/Documents/GitHub/Variational-Bayes/VB_mean paper/leukLaplace.RDS")
Leuk_G <- r
Leuk_VB <- r1
Leuk_L <- rr

par(mfrow = c(1,1),mar = c(4,2,0.5,0.1))
plot(Leuk_G$marginals.fixed$a0, col = "black", pch = 19, cex = 0.3, ylab = "", xlab = expression(paste(beta[0])))
abline(v = Leuk_G$summary.fixed$mean[1], lwd = 2)
lines(Leuk_VB$marginals.fixed$a0, col = "blue", lwd = 2)
abline(v = Leuk_VB$summary.fixed$mean[1], lwd = 2, col = "blue")
lines(Leuk_L$marginals.fixed$a0, col = "red", lwd = 4, lty = 2)
abline(v = Leuk_L$summary.fixed$mean[1], lwd = 4, col = "red", lty = 2)

plot(Leuk_G$marginals.fixed$age , col = "black", pch = 19, cex = 0.3, ylab = "", xlab = expression(paste(beta[1])))
abline(v = Leuk_G$summary.fixed$mean[3], lwd = 2)
lines(Leuk_VB$marginals.fixed$age, col = "blue", lwd = 2)
abline(v = Leuk_VB$summary.fixed$mean[3], lwd = 2, col = "blue")
abline(v = Leuk_L$summary.fixed$mean[3], lwd = 4, col = "red", lty = 2)
lines(Leuk_L$marginals.fixed$age, col = "red", lwd = 4, lty = 2)

plot(Leuk_G$marginals.fixed$wbc  , col = "black", pch = 19, cex = 0.3, ylab = "", xlab = expression(paste(beta[2])))
lines(Leuk_VB$marginals.fixed$wbc, col = "blue", lwd = 2)
lines(Leuk_L$marginals.fixed$wbc, col = "red", lwd = 4, lty = 2)
abline(v = Leuk_G$summary.fixed$mean[4], lwd = 2)
abline(v = Leuk_VB$summary.fixed$mean[4], lwd = 2, col = "blue")
abline(v = Leuk_L$summary.fixed$mean[4], lwd = 4, col = "red", lty = 2)

plot(Leuk_G$marginals.fixed$tpi , col = "black", pch = 19, cex = 0.3, ylab = "", xlab = expression(paste(beta[3])))
lines(Leuk_VB$marginals.fixed$tpi, col = "blue", lwd = 2)
lines(Leuk_L$marginals.fixed$tpi, col = "red", lwd = 4, lty = 2)
abline(v = Leuk_G$summary.fixed$mean[5], lwd = 2)
abline(v = Leuk_VB$summary.fixed$mean[5], lwd = 2, col = "blue")
abline(v = Leuk_L$summary.fixed$mean[5], lwd = 4, col = "red", lty = 2)

plot(Leuk_G$marginals.hyperpar$`Precision for baseline.hazard` , col = "black",type = "l", ylab = "", xlab = expression(tau))
plot(Leuk_G$marginals.hyperpar$`Range for spatial` , col = "black", type= "l", ylab = "", xlab = "r")
plot(Leuk_G$marginals.hyperpar$`Stdev for spatial` , col = "black", type = "l", ylab = "", xlab = expression(paste(sigma[u])))



plot(mesh, asp = 1, col = "grey", lwd = 0.1, main = "")
points(coo, col='red', pch = 20, cex = 0.5 ) 
lines(nwEngland@polygons[[1]] , col = "black", pch = 19, cex = 0.1)

plot(nwEngland)

bbnw <- bbox(nwEngland)
bbnw
r0 <- diff(range(bbnw[1, ])) / diff(range(bbnw[2, ]))
r0
prj <- inla.mesh.projector(mesh, xlim = bbnw[1, ], 
                           ylim = bbnw[2, ], dims = round(c(r0, 1)*300))

## NA's were not to plot 
spat.m <- inla.mesh.project(prj, r$summary.random$spatial$mean)
spat.sd <- inla.mesh.project(prj, r$summary.random$spatial$sd)
spat.mvb <- inla.mesh.project(prj, r1$summary.random$spatial$mean)
spat.minla <- inla.mesh.project(prj, rr$summary.random$spatial$mean)
ov <- over(SpatialPoints(prj$lattice$loc), nwEngland)
spat.sd[is.na(ov)] <- NA
spat.m[is.na(ov)] <- NA
spat.mvb[is.na(ov)] <- NA
spat.minla[is.na(ov)] <- NA

### plot the spatial risk
library(fields)
library(raster)
par(mfrow = c(1,1), mar = c(0, 0, 0, 5))
plot(nwEngland)
image.plot(x = prj$x, y = prj$y, z = spat.m, add=TRUE)
plot(nwEngland)
image.plot(x = prj$x, y = prj$y, z = spat.mvb, add=TRUE)
plot(nwEngland)
image.plot(x = prj$x, y = prj$y, z = spat.minla, add=TRUE)

plot(nwEngland)
image.plot(x = prj$x, y = prj$y, z = abs(spat.mvb-spat.m), add=TRUE)

plot(nwEngland)
image.plot(x = prj$x, y = prj$y, z = spat.sd, add = TRUE)




