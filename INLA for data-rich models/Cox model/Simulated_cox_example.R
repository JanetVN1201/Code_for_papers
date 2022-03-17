library(simsurv)
library(INLA)
haz <- function(t, x, betas, ...) {
betas[["shape"]] * (t ^ (betas[["shape"]] - 1)) *sin(betas[["alpha"]]*t/60) * exp(
    betas[["betaEvent_cov"]] * x[["Z1"]])
}

set.seed(1454)
N <- 100    # number of individuals
# Population (fixed effect) parameters
betas <- data.frame(
  shape                = rep(1.2, N),
  betaEvent_cov        = rep(0.1,N),
  alpha = rep(0.5,N)
  )
# Covariate data
covdat <- data.frame(
  Z1 = stats::rnorm(N, 0, 1)   # a continuous covariate
)
# Then we simulate the survival times based on the
# hazard function, covariates, and true parameter values
y <- simsurv(hazard = haz, x = covdat, betas = betas, maxt = 30)
y$x = covdat$Z1
hist(y$eventtime)
#Augmented data
expanded_df = inla.coxph(inla.surv(time = eventtime, event = status) ~ -1 + x,
                 data = y, 
                 control.hazard = list(model="rw2", n.intervals = 50))


res1 = inla(expanded_df$formula,
            data = c(as.list(expanded_df$data),expanded_df$data.list),
            family = expanded_df$family,
            E = expanded_df$E)

res2 = inla(expanded_df$formula,
            data = c(as.list(expanded_df$data),expanded_df$data.list),
            family = expanded_df$family,
            E = expanded_df$E,
            inla.mode = "experimental")

summary(res1)
summary(res2)

#Plot baseline hazard
plot(res1$summary.random$baseline.hazard[1:50,1], 
     exp(res1$summary.random$baseline.hazard[1:50,2]),
     lty = 1, lwd = 2, col = "black", type = "l",
     main = "Baseline hazard estimation")
lines(res2$summary.random$baseline.hazard[1:50,1], 
     exp(res2$summary.random$baseline.hazard[1:50,2]),
     lty = 2, lwd = 2, col = "blue", type = "l")
#Truth
lines(res1$summary.random$baseline.hazard[,1],
      betas$shape[1] * (res1$summary.random$baseline.hazard[,1] ^ (betas$shape[1] - 1))*
        sin(betas$alpha[1]*res1$summary.random$baseline.hazard[,1]/60),
      col = "red", lwd = 3, lty = 3)
legend(col = c("black", "blue", "red"),
       lty = c(1,2,3),
       lwd = c(2,2,2),
       legend = c("CoxPH family","Poisson family", "Truth"),
       x = 4, y = 0.5)
       


