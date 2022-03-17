library(INLA)
library(MASS)
data(Aids2)

Aids_data <- Aids2
Aids_data$time = Aids_data$death - Aids_data$diag
Aids_data$C = 0
Aids_data$C[Aids_data$status=="D"] = 1
Aids_data$time = Aids_data$time/1000
Aids_data$AZT = 0
Aids_data$AZT[Aids_data$diag>10043] = 1

expanded_df = inla.coxph(inla.surv(time = time, event = C) ~ 1 + AZT + f(age, model = "rw1") + T.categ,
                 data = Aids_data, 
                 control.hazard = list(model="rw1", n.intervals = 50))

#Old INLA
res1 = inla(expanded_df$formula,
            data = c(as.list(expanded_df$data),expanded_df$data.list),
            family = expanded_df$family,
            E = expanded_df$E)
tOld = res1$cpu.used[4]

#New INLA
res2 = inla(expanded_df$formula,
            data = c(as.list(expanded_df$data),expanded_df$data.list),
            family = expanded_df$family,
            E = expanded_df$E,
            inla.mode = "experimental")
summary(res2)

par(mar = c(4,4,1,1))

plot(res2$summary.random$baseline.hazard[,1]*1000,
      res2$summary.random$baseline.hazard[,2],
     type = "l",
     xlab = "Time",
     ylab = "Log Baseline hazard",
     ylim = c(-1.5,1.5),
     lwd = 2)
lines(res2$summary.random$baseline.hazard[,1]*1000,
               res2$summary.random$baseline.hazard[,6],
     type = "l",
     add = T,
     col = "blue",
     lwd = 2, lty = 2)
lines(res2$summary.random$baseline.hazard[,1]*1000,
               res2$summary.random$baseline.hazard[,4],
     type = "l",
     add = T,
     col = "blue",
     lwd = 2, lty = 2)

plot(res2$summary.random$age$ID,
     exp(res2$summary.random$age$mean),
     type = "l",
     ylab = "Hazard ratio",
     xlab = "Age",
     ylim = c(0,3),
     lwd = 2)
lines(res2$summary.random$age$ID,
     exp(res2$summary.random$age$`0.025quant`) ,
     lty = 2, lwd = 2, col = "blue")
lines(res2$summary.random$age$ID,
      exp(res2$summary.random$age$`0.975quant`) ,
      lty = 2, lwd = 2, col = "blue")
abline(h = 1, lty = 3)
abline(v = 49, lty = 3)

exp(res2$summary.fixed$mean)

tNew = res2$cpu.used[4]





