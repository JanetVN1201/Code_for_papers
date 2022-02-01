library(icenReg)
library(INLA)
data(IR_diabetes)

IR_diabetes$event = rep(3, nrow(IR_diabetes))
for (i in 1:nrow(IR_diabetes)){
if (IR_diabetes$right[i]==IR_diabetes$left[i]){
  IR_diabetes$event[i] = 1
}
}
IR_diabetes$ID = 1:nrow(IR_diabetes)

res1 = inla(inla.surv(time = IR_diabetes$left, time2 = IR_diabetes$right,
                  event = IR_diabetes$event) ~ as.factor(gender) + f(ID, model = "iid"),
            data = IR_diabetes,
            family = "loglogistic.surv",
            verbose = FALSE,
            control.compute = list(config = TRUE),
            control.family = list(variant = 1))

summary(res1)


