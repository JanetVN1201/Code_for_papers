library(icenReg)
library(INLA)
data(IR_diabetes)
fit <- ic_par(cbind(left, right) ~ gender, 
              data = IR_diabetes,
              model = "po",
              dist = "loglogistic")
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
inla.rerun(res1)
summary(res1)

sigma_mar = inla.tmarginal(function(x) sqrt(1/x), res1$marginals.hyperpar$`Precision for ID`)
sigma_sum = inla.zmarginal(sigma_mar)

par(mar = c(4,4,1,1))
plot(res1$marginals.fixed$`(Intercept)`, type = "l", lwd = 2,main = "",
     ylab = expression(paste(pi,"(", beta[0], "|y)")), xlab = expression(beta[0]))

plot(res1$marginals.fixed$`as.factor(gender)male`, type = "l", lwd = 2,main = "",
     ylab = expression(paste(pi,"(", beta[1], "|y)")), xlab = expression(beta[1]))

plot(sigma_mar, type = "l", lwd = 2,main = "",
     ylab = expression(paste(pi,"(", sigma, "|y)")), xlab = expression(sigma))

plot(res1$marginals.hyperpar$`alpha for loglogistic observations`, type = "l", lwd = 2,main = "",
     ylab = expression(paste(pi,"(", alpha, "|y)")), xlab = expression(alpha))

sample1 = inla.posterior.sample(res1, n = 100)
IR_diabetes = cbind(IR_diabetes, matrix(rep(NA,73100), byrow = TRUE, nrow = 731))
hyper_dia = rep(NA,100)
betas = matrix(rep(NA, 200), ncol = 2)
for (i in 1:100){
IR_diabetes[,i+3] = sample1[[i]]$latent[1:731]
hyper_dia[i] = sample1[[i]]$hyperpar
betas[i,1] = sample1[[i]]$latent[1463]
betas[i,2] = sample1[[i]]$latent[1464]
}


tseq = seq(0,50, length = 200)
st_male = matrix(rep(NA, 200*100), ncol = 100)
st_female = st_male

for (i in 1:100){
  st_male[,i] = 1/(1+(exp(betas[i,1] + betas[i,2])*tseq)^hyper_dia[i])
  st_female[,i] = 1/(1+(exp(betas[i,1])*tseq)^hyper_dia[i])
}

st_m = matrix(rep(NA, 3*200), ncol = 3)
st_f = matrix(rep(NA, 3*200), ncol = 3)

for (i in 1:200){
st_m[i,1] = quantile(st_male[i,], 0.025)
st_m[i,2] = quantile(st_male[i,], 0.5)
st_m[i,3] = quantile(st_male[i,], 0.975)
st_f[i,1] = quantile(st_female[i,], 0.025)
st_f[i,2] = quantile(st_female[i,], 0.5)
st_f[i,3] = quantile(st_female[i,], 0.975)
}

par(mar = c(5,5,1,1))
npmleFit <- ic_np(cbind(left, right) ~ gender, data = IR_diabetes)
plot(npmleFit, main = "", col = c("red", "blue"), lty = c(1,2))



lines(tseq, st_m[,1], col = "blue", lty = 3)
lines(tseq, st_m[,2], col = "blue", lty = 2, lwd = 2)
lines(tseq, st_m[,3], col = "blue", lty = 3)
lines(tseq, st_f[,1], col = "red", lty = 3)
lines(tseq, st_f[,2], col = "red", lty = 1, lwd = 2)
lines(tseq, st_f[,3], col = "red", lty = 3)
legend(x = 30, y = 1, legend = c("Female", "Male"), lty = c(1,2), col = c("red","blue"))

#median st
tseq[which.min(abs(st_m[,2]-0.5))]
tseq[which.min(abs(st_f[,2]-0.5))]

###epileptic
#SANAD - epileptic############################
library(joineR)
library(INLA)

mTime=max(max(epileptic$time),max(epileptic$with.time))
epileptic$time<-epileptic$time/mTime
epileptic$with.time <-epileptic$with.time/mTime

epileptic$interaction <- with(epileptic, time * (treat == "LTG"))
epileptic$interaction2 <- with(epileptic, time * (treat == "CBZ"))
epileptic$interaction[epileptic$interaction==0]<-NA
epileptic$interaction2[epileptic$interaction2==0]<-NA

data1<-epileptic
dataL<-data.frame(ID=data1$id,LDose=log(data1$dose+0.1),Time=data1$time,TimeSp=data1$time,Age=data1$age,Gender=data1$gender,LD=data1$learn.dis,Treatment=data1$treat,InteractionLTG=data1$interaction,InteractionCBZ=data1$interaction2,list(V=data1$id,W=data1$id),Dose=data1$dose)
dataS<-data.frame(ID=data1$id,Time=as.numeric(data1$with.time),Status=data1$with.status2,StatusISC=data1$with.status.isc,StatusUAE=data1$with.status.uae,Age=data1$age,Gender=data1$gender,LD=data1$learn.dis,Treatment=data1$treat)
dataS<-subset(dataS,!duplicated(dataS$ID))
summary(dataL)
summary(dataS)

#Visualize
plot(dataL$Time[dataL$ID==42],dataL$Dose[dataL$ID==42],type="n",xlim=c(0,1),ylim=c(0,10))
for (i in unique(dataL$ID)){
  lines(dataL$Time[dataL$ID==i],dataL$Dose[dataL$ID==i],col=i)
}

dataLCBZ=dataL[dataL$Treatment=="CBZ",]
dataLLTG=dataL[dataL$Treatment=="LTG",]

#Mean trajectories
modCBZ<-inla(formula=dataLCBZ$Dose~f(inla.group(dataLCBZ$Time,n=50),model="rw2",scale.model = TRUE,hyper = list(prec = list(prior="pc.prec", param=c(1, 0.01)))),data=dataLCBZ,family="gaussian",control.compute=list(dic=TRUE,config=TRUE))
modLTG<-inla(formula=dataLLTG$Dose~f(inla.group(dataLLTG$Time,n=50),model="rw2",scale.model = TRUE,hyper = list(prec = list(prior="pc.prec", param=c(1, 0.01)))),data=dataLLTG,family="gaussian",control.compute=list(dic=TRUE,config=TRUE))
#Plot means
plot(modCBZ$summary.random$`inla.group(dataLCBZ|S|Time, n = 50)`[,1]*2400,log(modCBZ$summary.fixed[1,1]+modCBZ$summary.random$`inla.group(dataLCBZ|S|Time, n = 50)`[,2]),type="l",lwd=2,ylim=c(0.5,1.2),xlab="Time",ylab="log(Dose)",xlim=c(0,2400))
lines(modLTG$summary.random$`inla.group(dataLLTG|S|Time, n = 50)`[,1]*2400,log(modLTG$summary.fixed[1,1]+modLTG$summary.random$`inla.group(dataLLTG|S|Time, n = 50)`[,2]),lty=2,lwd=2)
legend(x=2000,y=0.7,legend=c("CBZ","LTG"),lty=c(1,2),lwd=c(2,2))

#Joint model
#Pre-work for entire predictor
nL<-nrow(dataL)
nS<-nrow(dataS)


dataS$int = rep(1, nrow(dataS))
dataS$time.l = dataS$Time
dataS$time.r = dataS$Time

for (i in dataS$ID){
  if (max(dataL$Time[dataL$ID==i])<dataS$Time[dataS$ID==i])
  {dataS$int[dataS$ID==i] = 3
  dataS$time.l[dataS$ID==i] = max(dataL$Time[dataL$ID==i])}
}

#not work it should be 3 only if it is 1
dataS$d.uae = dataS$StatusUAE*dataS$int
dataS$d.isc = dataS$StatusISC*dataS$int

fixed.eff<-data.frame(mu=as.factor(c(rep(1,nL),rep(1,nL),rep(2,nS),rep(3,nS))),
                      ageL=c(dataL$Age,dataL$Age,rep(0,nS),rep(0,nS)),
                      ageUAE=c(rep(0,nL),rep(0,nL),dataS$Age,rep(0,nS)),
                      ageISC=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$Age),
                      treatmentL=as.factor(c(dataL$Treatment,dataL$Treatment,rep(NA,nS),rep(NA,nS))),
                      treatmentUAE=as.factor(c(rep(NA,nL),rep(NA,nL),dataS$Treatment,rep(NA,nS))),
                      treatmentISC=as.factor(c(rep(NA,nL),rep(NA,nL),rep(NA,nS),dataS$Treatment)),
                      genderL=as.factor(c(dataL$Gender ,dataL$Gender ,rep(NA,nS),rep(NA,nS))),
                      genderUAE=as.factor(c(rep(NA,nL),rep(NA,nL),dataS$Gender ,rep(NA,nS))),
                      genderISC=as.factor(c(rep(NA,nL),rep(NA,nL),rep(NA,nS),dataS$Gender)))
random.eff<-list(timeL_LTG=c(dataL$InteractionLTG,rep(NA,nL),rep(NA,nS),rep(NA,nS)),
                 timeL_CBZ=c(dataL$InteractionCBZ,rep(NA,nL),rep(NA,nS),rep(NA,nS)),
                 linpredL=c(rep(NA,nL),dataL$ID,rep(NA,nS),rep(NA,nS)),
                 linpredL2=c(rep(NA,nL),rep(-1,nL),rep(NA,nS),rep(NA,nS)),
                 betaUAE=c(rep(NA,nL),rep(NA,nL),dataS$ID,rep(NA,nS)),
                 betaISC=c(rep(NA,nL),rep(NA,nL),rep(NA,nS),dataS$ID),
                 frailtyUAE=c(rep(NA,nL),rep(NA,nL),dataS$ID,rep(NA,nS)),
                 frailtyISC=c(rep(NA,nL),rep(NA,nL),rep(NA,nS),dataS$ID))


jointdata<-c(fixed.eff,random.eff)
y.long <- c(dataL$Dose,rep(NA,nL),rep(NA, nS),rep(NA,nS))
y.eta<-c(rep(NA,nL),rep(0,nL),rep(NA,nS),rep(NA,nS))
y.survUAE <- inla.surv(time = c(rep(NA, nL), rep(NA,nL),dataS$time.l,rep(NA,nS)), 
                       event = c(rep(NA, nL),rep(NA,nL),dataS$d.uae,rep(NA,nS)),
                       time2 = c(rep(NA, nL), rep(NA,nL),dataS$time.r,rep(NA,nS)))
y.survISC <- inla.surv(time = c(rep(NA, nL), rep(NA,nL),rep(NA,nS),dataS$time.l), 
                       event = c(rep(NA, nL),rep(NA,nL),rep(NA,nS),dataS$d.isc),
                       time2 = c(rep(NA, nL), rep(NA,nL),rep(NA,nS),dataS$time.r))
y.joint<-list(y.long,y.eta,y.survUAE,y.survISC)

jointdata$Y=y.joint

#Model fit
formula.model=Y~-1 + mu + treatmentL+treatmentUAE+treatmentISC+f(inla.group(timeL_LTG,n=50),model="rw2",scale.model = TRUE,hyper = list(prec = list(prior="pc.prec", param=c(1, 0.01))))+
  f(inla.group(timeL_CBZ,n=50),model="rw2",scale.model = TRUE,hyper = list(prec = list(prior="pc.prec", param=c(1, 0.01))))+
  f(linpredL, linpredL2, model="iid", hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(betaUAE, copy="linpredL", hyper = list(beta = list(fixed = FALSE)))+
  f(betaISC, copy="linpredL", hyper = list(beta = list(fixed = FALSE)))


Jointmodel= inla(formula.model, family = c("gaussian","gaussian","weibullsurv","weibullsurv"),
                 data = jointdata, verbose=TRUE, control.compute=list(config = TRUE, dic=TRUE),
                 control.family = list(
                   list(),
                   list(hyper = list(prec = list(initial = 10, fixed=TRUE))),
                   list(variant = 1),
                   list(variant = 1)
                 ), inla.mode = "experimental"
)

summary(Jointmodel)

#Plot splines

plot(Jointmodel$summary.random$`inla.group(timeL_LTG, n = 50)`[,1]*mTime,
     Jointmodel$summary.random$`inla.group(timeL_LTG, n = 50)`[,2], type = "l", 
     ylim= c(-1.1,0.5), lwd = 3, xlab = "Time", ylab = "Y(t)")
lines(Jointmodel$summary.random$`inla.group(timeL_LTG, n = 50)`[,1]*mTime,
      Jointmodel$summary.random$`inla.group(timeL_LTG, n = 50)`[,4], lty = 1, lwd = 1)
lines(Jointmodel$summary.random$`inla.group(timeL_LTG, n = 50)`[,1]*mTime,
      Jointmodel$summary.random$`inla.group(timeL_LTG, n = 50)`[,6], lty = 1, lwd = 1)
lines(Jointmodel$summary.random$`inla.group(timeL_CBZ, n = 50)`[,1]*mTime,
      Jointmodel$summary.random$`inla.group(timeL_CBZ, n = 50)`[,2], lty = 2, lwd = 3)
lines(Jointmodel$summary.random$`inla.group(timeL_CBZ, n = 50)`[,1]*mTime,
      Jointmodel$summary.random$`inla.group(timeL_CBZ, n = 50)`[,4], lty = 2, lwd = 1)
lines(Jointmodel$summary.random$`inla.group(timeL_CBZ, n = 50)`[,1]*mTime,
      Jointmodel$summary.random$`inla.group(timeL_CBZ, n = 50)`[,6], lty = 2, lwd = 1)
legend(legend = c("LTG","CBZ"), lty = c(1,2), x = 1500, y = -0.6)

#Survival curves
sample1 = inla.posterior.sample(Jointmodel, n = 100)
betas = matrix(NA, nrow = 100, ncol = 6)
hyperpar = matrix(NA, nrow = 100, ncol = 4)

linpreds = matrix(NA, ncol = 100, nrow = 2*nS)
for (i in 1:100){
  betas[i,1] = sample1[[i]]$latent[8717]
  betas[i,2] = sample1[[i]]$latent[8718]
  betas[i,3] = sample1[[i]]$latent[8719]
  betas[i,4] = sample1[[i]]$latent[8720]
  betas[i,5] = sample1[[i]]$latent[8721]
  betas[i,6] = sample1[[i]]$latent[8722]
  hyperpar[i,1] = sample1[[i]]$hyperpar[2]
  hyperpar[i,2] = sample1[[i]]$hyperpar[3]
  hyperpar[i,3] = sample1[[i]]$hyperpar[6]
  hyperpar[i,4] = sample1[[i]]$hyperpar[7]
}

tseq_LTG = Jointmodel$summary.random$`inla.group(timeL_LTG, n = 50)`[,1]
tseq_CBZ = Jointmodel$summary.random$`inla.group(timeL_CBZ, n = 50)`[,1]
st_uae_ltg = matrix(rep(NA, 50*100), ncol = 100)
st_uae_cbz = matrix(rep(NA, 47*100), ncol = 100, nrow = 47)
st_isc_ltg = matrix(rep(NA, 50*100), ncol = 100)
st_isc_cbz = matrix(rep(NA, 47*100), ncol = 100)

for (i in 1:100){
  spline_ltg = sample1[[i]]$latent[6805:6854]
  spline_cbz = sample1[[i]]$latent[6855:6901]
  st_uae_ltg[,i] = exp(-(exp(betas[i,2]+betas[i,5]+hyperpar[i,3]*(spline_ltg + betas[i,1]+betas[i,4]))*
                           tseq_LTG)^hyperpar[i,1])
  st_uae_cbz[,i] = exp(-(exp(betas[i,2]+hyperpar[i,3]*(spline_cbz + betas[i,1]))*
                           tseq_CBZ)^hyperpar[i,1])
  
  
  st_isc_ltg[,i] = exp(-(exp(betas[i,3]+ betas[i,6]+hyperpar[i,4]*(spline_ltg + betas[i,1]+betas[i,4]))*
                           tseq_LTG)^hyperpar[i,2])
  st_isc_cbz[,i] = exp(-(exp(betas[i,3]+hyperpar[i,4]*(spline_cbz + betas[i,1]))*
                           tseq_CBZ)^hyperpar[i,2])
  
}

st_uae_ltg1 = matrix(rep(NA, 3*50), ncol = 3)
st_uae_cbz1 = matrix(rep(NA, 3*47), ncol = 3)

st_isc_ltg1 = matrix(rep(NA, 3*50), ncol = 3)
st_isc_cbz1 = matrix(rep(NA, 3*47), ncol = 3)

for (i in 1:50){
  st_uae_ltg1[i,1] = quantile(st_uae_ltg[i,], 0.025)
  st_uae_ltg1[i,2] = quantile(st_uae_ltg[i,], 0.5)
  st_uae_ltg1[i,3] = quantile(st_uae_ltg[i,], 0.975)
  st_isc_ltg1[i,1] = quantile(st_isc_ltg[i,], 0.025)
  st_isc_ltg1[i,2] = quantile(st_isc_ltg[i,], 0.5)
  st_isc_ltg1[i,3] = quantile(st_isc_ltg[i,], 0.975)
}
for (i in 1:47){
  st_uae_cbz1[i,1] = quantile(st_uae_cbz[i,], 0.025)
  st_uae_cbz1[i,2] = quantile(st_uae_cbz[i,], 0.5)
  st_uae_cbz1[i,3] = quantile(st_uae_cbz[i,], 0.975)
  st_isc_cbz1[i,1] = quantile(st_isc_cbz[i,], 0.025)
  st_isc_cbz1[i,2] = quantile(st_isc_cbz[i,], 0.5)
  st_isc_cbz1[i,3] = quantile(st_isc_cbz[i,], 0.975)
}

par(mar = c(4,4,0,0))

plot(tseq_LTG*mTime, st_uae_ltg1[,1],  lty = 1, type = "l",
     ylab = expression(paste(S [UAE] (t))),
     xlab = "Time", ylim = c(0.2,1))
lines(tseq_LTG*mTime, st_uae_ltg1[,2],  lty = 1, lwd = 3)
lines(tseq_LTG*mTime, st_uae_ltg1[,3], lty = 1)
lines(tseq_CBZ*mTime, st_uae_cbz1[,1],  lty = 2)
lines(tseq_CBZ*mTime, st_uae_cbz1[,2], lty = 2, lwd = 3)
lines(tseq_CBZ*mTime, st_uae_cbz1[,3],  lty = 2)
legend(x = 1800, y = 1, legend = c("LTG", "CBZ"), lty = c(1,2), lwd = c(2,2))


plot(tseq_LTG*mTime, st_isc_ltg1[,1], lty = 1, type = "l", ylab = expression(paste(S [ISC] (t))),
     xlab = "Time", ylim = c(0.2,1))
lines(tseq_LTG*mTime, st_isc_ltg1[,2], lty = 1, lwd = 3)
lines(tseq_LTG*mTime, st_isc_ltg1[,3], lty = 1)
lines(tseq_CBZ*mTime, st_isc_cbz1[,1], lty = 2)
lines(tseq_CBZ*mTime, st_isc_cbz1[,2], lty = 2, lwd = 3)
lines(tseq_CBZ*mTime, st_isc_cbz1[,3], lty = 2)
legend(x = 1800, y = 1, legend = c("LTG", "CBZ"), lty = c(1,2), lwd = c(2,2))


#Mean curves


surv_UAE_LTG = exp(-(exp(-3.634-1.247+1.038*(Jointmodel$summary.random$`inla.group(timeL_LTG, n = 50)`[,2] + 2.434+0.310))*
                     Jointmodel$summary.random$`inla.group(timeL_LTG, n = 50)`[,1])^0.66)
surv_UAE_CBZ = exp(-(exp(-3.634+1.038*(Jointmodel$summary.random$`inla.group(timeL_CBZ, n = 50)`[,2] + 2.434))*
                       Jointmodel$summary.random$`inla.group(timeL_CBZ, n = 50)`[,1])^0.66)



plot(Jointmodel$summary.random$`inla.group(timeL_LTG, n = 50)`[,1],surv_UAE_LTG, type = "l")
lines(Jointmodel$summary.random$`inla.group(timeL_CBZ, n = 50)`[,1],surv_UAE_CBZ, type = "l", lty = 2)

surv_ISC_LTG = exp(-(exp(-2.947-0.305+1.055*(Jointmodel$summary.random$`inla.group(timeL_LTG, n = 50)`[,2] + 2.434+0.310))*
                       Jointmodel$summary.random$`inla.group(timeL_LTG, n = 50)`[,1])^1.01)
surv_ISC_CBZ = exp(-(exp(-2.947+1.055*(Jointmodel$summary.random$`inla.group(timeL_CBZ, n = 50)`[,2] + 2.434))*
                       Jointmodel$summary.random$`inla.group(timeL_CBZ, n = 50)`[,1])^1.01)

plot(Jointmodel$summary.random$`inla.group(timeL_LTG, n = 50)`[,1],surv_ISC_LTG, type = "l")
lines(Jointmodel$summary.random$`inla.group(timeL_CBZ, n = 50)`[,1],surv_ISC_CBZ, type = "l", lty = 2)


#########################################################################



#teeth data
library(bayesSurv)
data("tandmob2")
dataT = tandmob2
load("/Users/vanniej/Downloads/tooth24.RData")

dataT$d = rep(3, nrow(dataT))

for (i in 1:nrow(dataT)){
  if (is.na(dataT$EEND.24[i]))
    {dataT$d[i] = 0} else {
      if (is.na(dataT$EBEG.24[i]))
      {dataT$d[i] = 2} else
  if (dataT$EBEG.24[i] == dataT$EEND.24[i])
    dataT$d[i] = 1}
  
}
hist(dataT$d)
dataT$EBEG5.24 = dataT$EBEG.24-5
dataT$EEND5.24 = dataT$EEND.24-5

y = inla.surv(time = dataT$EBEG5.24, event = dataT$d, time2 = dataT$EEND5.24)
dataT$T24.DMF = rep(NA, nrow(dataT))
dataT$T24.DMF[tooth24$id] = tooth24$dmf

resT24 = inla(y ~ GENDER + T24.DMF,
              data = dataT,
              family = "weibullsurv",
              control.family = list(variant = 1),
              inla.mode = "experimental",
              control.compute = list(config = TRUE))
summary(resT24)


par(mar = c(5,5,1,1))
dataT1 = dataT[!is.na(dataT$EBEG5.24),]
npmleFit <- ic_np(cbind(EBEG5.24, EEND5.24) ~ GENDER, data = dataT1)
plot(npmleFit, main = "", col = c("red", "blue"), lty = c(1,2))


sample1 = inla.posterior.sample(resT24, n = 100)
dataT = cbind(dataT, matrix(rep(NA,nrow(dataT)*100), byrow = TRUE, nrow = nrow(dataT)))
hyper_dia = rep(NA,100)
betas = matrix(rep(NA, 300), ncol = 3)
for (i in 1:100){
  hyper_dia[i] = sample1[[i]]$hyperpar
  betas[i,1] = sample1[[i]]$latent[nrow(dataT)+1]
  betas[i,2] = sample1[[i]]$latent[nrow(dataT)+2]
  betas[i,3] = sample1[[i]]$latent[nrow(dataT)+3]
}


tseq = seq(0,8, length = 200)
st_male = matrix(rep(NA, 200*100), ncol = 100)
st_female = st_male
st_male_dmf1 = matrix(rep(NA, 200*100), ncol = 100)
st_female_dmf1 = st_male

for (i in 1:100){
  st_female[,i] = exp(-(exp(betas[i,1] + betas[i,2])*tseq)^hyper_dia[i])
  st_male[,i] = exp(-(exp(betas[i,1])*tseq)^hyper_dia[i])
  st_female_dmf1[,i] = exp(-(exp(betas[i,1] + betas[i,2] + betas[i,3])*tseq)^hyper_dia[i])
  st_male_dmf1[,i] = exp(-(exp(betas[i,1] + betas[i,3])*tseq)^hyper_dia[i])
}

st_m = matrix(rep(NA, 3*200), ncol = 3)
st_f = matrix(rep(NA, 3*200), ncol = 3)
st_m1 = matrix(rep(NA, 3*200), ncol = 3)
st_f1 = matrix(rep(NA, 3*200), ncol = 3)

for (i in 1:200){
  st_m[i,1] = quantile(st_male[i,], 0.025)
  st_m[i,2] = quantile(st_male[i,], 0.5)
  st_m[i,3] = quantile(st_male[i,], 0.975)
  st_f[i,1] = quantile(st_female[i,], 0.025)
  st_f[i,2] = quantile(st_female[i,], 0.5)
  st_f[i,3] = quantile(st_female[i,], 0.975)
  st_m1[i,1] = quantile(st_male_dmf1[i,], 0.025)
  st_m1[i,2] = quantile(st_male_dmf1[i,], 0.5)
  st_m1[i,3] = quantile(st_male_dmf1[i,], 0.975)
  st_f1[i,1] = quantile(st_female_dmf1[i,], 0.025)
  st_f1[i,2] = quantile(st_female_dmf1[i,], 0.5)
  st_f1[i,3] = quantile(st_female_dmf1[i,], 0.975)
}



yy = stepfun(npmleFit$scurves$girl$Tbull_ints[1:49],npmleFit$scurves$girl$S_curves$baseline)
zz = stepfun(npmleFit$scurves$boy$Tbull_ints[1:49],npmleFit$scurves$boy$S_curves$baseline)
plot(yy, do.points = FALSE, main = "", ylab = "S(t)", xlab = "Time (t)", col = "red", lwd = 2)
plot(zz, do.points = FALSE, add = TRUE, lty = 2, col = "blue", lwd = 2)

lines(tseq, st_m[,1], col = "blue", lty = 3)
lines(tseq, st_m[,2], col = "blue", lty = 2, lwd = 2)
lines(tseq, st_m[,3], col = "blue", lty = 3)
lines(tseq, st_f[,1], col = "red", lty = 3)
lines(tseq, st_f[,2], col = "red", lty = 1, lwd = 2)
lines(tseq, st_f[,3], col = "red", lty = 3)
lines(tseq, st_m1[,1], col = "blue", lty = 3)
lines(tseq, st_m1[,2], col = "blue", lty = 2, lwd = 2)
lines(tseq, st_m1[,3], col = "blue", lty = 3)
lines(tseq, st_f1[,1], col = "red", lty = 3)
lines(tseq, st_f1[,2], col = "red", lty = 1, lwd = 2)
lines(tseq, st_f1[,3], col = "red", lty = 3)
legend(x = 6.5, y = 1, legend = c("Female", "Male"), lty = c(1,2), col = c("red","blue"))

