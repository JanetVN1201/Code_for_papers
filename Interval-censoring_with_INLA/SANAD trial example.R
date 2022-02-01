library(joineR)
library(INLA)

#Scaling of time to (0;1) to prevent numerical overflow
mTime=max(max(epileptic$time),max(epileptic$with.time))
epileptic$time<-epileptic$time/mTime
epileptic$with.time <-epileptic$with.time/mTime

#Define the terms used for the splines for LTG and CBZ
epileptic$interaction <- with(epileptic, time * (treat == "LTG"))
epileptic$interaction2 <- with(epileptic, time * (treat == "CBZ"))
epileptic$interaction[epileptic$interaction==0]<-NA
epileptic$interaction2[epileptic$interaction2==0]<-NA

data1<-epileptic

#Create longitudinal dataset
dataL<-data.frame(ID=data1$id,LDose=log(data1$dose+0.1),
                  Time=data1$time,TimeSp=data1$time,Age=data1$age,
                  Gender=data1$gender,LD=data1$learn.dis,
                  Treatment=data1$treat,
                  InteractionLTG=data1$interaction,
                  InteractionCBZ=data1$interaction2,list(V=data1$id,W=data1$id),
                  Dose=data1$dose)

#Create survival dataset
dataS<-data.frame(ID=data1$id,Time=as.numeric(data1$with.time),
                  Status=data1$with.status2,StatusISC=data1$with.status.isc,
                  StatusUAE=data1$with.status.uae,Age=data1$age,
                  Gender=data1$gender,LD=data1$learn.dis,Treatment=data1$treat)
dataS<-subset(dataS,!duplicated(dataS$ID))

nL<-nrow(dataL)
nS<-nrow(dataS)

#Interval cesnoring where the last follow-up is less that the event time
dataS$int = rep(1, nrow(dataS))
dataS$time.l = dataS$Time
dataS$time.r = dataS$Time

for (i in dataS$ID){
  if (max(dataL$Time[dataL$ID==i])<dataS$Time[dataS$ID==i])
  {dataS$int[dataS$ID==i] = 3
  dataS$time.l[dataS$ID==i] = max(dataL$Time[dataL$ID==i])}
}

dataS$d.uae = dataS$StatusUAE*dataS$int
dataS$d.isc = dataS$StatusISC*dataS$int

#####INLA work
#Define the coavriates for the fixed and random effects for the 3 likelihoods
#and the copy of the longitudinal predictor

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

#Create a dataset that contains all the covariates and 3 response variables
jointdata<-c(fixed.eff,random.eff)
y.long <- c(dataL$Dose,rep(NA,nL),rep(NA, nS),rep(NA,nS))
#For the linear predictor to be shared
y.eta<-c(rep(NA,nL),rep(0,nL),rep(NA,nS),rep(NA,nS))
#Create the survival object for UAE
y.survUAE <- inla.surv(time = c(rep(NA, nL), rep(NA,nL),dataS$time.l,rep(NA,nS)), 
                       event = c(rep(NA, nL),rep(NA,nL),dataS$d.uae,rep(NA,nS)),
                       time2 = c(rep(NA, nL), rep(NA,nL),dataS$time.r,rep(NA,nS)))
#Create the survival object for UAE
y.survISC <- inla.surv(time = c(rep(NA, nL), rep(NA,nL),rep(NA,nS),dataS$time.l), 
                       event = c(rep(NA, nL),rep(NA,nL),rep(NA,nS),dataS$d.isc),
                       time2 = c(rep(NA, nL), rep(NA,nL),rep(NA,nS),dataS$time.r))
y.joint<-list(y.long,y.eta,y.survUAE,y.survISC)
jointdata$Y=y.joint

#Model formula
formula.model=Y~-1 + mu + treatmentL+treatmentUAE+treatmentISC+
f(inla.group(timeL_LTG,n=50),model="rw2",scale.model = TRUE,
   hyper = list(prec = list(prior="pc.prec", param=c(1, 0.01))))+
  f(inla.group(timeL_CBZ,n=50),model="rw2",scale.model = TRUE,
    hyper = list(prec = list(prior="pc.prec", param=c(1, 0.01))))+
  f(linpredL, linpredL2, model="iid", 
    hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(betaUAE, copy="linpredL", hyper = list(beta = list(fixed = FALSE)))+
  f(betaISC, copy="linpredL", hyper = list(beta = list(fixed = FALSE)))	

#Fit the model with the inla call
Jointmodel= inla(formula.model, family = c("gaussian","gaussian",
                                           "weibullsurv","weibullsurv"),
                 data = jointdata, verbose=TRUE, 
                 control.compute=list(config = TRUE),
                 control.family = list(
                   list(),
                   list(hyper = list(prec = list(initial = 10, fixed=TRUE))),
                   list(variant = 1),
                   list(variant = 1)
                 ), inla.mode = "experimental"
)
summary(Jointmodel)
