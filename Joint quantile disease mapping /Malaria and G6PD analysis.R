

# first run the data file 
nb <- poly2nb(map_2_c)
nb2INLA("map.adj", nb)
g <- inla.read.graph(filename = "map.adj")



# The Malaria quantile only

alpha = 0.2
r.m <- inla(formula = M.cases ~ 1 + offset(log(M.E)) +
              f(b1,model = "bym2", scale.model = T,
                graph = g)   ,
            family = "poisson",
            data =  merg.data,
            control.predictor = list(compute = T),
            control.compute = list(dic = T, waic = T, cpo = T),
            control.family = list(control.link = list(model = "quantile",
                                                      quantile = alpha)),
            verbose = F)



summary(r.m)
round(r.m$summary.hyperpar[,c(1,3,5,6)],3)


# the quantile for G6PD only

alpha = 0.8
r.d <- inla(formula = D.cases ~ 1 + offset(log(D.E))+
              f(b2,model = "bym2", scale.model = T,
                graph = g)  ,
            family = "poisson",
            data =  merg.data,
            control.predictor = list(compute = T),
            control.compute = list(dic = T, waic = T, cpo = T),
            control.family = list(control.link = list(model = "quantile",
                                                      quantile = alpha)),
            verbose = F)


summary(r.d)
round(r.d$summary.hyperpar[,c(1,3,5,6)],3)

# The joint quantile setup
b = length(merg.data$M.cases)
y1 = c(merg.data$M.cases , rep(NA,b))
y2 = c(rep(NA,b) , merg.data$D.cases)

E1 = c(merg.data$M.E, rep(NA,b))
E2 = c(rep(NA,b), merg.data$D.E)

yy = cbind(y1,y2)
EE = cbind(E1, E2)
mu = c(rep(1, b) , rep(2,b))
b1 = c(1:b,rep(NA ,b))
b2 = c(rep(NA,b),1:b)

besagproper = c(1:b , rep(NA,b))

Besag.c = c(rep(NA,b), 1:b)

m = as.factor(mu)
d = data.frame(yy, m, besagproper, Besag.c, b1, b2, EE)  

formula = yy ~ -1 + m + offset(log(E1)) + offset(log(E2)) +
  f(b1, model = "bym2", graph=g, scale.model=TRUE) +
  f(b2, model = "bym2", graph=g, scale.model=TRUE) +
f(besagproper, model = "besagproper", graph=g, hyper = list(theta1 = list(param = c(1,0.1)))) +
  f(Besag.c, copy="besagproper", hyper = list(beta = list(fixed = FALSE, param = c(0,3))))


alpha = 0.2
r1 <- inla(formula,
           family = c("poisson","poisson"),
           control.family = list(list(control.link = list(model = "quantile", quantile = alpha)),
                                 list(control.link = list(model = "quantile",
                                                          quantile = 1-alpha))), 
           data = d,
           control.compute = list(dic = T, waic = T, cpo = T),
           control.predictor = list(compute = T),
           verbose = F)

summary(r1)

#Other quantiles
r1a <- inla(formula,
            family = c("poisson","poisson"),
            control.family = list(list(control.link = list(model = "quantile",quantile = 0.8)),
                                  list(control.link = list(model = "quantile",
                                                           quantile = 0.3))),
            data = d,
            verbose = F,
            control.compute = list(dic = T, waic = T, cpo = T),
            control.predictor = list(compute = T))

##Joint mean model
r2 <- inla(formula,
           family = c("poisson","poisson"),
           data = d,
           verbose = F,
           control.compute = list(dic = T, waic = T, cpo = T),
           control.predictor = list(compute = T),
           inla.mode = "experimental")
summary(r2)






