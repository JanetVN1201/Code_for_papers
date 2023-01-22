

# first run the data file 
nb <- poly2nb(map_2)
nb2INLA("map.adj", nb)
g <- inla.read.graph(filename = "map.adj")



# The Malaria quantile only

alpha = 0.2
r.m <- inla(formula = map_2$M.cases ~ 1 + offset(log(map_2$M.E)) +
              f(b1,model = "bym2", scale.model = T,
                graph = g)   ,
            family = "poisson",
            data =  merg.data,
            control.predictor = list(compute = T),
            control.compute = list(dic = T, waic = T, cpo = T),
            control.family = list(control.link = list(model = "quantile",
                                                      quantile = alpha)),
            verbose = TRUE)



round(r.m$summary.hyperpar[,c(1,3,5,6)],3)


# the quantile for G6PD only

alpha = 0.8
r.d <- inla(formula = map_2$D.cases ~ 1 + offset(log(map_2$D.E))+
              f(b2,model = "bym2", scale.model = T,
                graph = g)  ,
            family = "poisson",
            data =  merg.data,
            control.predictor = list(compute = T),
            control.compute = list(dic = T, waic = T, cpo = T),
            control.family = list(control.link = list(model = "quantile",
                                                      quantile = alpha)),
            verbose = TRUE)



round(r.d$summary.hyperpar[,c(1,3,5,6)],3)


# The joint quantile setup

b = length(merg.data$M.cases)
y1 = c(merg.data$M.cases , rep(NA,b))
y2 = c(rep(NA,b) , merg.data$D.cases)


E1 = c(merg.data$M.E,rep(NA,b))
E2 = c(rep(NA,b) , merg.data$D.E)

yy = cbind(y1,y2)
EE = cbind(E1, E2)
mu = c(rep(1, b) , rep(2,b))
b1 = c(1:b,rep(NA ,b))
b2= c(rep(NA,b),1:b)

besagproper = c(1:b , rep(NA ,b))

Besag.c = c(rep(NA , b), 1:b)

m = as.factor(mu)
d = data.frame(yy, m , besagproper , Besag.c,b1,b2)  


formula = yy ~ -1 +m + offset(log(E1))+ offset(log(E2))+
  f(besagproper, model = "besagproper", graph=g)+
  f(Besag.c, copy="besagproper", hyper = list(beta = list(fixed = FALSE)))+
  f(b1 , model = "bym2", graph=g, scale.model=TRUE)+
  f(b2, model = "bym2", graph=g, scale.model=TRUE)




alpha = 0.2
r1 <- inla(formula,
           family = c("poisson","poisson"),
           control.family = list(list(control.link = list(model = "quantile",quantile = alpha)),
                                 list(control.link = list(model = "quantile",
                                                          quantile = 1-alpha))), 
           data = d,
           verbose = F,
           control.compute = list(dic = T, waic = T, cpo = T),
           control.predictor = list(compute = T),
           inla.mode = "experimental")


#Other quantiles
r1a <- inla(formula,
            family = c("poisson","poisson"),
            control.family = list(list(control.link = list(model = "quantile",quantile = 0.8)),
                                  list(control.link = list(model = "quantile",
                                                           quantile = 0.3))),
            data = d,
            verbose = F,
            control.compute = list(dic = T, waic = T, cpo = T),
            control.predictor = list(compute = T),
            inla.mode = "experimental")

##Joint mean model
r2 <- inla(formula,
           family = c("poisson","poisson"),
           data = d,
           verbose = F,
           control.compute = list(dic = T, waic = T, cpo = T),
           control.predictor = list(compute = T),
           inla.mode = "experimental")







