
library(INLA)

inla.setOption(
    pardiso.license='~/.pardiso.lic',
    smtp='pardiso')

### control.inla settings
ctri <- list(
    control.vb = list(enable = TRUE),
    strategy='gaussian',
    int.strategy='ccd')

### P(y_it = 1) = 1/(1 + exp(-eta_{ij}))
### where
###  \eta_{ij} = \theta_j - alpha_i
###  \theta_j : ability (ideal point) of subject j
###  \alpha_i : dificulty of item i

### load the data
load('ymat.RData')

n <- nrow(ymat) ### students
m <- ncol(ymat) ### items
print(c(students=n, items=m, ndata=n*m))

ldat <- list(j=rep(1:n, m),
             i=rep(1:m, each=n))
ldat$wi <- rep(-1, n*m)
ldat$y <- as.vector(ymat)
str(ldat)

### model settting
thfix <- list(theta=list(initial=log(1), fixed=TRUE))
ff <- y ~ 0 +
    f(j, model='iid', hyper=thfix, vb.correct=FALSE) +
    f(i, wi, model='iid', hyper=thfix, vb.correct=TRUE)

### control lp scale
betapriors <- vector('list', m)
for(k in 1:m)
    betapriors[[k]] <- list(prior='normal', param=c(1, 4), initial=1)
names(betapriors) <- paste0('theta', 1:m)

### fit the rash model
enem1pl <- inla(
    ff, 'binomial', data=ldat,
    control.predictor=list(link=1),
    verbose=TRUE,
    inla.mode = 'experimental', 
    control.inla = ctri,
    num.threads=paste0('1:-', detectCores()))

print(enem1pl$misc$nfunc)

save(list='enem1pl',
     file='enem1pl.RData')

### fit the 2pl model
ldat$b <- ldat$i
enem2pl <- inla(
    ff, 'binomial', data=ldat,
    control.predictor=list(link=1),
    verbose=TRUE,
    lp.scale=b,
    control.lp.scale=list(hyper=betapriors),
    inla.mode = 'experimental', 
    control.inla = ctri,
    num.threads=paste0(detectCores(), ':-1'))

print(enem2pl$misc$nfunc)

save(list='enem2pl',
     file='enem2pl.RData')


