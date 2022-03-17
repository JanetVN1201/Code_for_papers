
library(INLA)

### P(y_it = 1) = 1/(1 + exp(-eta_{ij}))
### where
###  \eta_{ij} = \theta_j - alpha_i
###  \theta_j : ability (ideal point) of subject j
###  \alpha_i : dificulty of item i

### model settting
thfix <- list(theta=list(initial=log(1), fixed=TRUE))
ff <- y ~ 0 +
    f(j, model='iid', hyper=thfix, vb.correct=FALSE) +
    f(i, wi, model='iid', hyper=thfix, vb.correct=TRUE)

### number of items
ni <- 20

### control lp scale
betapriors <- vector('list', ni)
for(k in 1:ni)
    betapriors[[k]] <- list(prior='normal', param=c(1, 4), initial=1)
names(betapriors) <- paste0('theta', 1:ni)

### number of students
nj.s <- c(50, 100, 1000, 10000, 20000)
nj.s

###Times
times = matrix(NA, nrow = length(nj.s), ncol = 2)

s = 5
### model parameters
    set.seed(s)
    beta <- rgamma(ni, 20, 20)
    alpha <- rnorm(ni, 0, 1)
    theta <- rnorm(max(nj.s), 0, 1)

    for(k in 1:length(nj.s)) {

### set the abilities and data structure
        nj <- nj.s[k]
        longData <- list(
            i=rep(1:ni, nj),
            j=rep(1:nj, each=ni))
        longData$wi <- rep(-1, ni*nj)
        longData$b <- longData$i
        
### simulate the data     
        eta <- theta[longData$j] -alpha[longData$i]
        longData$y <- rbinom(
            n=ni*nj,
            size = 1,
            prob = 1/(1+exp(-eta*beta[longData$b])))

ind = TRUE
if (nj>20000){
    ind = FALSE
    } 

if (ind == TRUE){
res1 = inla(formula=ff,
                    family='binomial',
                    data=longData,
                    lp.scale=b,
                    control.lp.scale=list(hyper=betapriors),
            control.inla=list(cmin=1e-5))
}
res2 = inla(formula=ff,
            family='binomial',
            data=longData,
            lp.scale=b,
            control.lp.scale=list(hyper=betapriors),
            inla.mode = "experimental",
            control.inla=list(cmin=1e-5)) 
if (ind == TRUE){
times[k,1] = res1$cpu.used[4]
}
times[k,2] = res2$cpu.used[4]

    }

    



