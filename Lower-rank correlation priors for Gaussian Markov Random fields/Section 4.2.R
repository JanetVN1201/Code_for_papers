
setwd(here::here("examples"))
getwd()

install.packages("INLAtools")
library(INLAtools)
library(graphpcor)

library(INLA)
library(ggplot2)
library(GGally)

glb <- INLA::inla.external.lib("graphpcor")




## STEP 1: data model definition

##   A data frame with 462 observations on the following 10 variables.
##  sbp systolic blood pressure
##  tobacco cumulative tobacco (kg)
##  ldl low density lipoprotein cholesterol
##  adiposity a numeric vector
##  famhist family history of heart disease,
##       a factor with levels ‘"Absent"’ and ‘"Present"’
##  typea type-A behavior
##  obesity a numeric vector
##  alcohol current alcohol consumption
##  age age at onset
##  chd response, coronary heart disease

## The chd variable is the outcome
## There are 8 continuous variables, but two are "non-Gaussian"
## and one binary (famhist) as covariates

## model idea:
##   Model chd (the response) is modeled as
##     chd[i] ~ bernoulli(inv_logit(alpha_1 + beta[j] * x[i,j]))
##       for i = 1, ..., n and j = 1, ..., p
##   where x is a latent variable modeled from the p covariates
##     x[i,] ~ N(0,R), R is a p = 9 dimensional correlation matrix
##   Model the first p-3 continuous variables as
##     y[i,j] = x[i,j], j = 1, ..., p-3
##   Model the next two (p-2 to p-1) contiuous variables as a Gamma(a, b=1)
##     y[i,j] ~ G1(a=exp(x[i,j]), j = p-2, p-1
##   Model famhist, coded as 0 (if 'Absent') or 1 (if 'Present') as
##     y[i,p] ~ bernoulli(inv_logit(alpha_2 + alpha_3 * x[i,p]))
## NOTE 1: y[,j], for j in 1, ..., p-2 were standardized
##         y[i,j] = x[i,j] for j in 1, ..., p-3
##         log(E(y[i,j])) = x[i,j] for j in p-2 and p-1
##         E(y[i,p]) = p_{ij} = 1/(1+exp(-alpha_2 - alpha_3 * x[i,p]))
## Note 2: I am not sure if we need alpha_3
##   (is there because famHist is the only covariate that is not scaled)
##   (fix it to one may be enough, its posterior is just above 1)
## Note 3: extension to allow error for the first p-1 covariates

data(SAheart, package = "msos")

str(SAheart)

(n <- nrow(SAheart))
(p <- ncol(SAheart)-1)

jj0 <- c(1, 3:4, 6:7, p)
jj1 <- c(2, 8)
jj2 <- c(jj0, jj1)

pm <- ggpairs(SAheart, columns = jj2,
              ggplot2::aes(colour = factor(chd), alpha = 0.5))

(pm)

xdata <- data.frame(
    scale(SAheart[, jj0]), ### first p-2 continuous (standardized)
    (SAheart[, jj1]<0.01)*0.01+SAheart[, jj1], ## p-2 and p-1
    famHist = (SAheart$famhist=="Present")+0L  ## p
); xnames <- colnames(xdata)
for(j in -2:-1+p) ## standardize p-2 and p-1 so its variance = 1
    xdata[, j] <- xdata[, j]/sd(xdata[, j])

sapply(xdata, sd)

cor_ci <- function(data, mapping, method = "pearson",
                   use = "complete.obs", ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)

  test <- cor.test(x, y, method = method)

  label <- sprintf("r = %.2f\nCI = [%.2f, %.2f]",
                   test$estimate,
                   test$conf.int[1],
                   test$conf.int[2])
  scol <- ifelse(test$conf.int[2]<0, "blue",
                 ifelse(test$conf.int[1]>0, "red", "gray"))
  ggplot(data = data, mapping = mapping) +
    annotate("text",
             x = mean(range(x, na.rm = TRUE)),
             y = mean(range(y, na.rm = TRUE)),
             label = label,
             color = scol,
             size = 4) +
    theme_void()
}

## Create customized pairs plot
ggpairs(
  SAheart[, c(jj2, 5)],
  upper = list(continuous = wrap("points", alpha = 0.6, size = 0.7)),
  diag  = list(continuous = wrap("barDiag", bins = 15, fill = "gray70")),
  lower = list(continuous = cor_ci)
)

round((cc <- cor(xdata))*100)

lcc <- chol(cc)
qc <- chol2inv(lcc)

## partial correlation matrix
pC <- cov2cor(qc)
dimnames(pC) <- dimnames(qc) <- dimnames(cc) <-
    list(xnames, xnames)
round(pC*100)

## define a graphpcor from a minimum spanning tree
nb <- lapply(1:p, function(i)
    setdiff(1:p,i)); class(nb) <- 'nb'
nbc <- lapply(1:p, function(i) 1-abs(pC[i, -i]))

library(spdep)
nbw <- nb2listw(nb, nbc, style="B")
mst <- mstree(nbw)

G0 <- matrix(0, p, p, dimnames = dimnames(pC))
for(i in 1:nrow(mst)) {
    G0[mst[i,1], mst[i,2]] <- 1
    G0[mst[i,2], mst[i,1]] <- 1
}
G0

library(graphpcor)
g0 <- graphpcor(G0)
(dg0 <- dim(g0))

p*(p-1)/2

tanh(0 + qnorm(c(0.025, 0.975)) / sqrt(n-3))

iil <- which(lower.tri(cc))

table(abs(pC[iil])>0.091)
table(abs(pC[iil])>0.09)
table(abs(pC[iil])>0.08)
table(abs(pC[iil])>0.07)
table(abs(pC[iil])>0.06)

g1 <- graphpcor(abs(pC)>0.07)
G1 <- attr(g1, "graph")
(dg1 <- dim(g1))

c(p=p, n=n)
g0
g1

lxy1 <- cbind(
    c(50,64,16,32,17,32,40,53,5),
    c(61,65,65,61,50,52,70,48,45)
)
lxy1 <- cbind(
    c(95, 33, 57, 78, 67, 61, 47, 55, 52),
    c(63, 55, 72, 63, 60, 53, 55, 47, 63)
)

lxy0 <- cbind(
    c(45, 15, 25, 35, 30, 30, 20, 20, 35),
    c(50, 50, 70, 60, 50, 40, 30, 20, 30)
)

par(mfrow = c(1, 2), mar = c(0,0,0,0))
plot(g0, Rgraphviz = TRUE)
##layout = lxy0, asp = 1,
  ##   edge.arrow.mode = 0, edge.color = "red", edge.width = 2,
    ## vertex.size = 25, ##vertex.shape = "sphere",
##     vertex.color = "lightblue", vertex.frame.color = "blue",
  ##   vertex.label.color = 'black', vertex.label.cex = 1.75,
    ## vertex.label.family = "Hershey", vertex.label.font = 1
##)
plot(g1, Rgraphviz=TRUE)
##     layout = lxy1,
  ##   edge.arrow.mode = 0, edge.color = "red",
    ## edge.width = 1+c(1,3,1,3,3,3,3,1,3,1,1,3,1,3),
##     vertex.size = 25, ##vertex.shape = "sphere",
  ##   vertex.color = "lightblue", vertex.frame.color = "blue",
    ## vertex.label.color = 'black', vertex.label.cex = 1.75,
     ##vertex.label.family = "Hershey", vertex.label.font = 1
##)

as.matrix(attr(g1, "graph"))

cg0 <- cgeneric(g0, lambda = 10,
                base = rep(0, dg0[2]),
                useINLAprecomp = TRUE)
#cg0$f$cgeneric$data$characters$shlib <- glb

cg1 <- cgeneric(g1, lambda = 10,
                base = rep(0, dg1[2]),
                useINLAprecomp = TRUE)

#library(INLA)

## create i and r copies
head(ic <- matrix(rep(1:p,each=n), n))
head(rc <- matrix(rep(1:n,p), n))
colnames(ic) <- paste("i", 1:p)
colnames(rc) <- paste("r", 1:p)

sdata <- inla.stack(
    inla.stack(
        tag = 'chd',
        data = list(chd = SAheart[, p+1]),
        effects = list(data.frame(
            alpha1 = 1, ic, rc)),
        A = list(1)),
    inla.stack(
        tag = "yx1",
        data = list(yx1 = unlist(xdata[, 1:(p-3)])),
        effects = list(data.frame(
            i = rep(1:(p-3), each = n),
            r = rep(1:n, p-3))),
        A = list(1)),
    inla.stack(
        tag = "yx2",
        data = list(yx2 = unlist(xdata[, -2:-1 + p])),
        effects = list(data.frame(
            i = rep((p-2):(p-1), each = n),
            r = rep(1:n, 2))),
        A = list(1)),
    inla.stack(
        tag = "famHist",
        data = list(famHist = xdata[, p]),
        effects = list(data.frame(
            alpha2 = rep(1, n),
            a3 = rep(p, n), ## as copy to fit alpha3
            a3r = 1:n, p)),
        A = list(1))
)

str(inla.stack.data(sdata))

ff0 <- paste(
    "list(chd, yx1, yx2, famHist)~0+alpha1+alpha2+",
    "f(i, model = cg0, replicate = r, n = 9, values = 1:9) +",
    paste(paste0("f(i.", 1:p, ", copy = 'i', replicate = r.",
                       1:p, ", fixed = FALSE)"), collapse = "+"),
    " + f(a3, copy = 'i', replicate = a3r, fixed = FALSE)")
ff0

ffg0 <- as.formula(ff0)
ffg1 <- as.formula(gsub("cg0", "cg1", ff0))
ffg1

cfam <- list(
    list(),
    list(hyper = list(prec = list(intial = 10, fixed = TRUE))),
    list(),
    list()
)

library(INLA)

if(file.exists("fit0.rds")) {
    fit0 <- readRDS("fit0.rds")
} else {
    fit0 <- inla(
        formula = ffg0,
        family = c("binomial", "gaussian", "gammajw", "binomial"),
        data = inla.stack.data(sdata),
        control.predictor = list(A = inla.stack.A(sdata)),
        control.family = cfam
    )
    saveRDS(fit0, "fit0.rds")
}

if(file.exists("fit1.rds")) {
    fit1 <- readRDS("fit1.rds")
} else {
    fit1 <- inla(
        formula = ffg1,
        family = c("binomial", "gaussian", "gammajw", "binomial"),
        data = inla.stack.data(sdata),
        control.predictor = list(A = inla.stack.A(sdata)),
        control.family = cfam
    )
    saveRDS(fit1, "fit1.rds")
}

rbind(fit0$cpu.used, fit1$cpu.used)

c(fit0$misc$nfunc, fit1$misc$nfunc)

fit0$summary.hy

summary(m0.chd <- glm(SAheart$chd~as.matrix(xdata), binomial))
sm0 <- coef(summary(m0.chd))
sm0 <- data.frame(sm0,
                  low = sm0[, 1] - 1.96 * sm0[, 2],
                  upp = sm0[, 1] + 1.96 * sm0[, 2])

sg0 <- fit0$summary.hyperpar[dg0[2] + 1:p,]
sg1 <- fit1$summary.hyperpar[dg1[2] + 1:p,]

## compare the regression coefficients with the simlple GLM: "M0"
par(mfrow = c(1, 1), mar = c(4,4,1,1), mgp = c(2,0.5,0), bty = 'n')
plot(1:p-0.2, coef(m0.chd)[-1], axes = FALSE, xlab = '',
     ylab = expression(beta[j]), xlim = c(0.5, p+0.5),
     ylim = range(sm0[-1, 5:6], sg0[, c(3,5)], sg1[, c(3,5)]), pch = 19)
segments(1:p-0.2, sm0[-1, 5], 1:p-0.2, sm0[-1, 6], lty = 1, lwd = 2)
abline(h = 0, lty = 2)
axis(2)
axis(1, 1:p, names(xdata))
points(1:p, sg0[, 1], pch = 19, col = 2)
segments(1:p, sg0[, 3], 1:p, sg0[, 5], col = 2, lwd = 2)
points(1:p+0.2, sg1[, 1], pch = 19, col = 4)
segments(1:p+0.2, sg1[, 3], 1:p+0.2, sg1[, 5], col = 4, lwd = 2)
abline(h = 0, lty = 2)
axis(2)
axis(1, 1:p, names(xdata))
legend("topleft", c("M0", "G0", "G1"), bty = "n",
       lty = 1, col = c(1,2,4), lwd = 2)

c0fit <- vcov(g0, theta = fit0$mode$theta[1:dg0[2]])
c1fit <- vcov(g1, theta = fit1$mode$theta[1:dg1[2]])

round(cc*100)
round(c0fit*100)
round(c1fit*100)


nmc <- 10000
h0sampls <- inla.hyperpar.sample(
    n = nmc, result = fit0, intern = TRUE
)
h1sampls <- inla.hyperpar.sample(
    n = nmc, result = fit1, intern = TRUE
)

iup <- which(upper.tri(cc))
c0sampls <- t(sapply(1:nmc, function(i) {
    vcov(g0, theta = h0sampls[i, 1:dg0[2]])[iup]
}))
c1sampls <- t(sapply(1:nmc, function(i) {
    vcov(g1, theta = h1sampls[i, 1:dg1[2]])[iup]
}))

lG0 <- upperPadding(G0)
lG1 <- upperPadding(G1)

fcols <- c(gray(0.35), rgb(0,0,1,0.7))

png("SAheartResultJoint.png", width = 4000, height = 3000, res = 300)
par(mfrow = c(p,p), mar = c(1.6, 1.6, 0.1, 0.1),
    mgp = c(1, 0.5, 0), bty = 'n')
kc <- k2 <- k1 <- k0 <- 0
for(i in 1:p) {
    for(j in 1:p) {
        if(i==j) {
            plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
            text(0, 0, xnames[j], cex = 2)
        }
        if(j>i) {
            iij0 <- which((i == (lG0@i+1) & (j == (lG0@j+1))))
            iij1 <- which((i == (lG1@i+1) & (j == (lG1@j+1))))
            if(length(iij0)>0) {
                k0 <- k0 + 1
                m0 <- inla.smarginal(fit0$internal.marginals.hyperpar[[k0]])
                h0 <- hist(h0sampls[, k0], 100, plot = FALSE)
                h0$ok <- TRUE
            } else {
                h0 <- list(ok=FALSE)
            }
            if(length(iij1)>0) {
                k1 <- k1 + 1
                m1 <- inla.smarginal(fit1$internal.marginals.hyperpar[[k1]])
                h1 <- hist(h1sampls[, k1], 100, plot = FALSE)
                h1$ok <- TRUE
            } else {
                h1 <- list(ok = FALSE)
            }
            if(h0$ok) {
                plot(m0, type = "l",
                     main = '', xlab = '', ylab = '',
                     xlim = range(0, h0$breaks, h1$breaks),
                     ylim = range(h0$dens, h1$dens),
                     col = fcols[1], border = 'transparent', lty = c(2), lwd = c(2))
                if(h1$ok) {
                    lines(m1, ##plot(h1, freq = FALSE, add = TRUE,
                         col = fcols[2], border = 'transparent', lwd = 2)
                    legend("topright", bty = "n",
                           as.expression(lapply(c(k0,k1), function(i) bquote(theta[.(i)]))),
                           col = fcols, border = 'transparent', lty = c(2,1), lwd = c(2,2))
                } else {
                    legend("topright", bty = "n",
                           as.expression(lapply(c(k0), function(i) bquote(theta[.(i)]))),
                           col = c(fcols[1]), border = 'transparent', lty = c(2,1), lwd = c(2,2))
                }
                rug(0, 0.1, lty = 3, lwd = 2)
            } else {
                if(h1$ok) {
                    plot(m1, type = "l", ##h1, freq = FALSE,
                         xlim = range(0, h1$breaks),
                         ylim = range(h0$dens, h1$dens),
                         main = '', xlab = '', ylab = '',
                         col = fcols[2], border = 'transparent')
                    legend("topleft", bty = "n",
                           as.expression(lapply(c(k1), function(i) bquote(theta[.(i)]))),
                           col = fcols[2], border = 'transparent')
                    rug(0, 0.1, lty = 3, lwd = 2)
                } else {
                    plot(0, type = 'n', axes = FALSE, xlab = '', ylab = '')
                    rug(0, 0.1, lty = 3, lwd = 2)
                }
            }
        }
        if(j<i) {
            k2 <- k2 + 1
            c_obs <- cc[iup[k2]]
            ic_obs <- tanh(c_obs + qnorm(c(0.025, 0.975)) / sqrt(n-3))
            h0 <- hist(c0sampls[, k2], 100, plot = FALSE)
            h1 <- hist(c1sampls[, k2], 100, plot = FALSE)
            plot(h0, freq = FALSE,
                 xlim = range(ic_obs, h0$breaks, h1$breaks),
                 ylim = range(h0$dens, h1$dens),
                 main = '', xlab = '', ylab = '',
                 col = fcols[1], border = 'transparent', density = 60, angle = 45)
            plot(h1, freq = FALSE, add = TRUE,
                 col = fcols[2], border = 'transparent')
            rug(c(c_obs, ic_obs), 0.1, lty = 1, lwd = 2, col = 6)
        }
      if(F){
        if((i==1) & (j==2)) {
            legend("bottom", "|", title = "Obs. IC", bty = "n",
                   lty = c(0), lwd = c(2), col = c(1), text.col = 6)
        }
        if((i==1) & (j==2)) {
            legend("top", title = "prior", bty = "n",
                   c("MST", "MST+"), lty = 2, lwd = 2, col = c(2, 4))
        }}
    }
}
dev.off()


if(FALSE) {
    system("eog SAheartResultJoint.png &")

}

