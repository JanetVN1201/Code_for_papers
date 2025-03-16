library(runjags)
library(coda)

# Compute the q-value used in the truncated DW distribution
QValTDW <- function(mStar_vec, alpha, c) {
  exp(log(0.5)/(mStar_vec^(1/alpha) - c^(1/alpha)))
}

# Invert the truncated DW CDF
QTDW <- function(u, mStar_vec, alpha, c) {
  qVal <- QValTDW(mStar_vec, alpha, c)
  inside <- c^(1/alpha) + (log(1 - u)/log(qVal))
  val <- inside^alpha - 1
  y <- ceiling(val)
  y[y < c] <- c
  y
}

# Generate samples from truncated DW distribution
rTDW <- function(mStar_vec, alpha, c) {
  N <- length(mStar_vec)
  u <- runif(N)
  QTDW(u, mStar_vec, alpha, c)
}

# Generate samples from a contaminated truncated DW distribution
rcTDW <- function(mStar_vec, alpha, alpha2, delta, c) {
  N <- length(mStar_vec)
  picks <- runif(N) < delta
  out <- integer(N)
  out[picks] <- rTDW(mStar_vec[picks], alpha, c)
  out[!picks] <- rTDW(mStar_vec[!picks], alpha2, c)
  out
}

# Simulate regression data from cTDW distribution
SimRegcTDW <- function(N, c, beta0, beta1, alphaTrue, etaTrue, deltaTrue) {
  X <- rnorm(N, mean = 0, sd = 1)
  linpred <- beta0 + beta1*X
  mStar <- c + exp(linpred)
  alpha2 <- alphaTrue*etaTrue
  Y <- rcTDW(mStar, alphaTrue, alpha2, deltaTrue, c)
  data.frame(X = X, Y = Y)
}

# JAGS model for single-component TDW, using zero-trick
TDWJAGSModel <- "
model {
  beta0 ~ dnorm(0, 0.001)
  beta1 ~ dnorm(0, 0.001)
  alpha ~ dgamma(0.001, 0.001)
  Cpow <- pow(c, 1/alpha)
  for(i in 1:N) {
    log_mStarMinusC[i] <- beta0 + beta1*X[i]
    mStar[i] <- c + exp(log_mStarMinusC[i])
    A[i] <- pow(Y[i], 1/alpha)
    B[i] <- pow(Y[i] + 1.0, 1/alpha)
    logq[i] <- log(0.5)/(pow(mStar[i], 1/alpha) - Cpow)
    logNum[i] <- A[i]*logq[i] + log(1 - exp((B[i] - A[i])*logq[i]))
    logDen[i] <- Cpow*logq[i]
    logPMF[i] <- logNum[i] - logDen[i]
    zeros[i] ~ dpois(1E10 - logPMF[i])
  }
}
"

# JAGS model for contaminated TDW, using zero-trick
cTDWJAGSModel <- "
model {
  beta0 ~ dnorm(0, 0.001)
  beta1 ~ dnorm(0, 0.001)
  alpha ~ dgamma(0.001, 0.001)
  eta ~ dgamma(0.001, 0.001) T(1, )
  delta ~ dunif(0, 0.5)
  alpha2 <- alpha*eta
  C1 <- pow(c, 1/alpha)
  C2 <- pow(c, 1/alpha2)
  for(i in 1:N) {
    log_mStarMinusC[i] <- beta0 + beta1*X[i]
    mStar[i] <- c + exp(log_mStarMinusC[i])
    A1[i] <- pow(Y[i], 1/alpha)
    B1[i] <- pow(Y[i] + 1.0, 1/alpha)
    logq1[i] <- log(0.5)/(pow(mStar[i], 1/alpha) - C1)
    logNum1[i] <- A1[i]*logq1[i] + log(1 - exp((B1[i] - A1[i])*logq1[i]))
    logDen1[i] <- C1*logq1[i]
    logPMF_narrow[i] <- logNum1[i] - logDen1[i]
    A2[i] <- pow(Y[i], 1/alpha2)
    B2[i] <- pow(Y[i] + 1.0, 1/alpha2)
    logq2[i] <- log(0.5)/(pow(mStar[i], 1/alpha2) - C2)
    logNum2[i] <- A2[i]*logq2[i] + log(1 - exp((B2[i] - A2[i])*logq2[i]))
    logDen2[i] <- C2*logq2[i]
    logPMF_heavy[i] <- logNum2[i] - logDen2[i]
    logMixA[i] <- log(delta) + logPMF_narrow[i]
    logMixB[i] <- log(1 - delta) + logPMF_heavy[i]
    maxAB[i] <- ifelse(logMixA[i] > logMixB[i], logMixA[i], logMixB[i])
    diffAB[i] <- abs(logMixA[i] - logMixB[i])
    logPMF_mix[i] <- maxAB[i] + log(1 + exp(-diffAB[i]))
    zeros[i] ~ dpois(1E10 - logPMF_mix[i])
  }
}
"

# Fit single-component TDW model
FitTDW <- function(dat, cVal, burnin, sample, thin, adapt, nChains, beta0Start,
                   beta1Start, alphaStart) {
  N <- nrow(dat)
  jagsData <- list(
    N = N,
    X = dat$X,
    Y = dat$Y,
    c = cVal,
    zeros = rep(0, N)
  )
  initList <- replicate(
    nChains,
    list(
      beta0 = beta0Start,
      beta1 = beta1Start,
      alpha = alphaStart,
      .RNG.name = "base::Mersenne-Twister",
      .RNG.seed = sample.int(1e6, 1)
    ),
    simplify = FALSE
  )
  monitorPars <- c("beta0", "beta1", "alpha")
  fit <- autorun.jags(
    model = TDWJAGSModel,
    data = jagsData,
    inits = initList,
    monitor = monitorPars,
    n.chains = nChains,
    adapt = adapt,
    startburnin = burnin,
    startsample = sample,
    thin = thin,
    method = "parallel",
    psrf.target = 1.05,
    max.time = "60m",
    silent.jags = TRUE
  )
  fit
}

# Fit contaminated TDW model
FitcTDW <- function(dat, cVal, burnin, sample, thin, adapt, nChains, beta0Start,
                    beta1Start, alphaStart, etaStart, deltaStart) {
  N <- nrow(dat)
  jagsData <- list(
    N = N,
    X = dat$X,
    Y = dat$Y,
    c = cVal,
    zeros = rep(0, N)
  )
  initList <- replicate(
    nChains,
    list(
      beta0 = beta0Start,
      beta1 = beta1Start,
      alpha = alphaStart,
      eta = etaStart,
      delta = deltaStart,
      .RNG.name = "base::Mersenne-Twister",
      .RNG.seed = sample.int(1e6, 1)
    ),
    simplify = FALSE
  )
  monitorPars <- c("beta0", "beta1", "alpha", "eta", "delta")
  fit <- autorun.jags(
    model = cTDWJAGSModel,
    data = jagsData,
    inits = initList,
    monitor = monitorPars,
    n.chains = nChains,
    adapt = adapt,
    startburnin = burnin,
    startsample = sample,
    thin = thin,
    method = "parallel",
    psrf.target = 1.05,
    max.time = "60m",
    silent.jags = TRUE
  )
  fit
}

# Run multiple simulations of each model and store results
RunModels <- function(scenarioName, nsim, N, cVal, beta0True, beta1True,
                      alphaTrue, etaTrue, deltaTrue, burnin, sample, thin,
                      adapt, nChains) {
  resultsList <- vector("list", nsim)
  for(rep_i in seq_len(nsim)) {
    message(sprintf("=== Running replicate %d/%d ===", rep_i, nsim))
    simDat <- SimRegcTDW(N, cVal, beta0True, beta1True, alphaTrue, etaTrue,
                         deltaTrue)
    fitTDW <- FitTDW(
      dat = simDat,
      cVal = cVal,
      burnin = burnin,
      sample = sample,
      thin = thin,
      adapt = adapt,
      nChains = nChains,
      beta0Start = beta0True,
      beta1Start = beta1True,
      alphaStart = alphaTrue
    )
    fitcTDW <- FitcTDW(
      dat = simDat,
      cVal = cVal,
      burnin = burnin,
      sample = sample,
      thin = thin,
      adapt = adapt,
      nChains = nChains,
      beta0Start = beta0True,
      beta1Start = beta1True,
      alphaStart = alphaTrue,
      etaStart = etaTrue,
      deltaStart = deltaTrue
    )
    posteriorSummary <- function(fitObj, paramVec) {
      mm <- as.mcmc(fitObj)
      sumQ <- summary(mm)$quantiles
      lapply(paramVec, function(pa) {
        c(median = sumQ[pa, "50%"],
          lower = sumQ[pa, "2.5%"],
          upper = sumQ[pa, "97.5%"])
      })
    }
    tdwPars <- c("beta0", "beta1", "alpha")
    ctdwPars <- c("beta0", "beta1", "alpha", "eta", "delta")
    tdw_list <- posteriorSummary(fitTDW, tdwPars)
    ctdw_list <- posteriorSummary(fitcTDW, ctdwPars)
    list2df <- function(modelName, paramNames, smallList) {
      df0 <- do.call(rbind, smallList)
      data.frame(
        param = paramNames,
        estimate = df0[, "median"],
        cred_lower = df0[, "lower"],
        cred_upper = df0[, "upper"],
        model = modelName,
        replicate = rep_i,
        stringsAsFactors = FALSE
      )
    }
    tdw_df <- list2df("TDW", tdwPars, tdw_list)
    ctdw_df <- list2df("cTDW", ctdwPars, ctdw_list)
    all_df <- rbind(tdw_df, ctdw_df)
    resultsList[[rep_i]] <- all_df
  }
  bigDF <- do.call(rbind, resultsList)
  bigDF$scenario <- scenarioName
  bigDF$N <- N
  bigDF$cVal <- cVal
  bigDF$beta0True <- beta0True
  bigDF$beta1True <- beta1True
  bigDF$alphaTrue <- alphaTrue
  bigDF$etaTrue <- etaTrue
  bigDF$deltaTrue <- deltaTrue
  bigDF
}

# Save scenario results to a file
SaveModels <- function(scenarioName, nsim, N, cVal, beta0True, beta1True,
                       alphaTrue, etaTrue, deltaTrue, burnin, sample, thin,
                       adapt, nChains) {
  bigDF <- RunModels(
    scenarioName, nsim, N, cVal, beta0True, beta1True,
    alphaTrue, etaTrue, deltaTrue, burnin, sample, thin, adapt, nChains
  )
  scenarioObj <- list(
    date = Sys.Date(),
    time = format(Sys.time(), "%H:%M:%S"),
    scenarioName = scenarioName,
    simParams = list(
      nsim = nsim,
      N = N,
      cVal = cVal,
      beta0True = beta0True,
      beta1True = beta1True,
      alphaTrue = alphaTrue,
      etaTrue = etaTrue,
      deltaTrue = deltaTrue,
      burnin = burnin,
      sample = sample,
      thin = thin,
      adapt = adapt,
      nChains = nChains
    ),
    resultsLong = bigDF
  )
  outFile <- sprintf(
    "Scenario_%s_%s_%s.rds",
    scenarioName, Sys.Date(), format(Sys.time(), "%H%M%S")
  )
  saveRDS(scenarioObj, outFile)
  cat("Saved scenario results to:", outFile, "\n")
  scenarioObj
}