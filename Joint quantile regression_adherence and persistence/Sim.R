library(INLA) # Bayesian inference with INLA
library(PermAlgo) # Permutation algorithm for survival times
library(mvtnorm)
library(tibble)
library(dplyr)
library(xtable)

# Working directory and read in datasets
#setwd('Quantile_joint')

# Solve alpha of qkumar
CalcAlpha <- function(Kappa, q, beta) {
  log(1 - (1 - q)^(1/beta))/log(Kappa)
}

# Solve beta of qkumar
CalcBeta <- function(psi, q) {
  log(1 - q)/log(1 - exp(-psi))
}

# Function to simulate from qkumar distribution
SimKUM <- function(n, Kappa, psi, q) {
  # Calculate alpha and beta
  beta <- CalcBeta(psi, q)
  alpha <- CalcAlpha(Kappa, q, beta)
  
  # Generate uniform samples
  u <- runif(n)
  
  # Transform using inverse CDF
  SimSample <- (1 - (1 - u)^(1/beta))^(1/alpha)
  
  return(SimSample)
}

# Function to perform simulation study
SimStudy <- function(ParamSet, iteration, scenario, cens0, q0) {
  print(paste0("Iteration ", iteration))
  
  params <- modifyList(SimParams, as.list(ParamSet))
  
  nsujet = params[[1]]
  psi = params[[2]]
  SD1 = params[[3]]
  SD2 = params[[4]]
  rho = params[[5]]
  gapLongi = params[[6]]
  gap = params[[7]]
 followup = params[[8]]
 q = q0
 beta_0 = params[[10]]
 beta_1 = params[[11]]
 beta_2 = params[[12]]
 beta_3 = params[[13]]
 gamma_1 = params[[14]]
 gamma_2 = params[[15]]
 phi = params[[16]]

  
  with(params, {
    Covar <- rho*SD1*SD2
    Sigma <- matrix(c(SD1^2, Covar, Covar, SD2^2), ncol = 2, nrow = 2)
    
    L <- t(chol(solve(Sigma)))
    diag(L) <- log(diag(L))

    theta_aa = log(phi)*0.1
    
    mestime = seq(0, followup, gap) # Measurement times
    timesLongi = mestime[which(round(mestime - round(mestime/gapLongi, 0)*gapLongi, 6) == 0)] # Visit times
    time = rep(mestime, nsujet) # Time column
    nmesindiv = followup/gap + 1 # Max number of individual measurements
    nmesy = nmesindiv*nsujet # Max total number of longitudinal measurements
    idY <- rep(1:nsujet, each = nmesindiv) # Individual ID
    
    MVnorm <- rmvnorm(nsujet, rep(0, 2), Sigma)
    b_int <- rep(MVnorm[, 1], each = nmesindiv) # Random intercept Y
    b_slo <- rep(MVnorm[, 2], each = nmesindiv) # Random slope Y
    
    binX = rbinom(nsujet, 1, 0.5) # Binary covariate for each subject
    ctsX = rnorm(nsujet, 1, 0.5) # Continuous covariate for each subject
    
    # Linear predictors
    #linPredY <- (beta_0 + b_int) + (beta_1 + b_slo)*time + beta_2*binX[idY] + beta_3*ctsX[idY]
    linPredY <- (beta_0 + b_int) + (beta_1 + b_slo)*time + beta_2*binX[idY] + beta_3*ctsX[idY]
    Y = SimKUM(nmesy, Kappa = plogis(linPredY), psi = psi, q = q)
    
    # Replicate binX and ctsX for each time point of each subject
    rep_binX <- rep(binX, each = nmesindiv)
    rep_ctsX <- rep(ctsX, each = nmesindiv)
    
    # Permutation algorithm to generate survival times dependent on linear predictors
    #DatTmp <- permalgorithm(nsujet, nmesindiv, Xmat = matrix(c(linPredY, rep_binX, rep_ctsX), nrow = nsujet*nmesindiv),
    DatTmp <- permalgorithm(nsujet, nmesindiv, Xmat = matrix(c(linPredY), nrow = nsujet*nmesindiv),
                            eventRandom = round(rexp(nsujet, cens0) + 1, 0), # ~40% events
                            censorRandom = runif(nsujet, 1, nmesindiv), # Uniform random censoring
                            XmatNames = c("linPredY"), # Predictor for association
                            betas = c(phi)) # Association parameter and survival regression parameters
    
    DatTmp2 = DatTmp[c(which(diff(DatTmp[,"Id"]) == 1), dim(DatTmp)[1]), c("Id", "Event", "Stop")]
    DatTmp2$eventTimes <- mestime[DatTmp2$Stop + 1] # Event times
    survDat <- DatTmp2[, c("Id", "eventTimes", "Event")]
    survDat$binX <- binX[survDat$Id] # Align binX with survival data
    survDat$ctsX <- ctsX[survDat$Id] # Align ctsX with survival data
    DatTmp$time <- mestime[DatTmp$Start + 1] # Measurement time of longitudinal outcomes
    DatTmp$Uid <- paste(DatTmp$Id, DatTmp$time) # Unique identifier to match covariates and longitudinal outcomes
    longDat3 <- merge(DatTmp[, c("Uid", "Id", "time")], 
                      cbind("Uid" = paste(idY, time), 
                            "binX" = rep_binX, 
                            "ctsX" = rep_ctsX, 
                            "Y" = Y), 
                      by = "Uid")
    longDat <- sapply(longDat3[longDat3$time %in% timesLongi, -1], as.numeric)
    longDat <- as.data.frame(longDat[order(longDat[, "Id"], longDat[, "time"]), ])
    
  #  print(paste0("Censoring percentage = ", round((1-sum(survDat$Event)/nsujet)*100,2)))
    # Plots to check feasibility of parameter combinations
    
    # Flag to control execution of plots
    runPlots <- FALSE
    
    if(runPlots) {
      
      library(ggplot2)
      
      # Profile plot of longitudinal outcomes
      Profplot <- ggplot(longDat, aes(x = time, y = Y, group = Id)) +
        geom_line() +
        theme_minimal()
      
      print(Profplot)
      
      # Histogram of survival outcomes
      print(hist(subset(survDat, Event == 0)$eventTimes))
      
    }
    
    # Model fit
    
    NL <- nrow(longDat)
    NS <- nrow(survDat)
    
    # Prepare data for INLA
    # Cox model structure for survival (with Bayesian smooth splines for baseline hazard, i.e., "RW2")
    
    #CoxExt <- inla.coxph(YS ~ -1 + Intercept + binX + ctsX, 
    CoxExt <- inla.coxph(YS ~ -1 + Intercept,  
                         control.hazard = list(model = "rw2", scale.model = TRUE, diagonal = 1e-2,
                                               constr = TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))),
                         data = list(YS = inla.surv(time = c(survDat$eventTimes), event = c(survDat$Event)),
                                     Intercept = rep(1, NS), 
                                     eta = survDat$Id, 
                                     binX = survDat$binX, 
                                     ctsX = survDat$ctsX))
    
    # Time weight for time-dependent covariates in survival (i.e., fixed and random slope)
    t.weight <- CoxExt$data$baseline.hazard.time + 0.5*CoxExt$data$baseline.hazard.length
    nsCox <- nrow(CoxExt$data) # Number of intervals for survival
    NY <- NL + nsCox
    IDcox <- CoxExt$data$expand..coxph # Random effects ID for survival
    CoxExt$data$eta <- 1:nsCox # Unique ID for shared linear predictor Y
    
    # Replicate binX and ctsX according to structure of CoxExt$data
    Rep_binX <- binX[CoxExt$data$expand..coxph]
    Rep_ctsX <- ctsX[CoxExt$data$expand..coxph]
    
    # Merge replicated covariates with CoxExt$data
    IDBC <- data.frame(CoxExt$data, binX = Rep_binX, ctsX = Rep_ctsX)
    
    covariates <- list(
      InteY = rep(1, NY), # Intercept Y
      TIMEY = c(longDat$time, t.weight), # Time Y
      binXY = c(longDat$binX, IDBC$binX), # Binary covariate X for Y
      ctsXY = c(longDat$ctsX, IDBC$ctsX), # Continuous covariate X for Y
      IDY = c(longDat$Id, IDcox), # Random intercept Y
      IDY_s = c(NS + longDat$Id, NS + IDcox), # Random slope Y
      WY_s = c(longDat$time, t.weight), # Weight random slope Y
      u1 = c(rep(NA, NL), CoxExt$data$eta), # Y association with survival
      w1 = c(rep(NA, NL), rep(-1, nsCox)) # Y weight association with survival
    )
    
    Y.joint <- list(
      Y = c(longDat$Y, rep(NA, nsCox)),
      Y.eta = c(rep(NA, NL), rep(0, nsCox)) # Y association with survival
    )
    
    jointdf <- data.frame(covariates, Y.joint)
    
    joint.DataCox <- c(as.list(inla.rbind.data.frames(jointdf, CoxExt$data)), CoxExt$data.list)
    Yjoint <- joint.DataCox[c("Y", "Y.eta", "y..coxph")] # Outcomes (longitudinal and survival)
    joint.DataCox$Y <- Yjoint
    
    # Update formula from Cox model structure to add longitudinal part
    #formulaJ <- update(CoxExt$formula, Yjoint ~ . -1 + InteY + TIMEY + binXY + ctsXY +
    formulaJ <- update(CoxExt$formula, Yjoint ~ . -1 + InteY + TIMEY + binXY + ctsXY +
                         f(IDY, model = "iidkd", order = 2, n = NS*2, constr = FALSE,
                           hyper = list(theta = list(param = c(7, 1, 1, 0)))) +
                         f(IDY_s, WY_s, copy ="IDY") +
                         f(u1, w1, model = "iid", hyper = list(prec = list(initial = -6, fixed = TRUE)), constr = FALSE) +
                         f(eta, copy = "u1", hyper = list(beta = list(fixed = FALSE, param = c(0, 1), initial = 0))))
    
    # INLA function call
    JMinla <- try(inla(verbose = FALSE,
                   formula = formulaJ,
                   family = c("qkumar", "Gaussian", CoxExt$family),
                   data = joint.DataCox,
                   control.family = list(
                     list(control.link = list(quantile = q)),
                     list(hyper = list(prec = list(initial = 12, fixed = TRUE))),
                     list()),
                   E = joint.DataCox$E..coxph,
                   control.inla = list(int.strategy = "eb")))
  
    if (!inherits(JMinla, "try-error")){
    # Extracting required summaries
    fixed <- JMinla$summary.fixed
    hyperpar <- JMinla$summary.hyperpar
    fixed <- cbind(Parameter = rownames(fixed), fixed)
    rownames(fixed) <- NULL
    hyperpar <- cbind(Parameter = rownames(hyperpar), hyperpar)
    rownames(hyperpar) <- NULL
    
    combined <- bind_rows(fixed, hyperpar)
    
    # xx <- inla.iidkd.sample(10^4, JMinla, "IDY")
    # 
    # VarCov <- matrix(unlist(xx), nrow = 2^2)
    # VarCovSD <- matrix(apply(VarCov,1,sd),2,2)
    # VarCov025 <- matrix(apply(VarCov,1,function(x) quantile(x, 0.025)),2,2)
    # VarCov975 <- matrix(apply(VarCov,1,function(x) quantile(x, 0.975)),2,2)
    # 
    # 
    # ## compute the mean
    # qq <- matrix(rowMeans(matrix(unlist(xx), nrow = 2^2)), 2, 2)
    # iSigma <- 1/sqrt(diag(Sigma))
    # Cor <- diag(iSigma) %*% Sigma %*% diag(iSigma)
    # round(dig = 3, cbind(inla = c(diag(qq), qq[lower.tri(qq)]),
    #                      true = c(sqrt(diag(Sigma)), Cor[lower.tri(Cor)])))
    
    
    # combined[combined$Parameter == "Theta1 for IDY","Parameter"] <- "SD of IDY (1)"
    # combined[combined$Parameter == "Theta2 for IDY","Parameter"] <- "SD of IDY (2)"
    # combined[combined$Parameter == "Theta3 for IDY","Parameter"] <- "rho of IDY"
    # 
    # combined[combined$Parameter == "SD of IDY (1)","mean"] <- sqrt(diag(qq)[1])
    # combined[combined$Parameter == "SD of IDY (2)","mean"] <- sqrt(diag(qq)[2])
    # combined[combined$Parameter == "rho of IDY","mean"] <- qq[lower.tri(qq)]
    # 
    # combined[combined$Parameter == "SD of IDY (1)","sd"] <-   VarCovSD[1,1]
    # combined[combined$Parameter == "SD of IDY (2)","sd"] <-   VarCovSD[2,2]
    # combined[combined$Parameter == "rho of IDY","sd"] <-   VarCovSD[1,2]
    # 
    # combined[combined$Parameter == "SD of IDY (1)","0.025quant"] <-   VarCov025[1,1]
    # combined[combined$Parameter == "SD of IDY (2)","0.025quant"] <-   VarCov025[2,2]
    # combined[combined$Parameter == "rho of IDY","0.025quant"] <-   VarCov025[1,2]
    # 
    # combined[combined$Parameter == "SD of IDY (1)","0.975quant"] <-   VarCov975[1,1]
    # combined[combined$Parameter == "SD of IDY (2)","0.975quant"] <-   VarCov975[2,2]
    # combined[combined$Parameter == "rho of IDY","0.975quant"] <-   VarCov975[1,2]
    # 
    
    ParamMap <- data.frame(
      Parameter = c("InteY", 
                    "TIMEY", 
                    "binXY", 
                    "ctsXY", 
                    "precision for qkumar observations",
                    "Theta1 for IDY", 
                    "Theta2 for IDY",
                    "Theta3 for IDY", 
                    "Beta for eta",
                    "binX",
                    "ctsX"),
      TableLabel = c("$\\beta_0$", 
                     "$\\beta_{\\text{time}}$", 
                     "$\\beta_{\\text{bin}}$", 
                     "$\\beta_{\\text{cont}}$", 
                     "$\\psi$", 
                     "$\\theta_{1}$", 
                     "$\\theta_{2}$", 
                     "$\\theta_{3}$", 
                     "$\\phi$",
                     "$\\gamma_{\\text{bin}}$",
                     "$\\gamma_{\\text{cont}}$"),
      Sort = 1:11
    )
    
    #combined <- bind_rows(fixed, hyperpar)
    combined <- merge(combined, ParamMap, by = "Parameter", all.x = TRUE)
    
    combined$TrueValue <- with(combined, case_when(
      Parameter == "InteY" ~ beta_0,
      Parameter == "TIMEY" ~ beta_1,
      Parameter == "binXY" ~ beta_2,
      Parameter == "ctsXY" ~ beta_3,
      Parameter == "precision for qkumar observations" ~ psi,
      Parameter == "Theta1 for IDY" ~ L[1,1],
      Parameter == "Theta2 for IDY" ~ L[2,2],
      Parameter == "Theta3 for IDY" ~ L[2,1], #Cor[lower.tri(Cor)]
      Parameter == "Beta for eta" ~ phi,
      Parameter == "binX" ~ gamma_1,
      Parameter == "ctsX" ~ gamma_2,
      TRUE ~ NA_real_
    ))
    
    # Scenario and iteration columns
    combined$Scenario <- paste("Scenario", scenario)
    combined$Iteration <- iteration
    
    return(list(combined = combined))}
    else {return(list(NA))}
  })
}

CombSim <- function(SimResult) {
  if (is.list(SimResult) && length(SimResult) > 0) {
    CombRes <- do.call(rbind, lapply(SimResult, function(x) x$combined))
    if (nrow(CombRes) > 0) {
      return(CombRes)
    } else {
      warning("Combined data frame is empty.")
      return(data.frame())
    }
  } else {
    warning("Input is not a list or is empty.")
    return(data.frame())
  }
}

# Function to identify non-static parameters
IdentifyNonStaticParams <- function(paramSets) {
  nonStaticParams <- names(paramSets)[sapply(paramSets, function(x) length(unique(x)) > 1)]
  return(nonStaticParams)
}

# Function to run simulations across all parameter sets
RunSimulations <- function(ParamSets, NSim, cens0, q0) {
  AllResult <- list()
  nonStaticParams <- IdentifyNonStaticParams(ParamSets)
  
  for (i in 1:nrow(ParamSets)) {
    set <- ParamSets[i, ]
    print(paste("Running simulation set:", i, "/", nrow(ParamSets)))
    
    SimResult <- lapply(seq_len(NSim), function(j) SimStudy(ParamSet = set, iteration = j, scenario = i, q0 = q0, cens0 = cens0))
    CombSimRes <- CombSim(SimResult)
    
    if (nrow(CombSimRes) > 0 && "TableLabel" %in% colnames(CombSimRes)) {
      CalcFreqRes <- do.call(rbind, lapply(split(CombSimRes, CombSimRes$TableLabel), function(df) {
        Bias <- mean(df$mean - df$TrueValue, na.rm = TRUE)
        RMSE <- sqrt(mean((df$mean - df$TrueValue)^2, na.rm = TRUE))
        Coverage <- mean(df$`0.025quant` <= df$TrueValue & df$`0.975quant` >= df$TrueValue, na.rm = TRUE)
        freqRes <- data.frame(
          Set = paste("Set", i),
          q = set$q,
          Sort = df$Sort[1],
          TableLabel = df$TableLabel[1],
          Value = df$TrueValue[1],
          Bias = Bias,
          RMSE = RMSE,
          Coverage = Coverage
        )
        for (col in nonStaticParams) {
          freqRes[[col]] <- set[[col]][1]
        }
        return(freqRes)
      }))
    } else {
      warning("CombSimRes is empty or does not contain expected columns.")
    }
    
    AllResult[[paste("Set", i)]] <- CalcFreqRes
  }
  
  return(AllResult)
}
