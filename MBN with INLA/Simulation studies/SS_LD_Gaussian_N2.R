library(bnlearn)
library(igraph)
library(INLA)
library(gridExtra)
library(ggplot2)
library(lme4)
library(dplyr)
library(mvtnorm)
library(caret)
library(graphpcor)

nsims = 100
## extending MBN using INLA to longtudenal data (correlated random effects $2 \times 2$ covariance matrix  )
## ======Custom score function MLE========##
fitME = function(node, parents, group, t, data) {
  
  # Remove 't' from parents to avoid duplication
  rhs = paste(c("1", setdiff(parents, t), t), collapse = " + ")  
  formula = as.formula(paste(node, "~", rhs, "+ (1 +", t, "|", group, ")"))
  
  model = try(lmer(formula, data = data, control = lmerControl(calc.derivs = FALSE)))
  
  if (is(model, "try-error"))
    return(NULL)
  
  return(model)
  
} #FITME

scoreME = function(node, parents, data, args = list()) {
  
  if (node == args$group)
    return(0)
  if (node == args$t)
    return(0)
  
  # mixed effects models with performance tuning to speed up learning.
  require(lme4)
  environment(nloptwrap)$defaultControl$maxeval = 5 * 1e3
  environment(nloptwrap)$defaultControl$xtol_abs = 1e-6
  environment(nloptwrap)$defaultControl$ftol_abs = 1e-6
  
  local.distribution =
    fitME(node = node, parents = setdiff(parents, c(args$group, args$t)),  # Ensure t is not duplicated
          group = args$group, t = args$t, data = data)
  
  if (is.null(local.distribution))
    return(-Inf)
  
  return(- BIC(local.distribution))
  
} #SCOREME

##=== Custom score function-- INLA default wishart prior ====##
fitME1 = function(node, parents, group, t, data) {
  nid <- length(unique(data[[group]]))  
  data$numid <- as.numeric(as.factor(data[[group]]))  
  data$slopeid <- data$numid +nid
  
  rhs = paste(c("1", setdiff(parents, t), t), collapse = " + ")
  
  formula = as.formula(paste0(node, "~", rhs, 
                              "+ f(numid, model = 'iidkd', n =", 2 * nid, ", order = 2)",
                              "+ f(slopeid, t, copy = 'numid')"))
  
  model = tryCatch(
    inla(formula, family = "gaussian", data = data,
         control.inla = list(int.strategy = "eb")),
    error = function(e) {
      message("Error fitting model for node ", node, ": ", e$message)
      return(NULL)
    }
  )
  return(model)
}

scoreME1 =  function(node, parents, data, args = list()) {
  
  if (node == args$group)
    return(0)
  if (node == args$t)
    return(0)
  
  require(INLA)
  
  local.distribution =
    fitME1(node = node, parents = setdiff(parents, c(args$group, args$t)),  # Ensure both are removed
           group = args$group, t = args$t, data = data)
  
  if (is.null(local.distribution))
    return(-Inf)
  
  return(local.distribution$mlik[1])
  
  
}

#== score function INLA_2 (user defind wishart prior)==#
fitME2 = function(node, parents, group, t, data) {
  nid <- length(unique(data[[group]]))  
  data$numid <- as.numeric(as.factor(data[[group]]))  
  data$slopeid <- data$numid +nid
  
  rhs = paste(c("1", setdiff(parents, t), t), collapse = " + ")
  
  formula = as.formula(paste0(
    node, "~", rhs, 
    "+ f(numid, model = 'iidkd', order = 2, n = ", 2 * nid, 
    ", constr = FALSE, hyper = list(theta1 = list(param = c(6,1,1,0))))",
    "+ f(slopeid, t, copy = 'numid')"
  ))
  
  model = tryCatch(
    inla(formula, family = "gaussian", data = data,
         control.inla = list(int.strategy = "eb")),
    error = function(e) {
      message("Error fitting model for node ", node, ": ", e$message)
      return(NULL)
    }
  )
  return(model)
}

scoreME2 =  function(node, parents, data, args = list()) {
  
  if (node == args$group)
    return(0)
  if (node == args$t)
    return(0)
  
  require(INLA)
  
  local.distribution =
    fitME2(node = node, parents = setdiff(parents, c(args$group, args$t)),  # Ensure both are removed
           group = args$group, t = args$t, data = data)
  
  if (is.null(local.distribution))
    return(-Inf)
  
  return(local.distribution$mlik[1])
}
#================== user defind prior 2
fitME3 = function(node, parents, group, t, data) {
  nid <- length(unique(data[[group]]))  
  data$numid <- as.numeric(as.factor(data[[group]]))  
  data$slopeid <- data$numid +nid
  
  rhs = paste(c("1", setdiff(parents, t), t), collapse = " + ")
  
  formula = as.formula(paste0(
    node, "~", rhs, 
    "+ f(numid, model = 'iidkd', order = 2, n = ", 2 * nid, 
    ", constr = FALSE, hyper = list(theta1 = list(param = c(4,1,1,0))))",
    "+ f(slopeid, t, copy = 'numid')"
  ))
  
  model = tryCatch(
    inla(formula, family = "gaussian", data = data,
         control.inla = list(int.strategy = "eb")),
    error = function(e) {
      message("Error fitting model for node ", node, ": ", e$message)
      return(NULL)
    }
  )
  return(model)
}

scoreME3 =  function(node, parents, data, args = list()) {
  
  if (node == args$group)
    return(0)
  if (node == args$t)
    return(0)
  
  require(INLA)
  
  local.distribution =
    fitME3(node = node, parents = setdiff(parents, c(args$group, args$t)),  # Ensure both are removed
           group = args$group, t = args$t, data = data)
  
  if (is.null(local.distribution))
    return(-Inf)
  
  return(local.distribution$mlik[1])
}

#== score function INLA- default LKJ prior)==#
fitME3 <- function(node, parents, group, t, data) {
  data$ID <- as.numeric(data$ID)
  data$rintercept <- rep(1, each = nrepeat * nindividual)
  data$rSlope <- rep(2, each = nrepeat * nindividual)
  n <- 2
  eta_values <- c(10, 30)  # Try these eta values in order
  
  cinla <- list(int.strategy = 'eb')
  rhs <- paste(c("1", setdiff(parents, t), t), collapse = " + ")
  
  for (eta in eta_values) {
    #message("Trying eta = ", eta, " for node ", node)
    cmodel <- cgeneric(model = "LKJ", n = n, eta = eta)
    
    formula <- as.formula(paste0(node, " ~ ", rhs,
                                 " + f(rintercept, model = cmodel, replicate = ID)",
                                 " + f(rSlope, ", t, ", copy = 'rintercept', replicate = ID)"))
    
    model <- try(
      inla(formula, family = "gaussian", data = data,
           control.mode = list(
             theta = c(4, 0),
             restart = TRUE),
           control.inla = cinla),
      silent = TRUE
    )
    
    # If the model is not an error object, return it
    if (!inherits(model, "try-error")) {
      return(model)
    } else {
      #message("INLA crashed with eta = ", eta, ". Trying next eta...")
    }
  }
  
  warning("All model attempts failed for node ", node)
  return(NULL)
}


scoreME3 =  function(node, parents, data, args = list()) {
  
  if (node == args$group)
    return(0)
  if (node == args$t)
    return(0)
  
  require(INLA)
  
  local.distribution =
    fitME3(node = node, parents = setdiff(parents, c(args$group, args$t)),  # Ensure both are removed
           group = args$group, t = args$t, data = data)
  
  if (is.null(local.distribution))
    return(-Inf)
  
  return(local.distribution$mlik[1])
  
  
}
#== score function INLA- user defind LKJ prior)==#
fitME4 <- function(node, parents, group, t, data) {
  data$ID <- as.numeric(data$ID)
  data$rintercept <- rep(1, each = nrepeat * nindividual)
  data$rSlope <- rep(2, each = nrepeat * nindividual)
  n <- 2
  eta_values <- c(13, 30)  # Try these eta values in order
  
  cinla <- list(int.strategy = 'eb')
  rhs <- paste(c("1", setdiff(parents, t), t), collapse = " + ")
  
  for (eta in eta_values) {
    #message("Trying eta = ", eta, " for node ", node)
    cmodel <- cgeneric(model = "LKJ", n = n, eta = eta)
    
    formula <- as.formula(paste0(node, " ~ ", rhs,
                                 " + f(rintercept, model = cmodel, replicate = ID)",
                                 " + f(rSlope, ", t, ", copy = 'rintercept', replicate = ID)"))
    
    model <- try(
      inla(formula, family = "gaussian", data = data,
           control.mode = list(
             theta = c(4, 0),
             restart = TRUE),
           control.inla = cinla),
      silent = TRUE
    )
    
    # If the model is not an error object, return it
    if (!inherits(model, "try-error")) {
      return(model)
    } else {
      # message("INLA crashed with eta = ", eta, ". Trying next eta...")
    }
  }
  
  warning("All model attempts failed for node ", node)
  return(NULL)
}


scoreME4 =  function(node, parents, data, args = list()) {
  
  if (node == args$group)
    return(0)
  if (node == args$t)
    return(0)
  
  require(INLA)
  
  local.distribution =
    fitME4(node = node, parents = setdiff(parents, c(args$group, args$t)),  # Ensure both are removed
           group = args$group, t = args$t, data = data)
  
  if (is.null(local.distribution))
    return(-Inf)
  
  return(local.distribution$mlik[1])
  
  
}

F1_score_structure1 <- function(TBN, LBN) {
  TBN = remove.node(TBN, "G")
  LBN = remove.node(LBN, "G")
  TBN = cpdag(TBN)
  LBN = cpdag(LBN)
  true_arcs <- arcs(TBN)
  learned_arcs <- arcs(LBN)
  
  arc_to_string <- function(x) apply(x, 1, paste, collapse = "~")
  
  true_set <- arc_to_string(true_arcs)
  learned_set <- arc_to_string(learned_arcs)
  
  TP <- sum(learned_set %in% true_set) 
  FP <- sum(!(learned_set %in% true_set))
  FN <- sum(!(true_set %in% learned_set))
  
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  
  F1 <- 2 * precision * recall / (precision + recall)
  
  return(F1)
}
#=====Parameter learning function lme======#
## Parameter learning

# MLE
ldistME_lme = function(node, parents, group, t, data) {
  require(lme4)
  
  environment(nloptwrap)$defaultControl$maxeval = 5 * 1e3
  environment(nloptwrap)$defaultControl$xtol_abs = 1e-6
  environment(nloptwrap)$defaultControl$ftol_abs = 1e-6
  
  model = fitME(node = node, parents = parents, group = group, t=t, data = data)
  ngroups = nlevels(data[[group]])
  
  ldist = list(
    coef = array(0, dim = c(length(parents) + 2, ngroups), 
                 dimnames = list(c("(Intercept)", "t", parents), NULL)),
    sd = rep(0, ngroups)
  )
  
  if (is.null(model) || isSingular(model, tol = 1e-4)) {
    environment(nloptwrap)$defaultControl$maxeval = 1e5
    environment(nloptwrap)$defaultControl$xtol_abs = 1e-5
    environment(nloptwrap)$defaultControl$ftol_abs = 1e-5
    
    model = fitME(node = node, parents = parents, group = group, t=t, data = data)
    
    environment(nloptwrap)$defaultControl$maxeval = 5 * 1e3
    environment(nloptwrap)$defaultControl$xtol_abs = 1e-6
    environment(nloptwrap)$defaultControl$ftol_abs = 1e-6
    
    if (is.null(model)) {
      warning("Model is NULL for node ", node)
      return(ldist)
    }
    
    if (isSingular(model, tol = 1e-4)) {
      warning("Singular fit for node ", node, " but extracting coefficients anyway.")
      print(summary(model))  
    }
  }
  
  fixefs = fixef(model)
  coef_names = names(fixefs)
  
  for (i in seq_len(ncol(ldist$coef))) {
    ldist$coef[, i] = 0 
    for (name in coef_names) {
      if (name %in% rownames(ldist$coef)) {
        ldist$coef[name, i] = fixefs[name]
      }
    }
  }
  
  ranefs = ranef(model)[[group]]
  group_levels = levels(data[[group]])
  
  for (i in seq_len(ngroups)) {
    ldist$coef["(Intercept)", i] = ldist$coef["(Intercept)", i] + as.numeric(ranefs[group_levels[i], "(Intercept)"])
    ldist$coef["t", i] = ldist$coef["t", i] + as.numeric(ranefs[group_levels[i], "t"])
  }
  resids = resid(model)
  for (i in seq_len(ngroups)) {
    ldist$sd[i] = sd(resids[data[[group]] == group_levels[i]])
  }
  
  # Handle NAs
  ldist$coef[is.na(ldist$coef)] = 0
  ldist$sd[is.na(ldist$sd)] = 0
  
  return(ldist)
}

#=== parameter learning function INLA (default prior) ===#

ldistME_inla = function(node, parents, group, t, data) {
  require(INLA)
  model = fitME1(node = node, parents = parents, group = group, t=t, data = data)
  ngroups = nlevels(data[[group]])
  
  ldist = list(
    coef = array(0, dim = c(length(parents) + 2, ngroups), 
                 dimnames = list(c("(Intercept)", "t", parents), NULL)),
    sd = rep(0, ngroups)
  )
  
  if (is.null(model)) {
    
    model = fitME1(node = node, parents = parents, group = group, t=t, data = data)
    if (is.null(model)) {
      warning("Model is NULL for node ", node)
      return(ldist)
    }
    
    if (any(is.na(model$summary.fixed))) {
      warning("Singular fit for node ", node, " but extracting coefficients anyway.")
      print(model$summary.fixed)  
    }
  }  
  
  fixefs = model$summary.fixed[, "mean"]
  coef_names = rownames(model$summary.fixed)
  names(fixefs) = coef_names
  coef_names = rownames(model$summary.fixed)
  
  
  for (i in seq_len(ncol(ldist$coef))) {
    ldist$coef[, i] = 0  # Initialize with zeros
    for (name in coef_names) {
      if (name %in% rownames(ldist$coef)) {
        ldist$coef[name, i] = fixefs[name]
      }
    }
  }
  
  ranefID <- model$summary.random$numid$mode[1:30]
  raneft  <- model$summary.random$numid$mode[31:60]
  
  ranefs <- data.frame(
    "(Intercept)" = ranefID,
    t = raneft,
    check.names = FALSE
  )
  group_levels = levels(data[[group]])
  
  for (i in seq_len(ngroups)) {
    ldist$coef["(Intercept)", i] = ldist$coef["(Intercept)", i] + as.numeric(ranefs[group_levels[i], "(Intercept)"])
    ldist$coef["t", i] = ldist$coef["t", i] + as.numeric(ranefs[group_levels[i], "t"])
    
  }
  
  pred = model$summary.fitted.values$mode
  resids = data[[node]] - pred
  
  for (i in seq_len(ngroups)) {
    ldist$sd[i] = sd(resids[data[[group]] == group_levels[i]])
  }
  
  # Handle NAs
  ldist$coef[is.na(ldist$coef)] = 0
  ldist$sd[is.na(ldist$sd)] = 0
  
  return(ldist)
}

#=== parameter learning function INLA (user defind prior) ===#

ldistME_inla1 = function(node, parents, group, t, data) {
  require(INLA)
  model = fitME2(node = node, parents = parents, group = group, t=t, data = data)
  ngroups = nlevels(data[[group]])
  
  ldist = list(
    coef = array(0, dim = c(length(parents) + 2, ngroups), 
                 dimnames = list(c("(Intercept)", "t", parents), NULL)),
    sd = rep(0, ngroups)
  )
  
  if (is.null(model)) {
    
    model = fitME2(node = node, parents = parents, group = group, t=t, data = data)
    if (is.null(model)) {
      warning("Model is NULL for node ", node)
      return(ldist)
    }
    
    if (any(is.na(model$summary.fixed))) {
      warning("Singular fit for node ", node, " but extracting coefficients anyway.")
      print(model$summary.fixed)  
    }
  }  
  
  fixefs = model$summary.fixed[, "mean"]
  coef_names = rownames(model$summary.fixed)
  names(fixefs) = coef_names
  coef_names = rownames(model$summary.fixed)
  
  
  for (i in seq_len(ncol(ldist$coef))) {
    ldist$coef[, i] = 0  # Initialize with zeros
    for (name in coef_names) {
      if (name %in% rownames(ldist$coef)) {
        ldist$coef[name, i] = fixefs[name]
      }
    }
  }
  
  ranefID <- model$summary.random$numid$mode[1:30]
  raneft  <- model$summary.random$numid$mode[31:60]
  
  ranefs <- data.frame(
    "(Intercept)" = ranefID,
    t = raneft,
    check.names = FALSE
  )
  group_levels = levels(data[[group]])
  
  for (i in seq_len(ngroups)) {
    ldist$coef["(Intercept)", i] = ldist$coef["(Intercept)", i] + as.numeric(ranefs[group_levels[i], "(Intercept)"])
    ldist$coef["t", i] = ldist$coef["t", i] + as.numeric(ranefs[group_levels[i], "t"])
    
  }
  
  pred = model$summary.fitted.values$mode
  resids = data[[node]] - pred
  
  for (i in seq_len(ngroups)) {
    ldist$sd[i] = sd(resids[data[[group]] == group_levels[i]])
  }
  
  # Handle NAs
  ldist$coef[is.na(ldist$coef)] = 0
  ldist$sd[is.na(ldist$sd)] = 0
  
  return(ldist)
}

#=== parameter learning function INLA (default lkj prior, eta = 10) ===#
ldistME_inla2 = function(node, parents, group, t, data) {
  require(INLA)
  model = fitME3(node = node, parents = parents, group = group, t=t, data = data)
  ngroups = nlevels(data[[group]])
  
  ldist = list(
    coef = array(0, dim = c(length(parents) + 2, ngroups), 
                 dimnames = list(c("(Intercept)", "t", parents), NULL)),
    sd = rep(0, ngroups)
  )
  
  if (is.null(model)) {
    
    model = fitME3(node = node, parents = parents, group = group, t=t, data = data)
    if (is.null(model)) {
      warning("Model is NULL for node ", node)
      return(ldist)
    }
    
    if (any(is.na(model$summary.fixed))) {
      warning("Singular fit for node ", node, " but extracting coefficients anyway.")
      print(model$summary.fixed)  
    }
  }  
  
  fixefs = model$summary.fixed[, "mean"]
  coef_names = rownames(model$summary.fixed)
  names(fixefs) = coef_names
  coef_names = rownames(model$summary.fixed)
  
  
  for (i in seq_len(ncol(ldist$coef))) {
    ldist$coef[, i] = 0  # Initialize with zeros
    for (name in coef_names) {
      if (name %in% rownames(ldist$coef)) {
        ldist$coef[name, i] = fixefs[name]
      }
    }
  }
  ranef = matrix(model$summary.random$rintercept$mean,
                 byrow = TRUE, nrow = nindividual)
  ranefID = ranef[, 1]
  raneft  = ranef[, 2]
  
  ranefs <- data.frame(
    "(Intercept)" = ranefID,
    t = raneft,
    check.names = FALSE
  )
  group_levels = levels(data[[group]])
  
  for (i in seq_len(ngroups)) {
    ldist$coef["(Intercept)", i] = ldist$coef["(Intercept)", i] + as.numeric(ranefs[group_levels[i], "(Intercept)"])
    ldist$coef["t", i] = ldist$coef["t", i] + as.numeric(ranefs[group_levels[i], "t"])
    
  }
  
  pred = model$summary.fitted.values$mode
  resids = data[[node]] - pred
  
  for (i in seq_len(ngroups)) {
    ldist$sd[i] = sd(resids[data[[group]] == group_levels[i]])
  }
  
  # Handle NAs
  ldist$coef[is.na(ldist$coef)] = 0
  ldist$sd[is.na(ldist$sd)] = 0
  
  return(ldist)
}

#=== parameter learning function INLA (user defind lkj prior, eta = 13) ===#
ldistME_inla3 = function(node, parents, group, t, data) {
  require(INLA)
  model = fitME4(node = node, parents = parents, group = group, t=t, data = data)
  ngroups = nlevels(data[[group]])
  
  ldist = list(
    coef = array(0, dim = c(length(parents) + 2, ngroups), 
                 dimnames = list(c("(Intercept)", "t", parents), NULL)),
    sd = rep(0, ngroups)
  )
  
  if (is.null(model)) {
    
    model = fitME4(node = node, parents = parents, group = group, t=t, data = data)
    if (is.null(model)) {
      warning("Model is NULL for node ", node)
      return(ldist)
    }
    
    if (any(is.na(model$summary.fixed))) {
      warning("Singular fit for node ", node, " but extracting coefficients anyway.")
      print(model$summary.fixed)  
    }
  }  
  
  fixefs = model$summary.fixed[, "mean"]
  coef_names = rownames(model$summary.fixed)
  names(fixefs) = coef_names
  coef_names = rownames(model$summary.fixed)
  
  
  for (i in seq_len(ncol(ldist$coef))) {
    ldist$coef[, i] = 0  # Initialize with zeros
    for (name in coef_names) {
      if (name %in% rownames(ldist$coef)) {
        ldist$coef[name, i] = fixefs[name]
      }
    }
  }
  ranef = matrix(model$summary.random$rintercept$mean,
                 byrow = TRUE, nrow = nindividual)
  ranefID = ranef[, 1]
  raneft  = ranef[, 2]
  
  ranefs <- data.frame(
    "(Intercept)" = ranefID,
    t = raneft,
    check.names = FALSE
  )
  group_levels = levels(data[[group]])
  
  for (i in seq_len(ngroups)) {
    ldist$coef["(Intercept)", i] = ldist$coef["(Intercept)", i] + as.numeric(ranefs[group_levels[i], "(Intercept)"])
    ldist$coef["t", i] = ldist$coef["t", i] + as.numeric(ranefs[group_levels[i], "t"])
    
  }
  
  pred = model$summary.fitted.values$mode
  resids = data[[node]] - pred
  
  for (i in seq_len(ngroups)) {
    ldist$sd[i] = sd(resids[data[[group]] == group_levels[i]])
  }
  
  # Handle NAs
  ldist$coef[is.na(ldist$coef)] = 0
  ldist$sd[is.na(ldist$sd)] = 0
  
  return(ldist)
}
#== KLD to compaire parametrs of mle and INLA with the true parameters ===#
KL = function(P, Q) {
  Pdag = remove.node(bn.net(P), "ID")
  Qdag = remove.node(bn.net(Q), "ID")
  klvalue = 0  
  
  for (group in seq_along(coef(P$ID))) {  
    
    Pldist = structure(vector(nnodes(Pdag), mode = "list"), names = nodes(Pdag))
    
    for (node in nodes(Pdag)) {
      if (is(P[[node]], "bn.fit.cgnode")) {
        Pldist[[node]] = list(coef = coef(P[[node]])[, group], 
                              sd = sigma(P[[node]])[group])
      } else {
        Pldist[[node]] = list(coef = coef(P[[node]]), sd = sigma(P[[node]]))
      }
    }
    
    Pgroup = custom.fit(Pdag, Pldist)
    
    Qldist = structure(vector(nnodes(Qdag), mode = "list"), names = nodes(Qdag))
    
    for (node in nodes(Qdag)) {
      if (is(Q[[node]], "bn.fit.cgnode")) {
        Qldist[[node]] = list(coef = coef(Q[[node]])[, group], 
                              sd = sigma(Q[[node]])[group])
      } else {
        Qldist[[node]] = list(coef = coef(Q[[node]]), sd = sigma(Q[[node]]))
      }
    }
    
    Qgroup = custom.fit(Qdag, Qldist)
    
    Pmvn = gbn2mvnorm(Pgroup)
    Qmvn = gbn2mvnorm(Qgroup)
    order <- names(Qmvn$mu)
    Pmvn$mu <- Pmvn$mu[order]
    Pmvn$sigma <- Pmvn$sigma[order, order]
    
    if (isTRUE(all.equal(Pmvn, Qmvn))) {
      next
    }
    
    svdQ = svd(Qmvn$sigma)
    det.sigmaQ = prod(svdQ$d)
    pos = svdQ$d > svdQ$d[1] * .Machine$double.eps
    inv.sigmaQ = svdQ$v[, pos, drop = FALSE] %*%
      diag(1 / svdQ$d[pos], nrow = length(svdQ$d[pos])) %*%
      t(svdQ$u[, pos, drop = FALSE])
    
    svdP = svd(Pmvn$sigma)
    det.sigmaP = prod(svdP$d)
    
    if ((det.sigmaP == 0) || (det.sigmaQ == 0)) {
      return(Inf)
    }
    
    klvalue = klvalue + 
      0.5 * as.numeric(
        sum(log(svdQ$d)) - sum(log(svdP$d)) + 
          sum(diag(inv.sigmaQ %*% Pmvn$sigma)) - nrow(Qmvn$sigma) + 
          t(Pmvn$mu - Qmvn$mu) %*% inv.sigmaQ %*% (Pmvn$mu - Qmvn$mu)
      )
  }
  
  return(klvalue)
}

F1_score = function(fitted) {
  
  pred = predict(fitted, node = "ID", data = data, method = "bayes-lw")
  obs = data$ID
  
  if (nlevels(pred) == 2)
    F1s = confusionMatrix(pred, obs)[["byClass"]]["F1"]
  else
    F1s = confusionMatrix(pred, obs)[["byClass"]][, "F1"]
  
  return(mean(F1s, na.rm = TRUE))
  
}

F1_score_structure <- function(TBN, LBN) {
  true_arcs <- arcs(TBN)
  learned_arcs <- arcs(LBN)
  
  arc_to_string <- function(x) apply(x, 1, paste, collapse = "~")
  
  true_set <- arc_to_string(true_arcs)
  learned_set <- arc_to_string(learned_arcs)
  
  TP <- sum(learned_set %in% true_set) 
  FP <- sum(!(learned_set %in% true_set))
  FN <- sum(!(true_set %in% learned_set))
  
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  
  F1 <- 2 * precision * recall / (precision + recall)
  
  return(F1)
}

F1_score_structure1 <- function(TBN, LBN) {
  TBN = remove.node(TBN, "ID")
  TBN = remove.node(TBN, "t")
  LBN = remove.node(LBN, "ID")
  LBN = remove.node(LBN, "t")
  TBN = cpdag(TBN)
  LBN = cpdag(LBN)
  true_arcs <- arcs(TBN)
  learned_arcs <- arcs(LBN)
  
  arc_to_string <- function(x) apply(x, 1, paste, collapse = "~")
  
  true_set <- arc_to_string(true_arcs)
  learned_set <- arc_to_string(learned_arcs)
  
  TP <- sum(learned_set %in% true_set) 
  FP <- sum(!(learned_set %in% true_set))
  FN <- sum(!(true_set %in% learned_set))
  
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  
  F1 <- 2 * precision * recall / (precision + recall)
  
  return(F1)
}
##===  Simulation study for longitudinal data with wishart prior ====##

# Random DAG used to simulate the data
TMBN = model2network("[ID][t][Y1|ID:t][Y3|ID:t][Y2|Y1:ID:t][Y5|Y3:ID:t][Y6|Y1:Y3:ID:t][Y4|Y2:ID:t][Y8|Y5:ID:t][Y7|Y4:ID:t][Y9|Y3:Y4:Y5:ID:t][Y10|Y4:ID:t]")

shd_lme <- list()
shd_inla <-list()
shd_inla1 <-list()
shd_inla2 <-list()
shd_inla3 <-list()

KL_lme <- list()
KL_inla <- list()
KL_inla1 <- list()
KL_inla2 <- list()
KL_inla3 <- list()

F1_lme = list()
F1_inla = list()
F1_inla1 = list()
F1_inla2 = list()
F1_inla3 = list()

F1_Slme <- list()
F1_Sinla <- list()
F1_Sinla1 <- list()
F1_Sinla2 <- list()
F1_Sinla3 <- list()

F1_Slme_1 <- list()
F1_Sinla_1 <- list()
F1_Sinla1_1 <- list()
F1_Sinla2_1 <- list()
F1_Sinla3_1 <- list()

for(k in 1:nsims) {
  nindividual = 30
  nrepeat = 5
  N = nrepeat*nindividual
  
  data = data.frame(
    "ID" = rep(1:nindividual, each = nrepeat),
    "t" = rep(1:nrepeat, nindividual)
  )
  #data$t <- data$t - 1
  #=== Y1|ID:t
  m= 2
  rho <- 0.5
  Sigma <- matrix(NA, m, m)
  diag(Sigma) <- c(1, 0.7)
  for(i in 1:m) {
    for (j in 1:m) {
      if (i != j) {
        Sigma[i, j] <- rho^abs(i-j) * sqrt(Sigma[i, i] * Sigma[j, j])
      }
    }
  }
  random_effects <- rmvnorm(n = nindividual,sigma = Sigma)
  
  data$u0 = rep(random_effects[, 1],
                each = nrepeat, times=1)
  
  data$u1 = rep(random_effects[, 2],
                each = nrepeat, times=1)
  
  sigmae1 = sqrt(0.5)
  data$e = rnorm(N, 0, sigmae1)
  
  data$b0 = 1
  data$b1 = 0.5
  
  data$Y1 = data$b0 + data$b1*data$t + data$u1*data$t + data$u0 + data$e
  
  ##
  data$b0i = data$b0 + data$u0
  data$b1i = data$b1 + data$u1
  
  b0 = list()
  for(i in 1:nindividual){
    b0[i] =  data$b0i[data$ID == i][1]
  }
  b0 = unlist(b0)
  
  b1 = list()
  for(i in 1:nindividual){
    b1[i] =  data$b1i[data$ID == i][1]
  }
  b1 =unlist(b1)
  
  sde = list()
  for(i in 1:nindividual){
    sde[i] =  sd(data$e[data$ID == i])
  }
  sde = unlist(sde)
  
  
  CPTY1 = list(
    coef = matrix(c(b0, b1), nrow = 2, ncol = nindividual, byrow = TRUE,
                  dimnames = list(c("(Intercept)", "t"), as.character(seq(0, nindividual-1)))),
    sd = sde
  )
  data <- data[, !(names(data) %in% c("u0","u1", "e", "b0","b1","b0i","b1i"))]
  
  
  #== Y2|Y1:ID:t
  m= 2
  rho <- 0.43
  Sigma <- matrix(NA, m, m)
  diag(Sigma) <- c(0.68, 0.52)
  for(i in 1:m) {
    for (j in 1:m) {
      if (i != j) {
        Sigma[i, j] <- rho^abs(i-j) * sqrt(Sigma[i, i] * Sigma[j, j])
      }
    }
  }
  random_effects <- rmvnorm(n = nindividual,sigma = Sigma)
  
  data$u0 = rep(random_effects[, 1],
                each = nrepeat, times=1)
  
  data$u1 = rep(random_effects[, 2],
                each = nrepeat, times=1)
  
  sigmae1 = sqrt(0.55)
  data$e = rnorm(N, 0, sigmae1)
  
  data$b0 = -1.03
  data$b1 = 0.5
  data$b2 = 0.35
  
  data$Y2 = data$b0 + data$b1*data$t + data$b2*data$Y1 + data$u1*data$t + data$u0 + data$e
  
  ##
  data$b0i = data$b0 + data$u0
  data$b1i = data$b1 + data$u1
  
  b0 = list()
  for(i in 1:nindividual){
    b0[i] =  data$b0i[data$ID == i][1]
  }
  b0 = unlist(b0)
  
  b1 = list()
  for(i in 1:nindividual){
    b1[i] =  data$b1i[data$ID == i][1]
  }
  b1 =unlist(b1)
  
  b2 = list()
  for(i in 1:nindividual){
    b2[i] =  data$b2[data$ID == i][1]
  }
  b2 =unlist(b2)
  
  sde = list()
  for(i in 1:nindividual){
    sde[i] =  sd(data$e[data$ID == i])
  }
  sde = unlist(sde)
  
  
  CPTY2 = list(
    coef = matrix(c(b0, b1, b2), nrow = 3, ncol = nindividual, byrow = TRUE,
                  dimnames = list(c("(Intercept)", "t", "Y1"), as.character(seq(0, nindividual-1)))),
    sd = sde
  )
  data <- data[, !(names(data) %in% c("u0","u1", "e", "b0","b1","b2","b0i","b1i"))]
  
  #=== Y3|ID:t
  m= 2
  rho <- -0.5
  Sigma <- matrix(NA, m, m)
  diag(Sigma) <- c(0.45, 0.35)
  for(i in 1:m) {
    for (j in 1:m) {
      if (i != j) {
        Sigma[i, j] <- rho^abs(i-j) * sqrt(Sigma[i, i] * Sigma[j, j])
      }
    }
  }
  random_effects <- rmvnorm(n = nindividual,sigma = Sigma)
  
  data$u0 = rep(random_effects[, 1],
                each = nrepeat, times=1)
  
  data$u1 = rep(random_effects[, 2],
                each = nrepeat, times=1)
  
  sigmae1 = sqrt(0.53) # 0.6
  data$e = rnorm(N, 0, sigmae1)
  
  data$b0 = 1.3
  data$b1 = 0.3
  #data$b2 = 0.4
  
  data$Y3 = data$b0 + data$b1*data$t +
            data$u1*data$t + data$u0 + data$e
  
  ##
  
  data$b0i = data$b0 + data$u0
  data$b1i = data$b1 + data$u1
  
  b0 = list()
  for(i in 1:nindividual){
    b0[i] =  data$b0i[data$ID == i][1]
  }
  b0 = unlist(b0)
  
  b1 = list()
  for(i in 1:nindividual){
    b1[i] =  data$b1i[data$ID == i][1]
  }
  b1 =unlist(b1)
  
  sde = list()
  for(i in 1:nindividual){
    sde[i] =  sd(data$e[data$ID == i])
  }
  sde = unlist(sde)
  
  CPTY3 = list(
    coef = matrix(c(b0, b1), nrow = 2, ncol = nindividual, byrow = TRUE,
                  dimnames = list(c("(Intercept)", "t"), as.character(seq(0, nindividual-1)))),
    sd = sde
  )
  data <- data[, !(names(data) %in% c("u0","u1", "e", "b0","b1","b0i","b1i", "b2"))]
  
  #=== Y4|Y2:ID:t
  m= 2
  rho <- 0.4
  Sigma <- matrix(NA, m, m)
  diag(Sigma) <- c(0.61, 0.48)
  for(i in 1:m) {
    for (j in 1:m) {
      if (i != j) {
        Sigma[i, j] <- rho^abs(i-j) * sqrt(Sigma[i, i] * Sigma[j, j])
      }
    }
  }
  random_effects <- rmvnorm(n = nindividual,sigma = Sigma)
  
  data$u0 = rep(random_effects[, 1],
                each = nrepeat, times=1)
  
  data$u1 = rep(random_effects[, 2],
                each = nrepeat, times=1)
  
  sigmae1 = sqrt(0.57) # 0.9
  data$e = rnorm(N, 0, sigmae1)
  
  data$b0 = -0.32
  data$b1 = 0.56
  data$b2 = 0.48 #1.02
  
  data$Y4 = data$b0 + data$b1*data$t + data$b2*data$Y2 + 
    data$u1*data$t + data$u0 + data$e
  
  ##
  data$b0i = data$b0 + data$u0
  data$b1i = data$b1 + data$u1
  
  b0 = list()
  for(i in 1:nindividual){
    b0[i] =  data$b0i[data$ID == i][1]
  }
  b0 = unlist(b0)
  
  b1 = list()
  for(i in 1:nindividual){
    b1[i] =  data$b1i[data$ID == i][1]
  }
  b1 =unlist(b1)
  
  b2 = list()
  for(i in 1:nindividual){
    b2[i] =  data$b2[data$ID == i][1]
  }
  b2 =unlist(b2)
  
  sde = list()
  for(i in 1:nindividual){
    sde[i] =  sd(data$e[data$ID == i])
  }
  sde = unlist(sde)
  
  CPTY4 = list(
    coef = matrix(c(b0, b1, b2), nrow = 3, ncol = nindividual, byrow = TRUE,
                  dimnames = list(c("(Intercept)", "t", "Y2"), as.character(seq(0, nindividual-1)))),
    sd = sde
  )
  data <- data[, !(names(data) %in% c("u0","u1", "e", "b0","b1","b0i","b1i", "b2"))]
  
  #===Y5|Y3:ID:t
  m= 2
  rho <- 0.6 # 0.4
  Sigma <- matrix(NA, m, m)
  diag(Sigma) <- c(0.62, 0.50)
  for(i in 1:m) {
    for (j in 1:m) {
      if (i != j) {
        Sigma[i, j] <- rho^abs(i-j) * sqrt(Sigma[i, i] * Sigma[j, j])
      }
    }
  }
  random_effects <- rmvnorm(n = nindividual, sigma = Sigma)
  
  data$u0 = rep(random_effects[, 1],
                each = nrepeat, times=1)
  
  data$u1 = rep(random_effects[, 2],
                each = nrepeat, times=1)
  
  sigmae1 = sqrt(0.5) #0.95
  data$e = rnorm(N, 0, sigmae1)
  
  data$b0 = 2.63
  data$b1 = 1.19
  data$b2 = -0.5
  
  data$Y5 = data$b0 + data$b1*data$t + data$b2*data$Y3 + 
    data$u1*data$t + data$u0 + data$e
  
  ##
  data$b0i = data$b0 + data$u0
  data$b1i = data$b1 + data$u1
  
  b0 = list()
  for(i in 1:nindividual){
    b0[i] =  data$b0i[data$ID == i][1]
  }
  b0 = unlist(b0)
  
  b1 = list()
  for(i in 1:nindividual){
    b1[i] =  data$b1i[data$ID == i][1]
  }
  b1 =unlist(b1)
  
  b2 = list()
  for(i in 1:nindividual){
    b2[i] =  data$b2[data$ID == i][1]
  }
  b2 =unlist(b2)
  
  sde = list()
  for(i in 1:nindividual){
    sde[i] =  sd(data$e[data$ID == i])
  }
  sde = unlist(sde)
  
  
  CPTY5 = list(
    coef = matrix(c(b0, b1, b2), nrow = 3, ncol = nindividual, byrow = TRUE,
                  dimnames = list(c("(Intercept)", "t", "Y3"), as.character(seq(0, nindividual-1)))),
    sd = sde
  )
  data <- data[, !(names(data) %in% c("u0","u1", "e", "b0","b1","b0i","b1i", "b2"))]
  
  
  #===Y6|Y1,Y3:t:ID
  m= 2
  rho <- 0.48
  Sigma <- matrix(NA, m, m)
  diag(Sigma) <- c(0.7, 0.53)
  for(i in 1:m) {
    for (j in 1:m) {
      if (i != j) {
        Sigma[i, j] <- rho^abs(i-j) * sqrt(Sigma[i, i] * Sigma[j, j])
      }
    }
  }
  random_effects <- rmvnorm(n = nindividual,sigma = Sigma)
  
  data$u0 = rep(random_effects[, 1],
                each = nrepeat, times=1)
  
  data$u1 = rep(random_effects[, 2],
                each = nrepeat, times=1)
  
  
  
  sigmae1 = sqrt(0.4) #0.4
  data$e = rnorm(N, 0, sigmae1)
  
  data$b0 = 1.03
  data$b1 = -0.27
  data$b2 = 0.77
  data$b3 = 1.25
  
  data$Y6 = data$b0 + data$b1*data$t + data$b2*data$Y1 + 
    data$b3*data$Y3 + data$u1*data$t + 
    data$u0 + data$e
  
  ##
  data$b0i = data$b0 + data$u0
  data$b1i = data$b1 + data$u1
  
  b0 = list()
  for(i in 1:nindividual){
    b0[i] =  data$b0i[data$ID == i][1]
  }
  b0 = unlist(b0)
  
  b1 = list()
  for(i in 1:nindividual){
    b1[i] =  data$b1i[data$ID == i][1]
  }
  b1 =unlist(b1)
  
  b2 = list()
  for(i in 1:nindividual){
    b2[i] =  data$b2[data$ID == i][1]
  }
  b2 =unlist(b2)
  
  b3 = list()
  for(i in 1:nindividual){
    b3[i] =  data$b3[data$ID == i][1]
  }
  b3 =unlist(b3)
  
  sde = list()
  for(i in 1:nindividual){
    sde[i] =  sd(data$e[data$ID == i])
  }
  sde = unlist(sde)
  
  
  CPTY6 = list(
    coef = matrix(c(b0, b1, b2, b3), nrow = 4, ncol = nindividual, byrow = TRUE,
                  dimnames = list(c("(Intercept)", "t", "Y1","Y3"), as.character(seq(0, nindividual-1)))),
    
    sd = sde
  )
  data <- data[, !(names(data) %in% c("u0","u1", "e", "b0","b1","b0i","b1i", "b2","b3"))]
  
  
  #===Y7|Y4:ID:t
  m= 2
  rho <- 0.54
  Sigma <- matrix(NA, m, m)
  diag(Sigma) <- c(0.75, 0.58)
  for(i in 1:m) {
    for (j in 1:m) {
      if (i != j) {
        Sigma[i, j] <- rho^abs(i-j) * sqrt(Sigma[i, i] * Sigma[j, j])
      }
    }
  }
  random_effects <- rmvnorm(n = nindividual, sigma = Sigma)
  
  data$u0 = rep(random_effects[, 1],
                each = nrepeat, times=1)
  
  data$u1 = rep(random_effects[, 2],
                each = nrepeat, times=1)
  
  sigmae1 = sqrt(0.57) # 0.87
  data$e = rnorm(N, 0, sigmae1)
  
  data$b0 = 1.28
  data$b1 = -0.43
  data$b2 = 0.53
 # data$b3 =  -0.8
  
  data$Y7 = data$b0 + data$b1*data$t + data$b2*data$Y4 + 
            data$u1*data$t + data$u0 + data$e
  
  ##
  data$b0i = data$b0 + data$u0
  data$b1i = data$b1 + data$u1
  
  b0 = list()
  for(i in 1:nindividual){
    b0[i] =  data$b0i[data$ID == i][1]
  }
  b0 = unlist(b0)
  
  b1 = list()
  for(i in 1:nindividual){
    b1[i] =  data$b1i[data$ID == i][1]
  }
  b1 =unlist(b1)
  
  b2 = list()
  for(i in 1:nindividual){
    b2[i] =  data$b2[data$ID == i][1]
  }
  b2 =unlist(b2)
  
  sde = list()
  for(i in 1:nindividual){
    sde[i] =  sd(data$e[data$ID == i])
  }
  sde = unlist(sde)
  
  
  CPTY7 = list(
    coef = matrix(c(b0, b1, b2), nrow = 3, ncol = nindividual, byrow = TRUE,
                  dimnames = list(c("(Intercept)", "t", "Y4"), as.character(seq(0, nindividual-1)))),
    sd = sde
  )
  data <- data[, !(names(data) %in% c("u0","u1", "e", "b0","b1","b0i","b1i", "b2","b3"))]
  
  
  #== Y8|Y5:ID:t
  m= 2
  rho <- 0.38
  Sigma <- matrix(NA, m, m)
  diag(Sigma) <- c(0.65, 0.48)
  for(i in 1:m) {
    for (j in 1:m) {
      if (i != j) {
        Sigma[i, j] <- rho^abs(i-j) * sqrt(Sigma[i, i] * Sigma[j, j])
      }
    }
  }
  
  random_effects <- rmvnorm(n = nindividual, sigma = Sigma)
  
  data$u0 = rep(random_effects[, 1],
                each = nrepeat, times=1)
  
  data$u1 = rep(random_effects[, 2],
                each = nrepeat, times=1)
  
  sigmae1 = sqrt(0.54) # 0.52
  data$e = rnorm(N, 0, sigmae1)
  
  data$b0 = 1
  data$b1 = 0.5
  data$b2 = -0.6
  data$Y8 = data$b0 + data$b1*data$t + data$b2*data$Y5  + data$u1*data$t + data$u0 + data$e
  
  ##
  data$b0i = data$b0 + data$u0
  data$b1i = data$b1 + data$u1
  
  b0 = list()
  for(i in 1:nindividual){
    b0[i] =  data$b0i[data$ID == i][1]
  }
  b0 = unlist(b0)
  
  b1 = list()
  for(i in 1:nindividual){
    b1[i] =  data$b1i[data$ID == i][1]
  }
  b1 =unlist(b1)
  
  b2 = list()
  for(i in 1:nindividual){
    b2[i] =  data$b2[data$ID == i][1]
  }
  b2 =unlist(b2)
  
  sde = list()
  for(i in 1:nindividual){
    sde[i] =  sd(data$e[data$ID == i])
  }
  sde = unlist(sde)
  
  CPTY8 = list(
    coef = matrix(c(b0, b1, b2), nrow = 3, ncol = nindividual, byrow = TRUE,
                  dimnames = list(c("(Intercept)", "t", "Y5"), as.character(seq(0, nindividual-1)))),
    sd = sde
  )
  data <- data[, !(names(data) %in% c("u0","u1", "e", "b0","b1","b2","b0i","b1i"))]
  
  #===Y9|Y1,Y4,Y8:ID:time
  m= 2
  rho <- 0.5
  Sigma <- matrix(NA, m, m)
  diag(Sigma) <- c(0.67, 0.56)
  for(i in 1:m) {
    for (j in 1:m) {
      if (i != j) {
        Sigma[i, j] <- rho^abs(i-j) * sqrt(Sigma[i, i] * Sigma[j, j])
      }
    }
  }
  random_effects <- rmvnorm(n = nindividual, sigma = Sigma)
  
  data$u0 = rep(random_effects[, 1],
                each = nrepeat, times=1)
  
  data$u1 = rep(random_effects[, 2],
                each = nrepeat, times=1)
  
  sigmae1 = sqrt(0.56) # 0.85
  data$e = rnorm(N, 0, sigmae1)
  
  data$b0 = 0.97
  data$b1 = 0.23
  data$b2 = 1.31
  data$b3 = -0.34
  data$b4 = -0.68
  
  data$Y9 = data$b0 +data$b1*data$t + data$b2*data$Y3 +data$b3*data$Y4 +data$b4*data$Y5 + 
    data$u1*data$t + data$u0 + data$e
  
  ##
  data$b0i = data$b0 + data$u0
  data$b1i = data$b1 + data$u1
  
  b0 = list()
  for(i in 1:nindividual){
    b0[i] =  data$b0i[data$ID == i][1]
  }
  b0 = unlist(b0)
  
  b1 = list()
  for(i in 1:nindividual){
    b1[i] =  data$b1i[data$ID == i][1]
  }
  b1 =unlist(b1)
  
  b2 = list()
  for(i in 1:nindividual){
    b2[i] =  data$b2[data$ID == i][1]
  }
  b2 =unlist(b2)
  
  b3 = list()
  for(i in 1:nindividual){
    b3[i] =  data$b3[data$ID == i][1]
  }
  b3 =unlist(b3)
  
  b4 = list()
  for(i in 1:nindividual){
    b4[i] =  data$b4[data$ID == i][1]
  }
  b4 =unlist(b4)
  
  sde = list()
  for(i in 1:nindividual){
    sde[i] =  sd(data$e[data$ID == i])
  }
  sde = unlist(sde)
  
  
  CPTY9 = list(
    coef = matrix(c(b0, b1, b2, b3, b4), nrow = 5, ncol = nindividual, byrow = TRUE,
                  dimnames = list(c("(Intercept)", "t", "Y3","Y4", "Y5"), as.character(seq(0, nindividual-1)))),
    sd = sde
  )
  data <- data[, !(names(data) %in% c("u0","u1", "e", "b0","b1","b0i","b1i", "b2","b3","b4"))]
  
  
  #===Y10|Y4:ID:time
  m= 2
  rho <- 0.53
  Sigma <- matrix(NA, m, m)
  diag(Sigma) <- c(0.56, 0.48)
  for(i in 1:m) {
    for (j in 1:m) {
      if (i != j) {
        Sigma[i, j] <- rho^abs(i-j) * sqrt(Sigma[i, i] * Sigma[j, j])
      }
    }
  }
  random_effects <- rmvnorm(n = nindividual, sigma = Sigma)
  
  data$u0 = rep(random_effects[, 1],
                each = nrepeat, times=1)
  
  data$u1 = rep(random_effects[, 2],
                each = nrepeat, times=1)
  
  
  sigmae1 = sqrt(0.55) # 0.75
  data$e = rnorm(N, 0, sigmae1)
  
  data$b0 = 0.34
  data$b1 = 0.87
  data$b2 = 1.12
  
  
  data$Y10 = data$b0 + data$b1*data$t + data$b2*data$Y4 + 
    data$u1*data$t + data$u0 + data$e
  
  ##
  data$b0i = data$b0 + data$u0
  data$b1i = data$b1 + data$u1
  
  b0 = list()
  for(i in 1:nindividual){
    b0[i] =  data$b0i[data$ID == i][1]
  }
  b0 = unlist(b0)
  
  b1 = list()
  for(i in 1:nindividual){
    b1[i] =  data$b1i[data$ID == i][1]
  }
  b1 =unlist(b1)
  
  b2 = list()
  for(i in 1:nindividual){
    b2[i] =  data$b2[data$ID == i][1]
  }
  b2 =unlist(b2)
  
  sde = list()
  for(i in 1:nindividual){
    sde[i] =  sd(data$e[data$ID == i])
  }
  sde = unlist(sde)
  
  CPTY10 = list(
    coef = matrix(c(b0, b1, b2), nrow = 3, ncol = nindividual, byrow = TRUE,
                  dimnames = list(c("(Intercept)", "t", "Y4"), as.character(seq(0, nindividual-1)))),
    sd = sde
  )
  data <- data[, !(names(data) %in% c("u0","u1", "e", "b0","b1","b0i","b1i", "b2","b3"))]
  
  #== CPT for ID&t
  CPTID = matrix(rep(1/nindividual, nindividual), ncol = nindividual,
                 dimnames = list(NULL, as.character(seq(0, nindividual-1))))
  
  CPTt = list(
    coef = c("(Intercept)" = mean(data$t)),
    sd = sd(data$t))
  
  
  #== True MBN ==#
  # True DAG
  TMBN = model2network("[ID][t][Y1|ID:t][Y3|ID:t][Y2|Y1:ID:t][Y5|Y3:ID:t][Y6|Y1:Y3:ID:t][Y4|Y2:ID:t][Y8|Y5:ID:t][Y7|Y4:ID:t][Y9|Y3:Y4:Y5:ID:t][Y10|Y4:ID:t]")
  
  # True parameters 
  cgfit = custom.fit(TMBN, dist = list(ID =CPTID, t = CPTt, Y1 = CPTY1, Y2=CPTY2, Y3 = CPTY3, Y4 =CPTY4, Y5 = CPTY5, Y6 = CPTY6, Y7 = CPTY7, Y8 = CPTY8, Y9=CPTY9, Y10 = CPTY10))
  
  #== Structure learning ==#
  data$ID<- as.factor(data$ID)
  data$t<- as.numeric(data$t)
  
  all_nodes <- names(data)
  
  whitelist.arcs <- data.frame(
    from = rep(c("ID", "t"), each = length(setdiff(all_nodes, c("ID", "t")))),
    to = rep(setdiff(all_nodes, c("ID", "t")), 2)
  )
  
  blacklist.arcs <- data.frame(
    from = rep(setdiff(all_nodes, c("ID", "t")), 2),
    to = rep(c("ID", "t"), each = length(setdiff(all_nodes, c("ID", "t"))))
  )
  
  
  
  dag.lme <- hc(data, score = "custom-score", fun = scoreME, 
                args = list(group = "ID", t= "t"), 
                whitelist = whitelist.arcs, 
                blacklist = blacklist.arcs)
  
  dag.inla <- hc(data, score = "custom-score", fun = scoreME1, 
                 args = list(group = "ID", t= "t"), 
                 whitelist = whitelist.arcs, 
                 blacklist = blacklist.arcs)
  
  dag.inla1 <- hc(data, score = "custom-score", fun = scoreME2, 
                  args = list(group = "ID", t= "t"), 
                  whitelist = whitelist.arcs, 
                  blacklist = blacklist.arcs)
  dag.inla2 <- hc(data, score = "custom-score", fun = scoreME3, 
                  args = list(group = "ID", t= "t"), 
                  whitelist = whitelist.arcs, 
                  blacklist = blacklist.arcs)
  dag.inla3 <- hc(data, score = "custom-score", fun = scoreME4, 
                  args = list(group = "ID", t= "t"), 
                  whitelist = whitelist.arcs, 
                  blacklist = blacklist.arcs)
  

  TMBN_cp <- cpdag(TMBN)
  dag.lme_cp  <- cpdag(dag.lme)
  dag.inla_cp <- cpdag(dag.inla)
  dag.inla1_cp <- cpdag(dag.inla1)
  dag.inla2_cp <- cpdag(dag.inla2)
  dag.inla3_cp <- cpdag(dag.inla3)
  #==== Structure comparison with the TRUE DAG====#
  shd_lme[[k]] = shd(TMBN_cp, dag.lme_cp) 
  shd_inla[[k]] = shd(TMBN_cp, dag.inla_cp)
  shd_inla1[[k]] = shd(TMBN_cp, dag.inla1_cp) 
  shd_inla2[[k]] = shd(TMBN_cp, dag.inla2_cp) 
  shd_inla3[[k]] = shd(TMBN_cp, dag.inla3_cp) 
  
  ##==== parameter learning_lme====##
  
  local.distributions = 
    structure(vector(nnodes(dag.lme), mode = "list"), names = nodes(dag.lme))
  
  for (node in nodes(dag.lme)) {
    
    if (node == "ID") {
      local.distributions[["ID"]] = prop.table(table(data$ID))
      
    } else if (node == "t") {
      local.distributions[["t"]] = list(
        coef = c("(Intercept)" = mean(data$t)),
        sd = sd(data$t)
      )
      
    } else {
      local.distributions[[node]] = ldistME_lme(
        node = node, 
        parents = setdiff(parents(dag.lme, node), c("ID", "t")),
        group = "ID", 
        t = "t",
        data = data
      )
    }
  }
  
  bn.lme = custom.fit(dag.lme, local.distributions) 
  
  #==== Parameter learning INLA_ default Wishart prior =====#
  local.distributions = 
    structure(vector(nnodes(dag.inla), mode = "list"), names = nodes(dag.inla))
  
  for (node in nodes(dag.inla)) {
    if (node == "ID") {
      local.distributions[["ID"]] = prop.table(table(data$ID))
      
    } else if (node == "t") {
      local.distributions[["t"]] = list(
        coef = c("(Intercept)" = mean(data$t)),
        sd = sd(data$t)
      )
      
    } else {
      local.distributions[[node]] = ldistME_inla(
        node = node, 
        parents = setdiff(parents(dag.inla, node), c("ID", "t")),
        group = "ID", 
        t = "t",
        data = data
      )
    }
  }
  
  bn.inla = custom.fit(dag.inla, local.distributions) 
  
  #==== Parameter learning INLA_ user defind Wishart prior =====#
  local.distributions = 
    structure(vector(nnodes(dag.inla1), mode = "list"), names = nodes(dag.inla1))
  
  for (node in nodes(dag.inla1)) {
    if (node == "ID") {
      local.distributions[["ID"]] = prop.table(table(data$ID))
      
    } else if (node == "t") {
      local.distributions[["t"]] = list(
        coef = c("(Intercept)" = mean(data$t)),
        sd = sd(data$t)
      )
      
    } else {
      local.distributions[[node]] = ldistME_inla1(
        node = node, 
        parents = setdiff(parents(dag.inla1, node), c("ID", "t")),
        group = "ID", 
        t = "t",
        data = data
      )
    }
  }
  
  bn.inla1 = custom.fit(dag.inla1, local.distributions)
  
  #==== Parameter learning INLA_ default LKJ prior =====#
  local.distributions = 
    structure(vector(nnodes(dag.inla2), mode = "list"), names = nodes(dag.inla2))
  
  for (node in nodes(dag.inla2)) {
    if (node == "ID") {
      local.distributions[["ID"]] = prop.table(table(data$ID))
      
    } else if (node == "t") {
      local.distributions[["t"]] = list(
        coef = c("(Intercept)" = mean(data$t)),
        sd = sd(data$t)
      )
      
    } else {
      local.distributions[[node]] = ldistME_inla2(
        node = node, 
        parents = setdiff(parents(dag.inla2, node), c("ID", "t")),
        group = "ID", 
        t = "t",
        data = data
      )
    }
  }
  
  bn.inla2 = custom.fit(dag.inla2, local.distributions)
  
  #==== Parameter learning INLA_ user defind LKJ prior =====#
  local.distributions = 
    structure(vector(nnodes(dag.inla3), mode = "list"), names = nodes(dag.inla3))
  
  for (node in nodes(dag.inla3)) {
    if (node == "ID") {
      local.distributions[["ID"]] = prop.table(table(data$ID))
      
    } else if (node == "t") {
      local.distributions[["t"]] = list(
        coef = c("(Intercept)" = mean(data$t)),
        sd = sd(data$t)
      )
      
    } else {
      local.distributions[[node]] = ldistME_inla2(
        node = node, 
        parents = setdiff(parents(dag.inla3, node), c("ID", "t")),
        group = "ID", 
        t = "t",
        data = data
      )
    }
  }
  
  bn.inla3 = custom.fit(dag.inla3, local.distributions)
  
  
  #=== parameter comparison ====#
  KL_lme[[k]] <- log(KL(cgfit,bn.lme))
  KL_inla[[k]] <-log(KL(cgfit, bn.inla))
  KL_inla1[[k]] <-log(KL(cgfit, bn.inla1))
  KL_inla2[[k]] <-log(KL(cgfit, bn.inla2))
  KL_inla3[[k]] <-log(KL(cgfit, bn.inla3))
  
  #== Classification accuracy_ clustering variable
  F1_lme[[k]] = F1_score(bn.lme)
  F1_inla[[k]] = F1_score(bn.inla)
  F1_inla1[[k]] = F1_score(bn.inla1)
  F1_inla2[[k]] = F1_score(bn.inla2)
  F1_inla3[[k]] = F1_score(bn.inla3)
  
  #== classification accuracy  for the structure ==#
  F1_Slme[[k]] <- F1_score_structure(TMBN_cp, dag.lme_cp)
  F1_Sinla[[k]] <- F1_score_structure(TMBN_cp, dag.inla_cp)
  F1_Sinla1[[k]] <- F1_score_structure(TMBN_cp, dag.inla1_cp)
  F1_Sinla2[[k]] <- F1_score_structure(TMBN_cp, dag.inla2_cp)
  F1_Sinla3[[k]] <- F1_score_structure(TMBN_cp, dag.inla3_cp)
  
  F1_Slme_1[[k]] <- F1_score_structure1(TMBN, dag.lme)
  F1_Sinla_1[[k]] <- F1_score_structure1(TMBN, dag.inla)
  F1_Sinla1_1[[k]] <- F1_score_structure1(TMBN, dag.inla1)
  F1_Sinla2_1[[k]] <- F1_score_structure1(TMBN, dag.inla2)
  F1_Sinla3_1[[k]] <- F1_score_structure1(TMBN, dag.inla3)
  
}

shd_lme<- unlist(shd_lme)
shd_inla<- unlist(shd_inla)
shd_inla1<- unlist(shd_inla1)
shd_inla2<- unlist(shd_inla2)
shd_inla3<- unlist(shd_inla3)

KL_lme  <- unlist (KL_lme)
KL_inla <- unlist(KL_inla)
KL_inla1 <- unlist(KL_inla1)
KL_inla2 <- unlist(KL_inla2)
KL_inla3 <- unlist(KL_inla3)

F1_lme <-unlist(F1_lme)
F1_inla <- unlist(F1_inla)
F1_inla1<-unlist(F1_inla1)
F1_inla2<-unlist(F1_inla2)
F1_inla3<-unlist(F1_inla3)

F1_Slme <-unlist(F1_Slme)
F1_Sinla<-unlist(F1_Sinla)
F1_Sinla1<-unlist(F1_Sinla1)
F1_Sinla2<-unlist(F1_Sinla2)
F1_Sinla3<-unlist(F1_Sinla3)

F1_Slme_1 <-unlist(F1_Slme_1)
F1_Sinla_1<-unlist(F1_Sinla_1)
F1_Sinla1_1<-unlist(F1_Sinla1_1)
F1_Sinla2_1<-unlist(F1_Sinla2_1)
F1_Sinla3_1<-unlist(F1_Sinla3_1)

stats<-data.frame(shd_lme, shd_inla, shd_inla1, shd_inla2, shd_inla3, 
                  KL_lme, KL_inla, KL_inla1, KL_inla2, KL_inla3, 
                  F1_lme, F1_inla, F1_inla1, F1_inla2, F1_inla3, 
                  F1_Slme, F1_Sinla, F1_Sinla1, F1_Sinla2, F1_Sinla3,
                  F1_Slme_1, F1_Sinla_1, F1_Sinla1_1, F1_Sinla2_1, F1_Sinla3_1)

setwd("BN_INLA/Output/27April")
out.file2 <- sprintf(
  "Gauss_Long2_%s_%s.rds",
  Sys.Date(), format(Sys.time(), "%H%M%S")
)
saveRDS(stats, file = out.file2)
#======
