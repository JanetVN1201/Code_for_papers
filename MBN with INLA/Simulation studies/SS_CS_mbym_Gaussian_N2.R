##=== extending MBN using INLA for hirarchecal (multilevel) data (correlated random effects $m \times m$ random effects ====#
library(bnlearn)
library(igraph)
library(INLA)
library(gridExtra)
library(ggplot2)
library(lme4)
library(dplyr)
library(mvtnorm)
library(psych)
library(caret)
library(graphpcor)

nsims = 100
##==== Custom score function MLE ====##
fitME = function(node, parents, group, data) {
  
  rhs = paste(c("1", parents), collapse = " + ")
  formula = paste(node, "~", rhs, "+", paste0("(", rhs, " | ", group, ")"))
  
  model = try(lmer(formula, data = data, control = lmerControl(calc.derivs = FALSE)))
  
  if (is(model, "try-error"))
    return(NULL)
  
  return(model)
  
}#FITME

scoreME =  function(node, parents, data, args = list()) {
  
  if (node == args$group)
    return(0)
  
  # mixed effects models with performance tuning to speed up learning.
  require(lme4)
  environment(nloptwrap)$defaultControl$maxeval = 5 * 1e3
  environment(nloptwrap)$defaultControl$xtol_abs = 1e-6
  environment(nloptwrap)$defaultControl$ftol_abs = 1e-6
  
  local.distribution =
    fitME(node = node, parents = setdiff(parents, args$group),
          group = args$group, data = data)
  
  if (is.null(local.distribution))
    return(-Inf)
  
  return(- BIC(local.distribution))
  
}#SCOREME

##=== Custom score function INLA_ default prior ===#

fitME1 = function(node, parents, group, data) {
  nid <- length(unique(data[[group]]))  
  data$numid <- as.numeric(as.factor(data[[group]]))  
  
  if (length(parents) == 0) {
    formula = as.formula(paste(node, "~ 1 + f(numid, model = 'iidkd', n =", nid, ", order = 2)"))
  } else {
    
    rhs = paste("1", paste(parents, collapse = " + "), sep = " + ")  
    m = length(parents) + 1  
    
    slope_names <- paste0("slopeid", seq_along(parents)) 
    data[[slope_names[1]]] <- data$numid + nid  
    
    if (length(parents) > 1) {
      for (i in 2:length(parents)) {
        data[[slope_names[i]]] <- data[[slope_names[i - 1]]] + nid  
      }
    }
    
    slope_terms <- paste0("f(", slope_names, ", ", parents, ", copy = 'numid')", collapse = " + ")
    
    formula = as.formula(paste(node, "~", rhs, "+ f(numid, model = 'iidkd', n =", m * nid, 
                               ", order =", m, ") +", slope_terms))
  }
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
  
  # mixed effects models with performance tuning to speed up learning.
  require(INLA)
  
  local.distribution =
    fitME1(node = node, parents = setdiff(parents, args$group),
           group = args$group, data = data)
  
  if (is.null(local.distribution))
    return(-Inf)
  
  return(local.distribution$mlik[1])
  
  
}

#== score function INLA- user defind wishart prior ===#
fitME2 = function(node, parents, group, data) {
  nid <- length(unique(data[[group]]))  
  data$numid <- as.numeric(as.factor(data[[group]]))  # Convert group to numeric IDs
  
  if (length(parents) == 0) {
    formula <- as.formula(paste(
      node, "~ 1 +",
      "f(numid, model = 'iidkd', n =", nid,
      ", order = 2, hyper = list(theta1 = list(param = c(10,1,1,0))))"
    ))
  } else {
    
    rhs = paste("1", paste(parents, collapse = " + "), sep = " + ")  
    m = length(parents) + 1  
    
    slope_names <- paste0("slopeid", seq_along(parents)) 
    data[[slope_names[1]]] <- data$numid + nid  
    
    if (length(parents) > 1) {
      for (i in 2:length(parents)) {
        data[[slope_names[i]]] <- data[[slope_names[i - 1]]] + nid  
      }
    }
    
    slope_terms <- paste0("f(", slope_names, ", ", parents, ", copy = 'numid')", collapse = " + ")
    
    if (length(parents) == 1) {
      formula <- as.formula(paste(
        node, "~", rhs, "+",
        "f(numid, model = 'iidkd', n =", m * nid,
        ", order =", m,
        ", constr = FALSE, hyper = list(theta1 = list(param = c(10,1,1,0))))",
        "+", slope_terms
      ))
    } else if (length(parents) == 2) {
      formula <- as.formula(paste(
        node, "~", rhs, "+",
        "f(numid, model = 'iidkd', n =", m * nid,
        ", order =", m,
        ", constr = FALSE, hyper = list(theta1 = list(param = c(10,1,1,1,0,0,0))))",
        "+", slope_terms
      ))
    } else if (length(parents) == 3) {
      formula <- as.formula(paste(
        node, "~", rhs, "+",
        "f(numid, model = 'iidkd', n =", m * nid,
        ", order =", m,
        ", constr = FALSE, hyper = list(theta1 = list(param = c(6,1,1,1,1,0,0,0,0,0,0))))",
        "+", slope_terms
      ))
    } else {
      formula <- as.formula(paste(
        node, "~", rhs, "+",
        "f(numid, model = 'iidkd', n =", m * nid,
        ", order =", m,
        ", constr = FALSE, hyper = list(theta1 = list(param = c(6,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0))))",
        "+", slope_terms
      ))
    }
  }
  
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
  
  # mixed effects models with performance tuning to speed up learning.
  require(INLA)
  
  local.distribution =
    fitME2(node = node, parents = setdiff(parents, args$group),
           group = args$group, data = data)
  
  if (is.null(local.distribution))
    return(-Inf)
  
  return(local.distribution$mlik[1])
}


#== score function INLA for default LKJ prior ===# 
fitME3 <- function(node, parents, group, data) {
  data[[group]] <- as.numeric(data[[group]])
  data$rintercept <- rep(1, each = ncluster * csize)
  
  eta_values <- c(10, 20) 
  
  if (length(parents) == 0) {
    m <- 2
  } else {
    m <- length(parents) + 1
    #if (length(parents) > 2) eta_values <- c(10, 20)
  }
  
  cinla <- list(int.strategy = 'eb')
  cfam <- list(hyper = list())
  
  # Build the formula based on number of parents
  if (length(parents) == 0) {
    formula_text <- paste0(node, " ~ 1 + f(rintercept, model = cmodel, replicate = ", group, ")")
  } else {
    rhs <- paste("1", paste(parents, collapse = " + "), sep = " + ")
    slope_names <- paste0("rSlope", seq_along(parents))
    
    data[[slope_names[1]]] <- rep(2, each = ncluster * csize)
    if (length(parents) > 1) {
      for (i in 2:length(parents)) {
        data[[slope_names[i]]] <- rep(i + 1, each = ncluster * csize)
      }
    }
    slope_terms <- paste0("f(", slope_names, ", ", parents, ", copy = 'rintercept', replicate = ", group, ")", collapse = " + ")
    
    formula_text <- paste0(node, " ~ ", rhs,
                           " + f(rintercept, model = cmodel, replicate = ", group, ")",
                           " + ", slope_terms)
  }
  
  for (eta in eta_values) {
    #message("Trying eta = ", eta, " for node ", node)
    
    cmodel <- cgeneric(
      model = "LKJ",
      n = m,
      eta = eta
    )
    
    formula <- as.formula(formula_text)
    model <- try(
      {
        inla(formula, family = "gaussian", data = data,
             control.inla = cinla,
             control.family = cfam)
      },
      silent = TRUE
    )
    
    if (!inherits(model, "try-error")) {
      return(model)
    } else {
      #message("INLA crashed for node ", node, " with eta = ", eta, ". Trying next eta...")
    }
  }
  
  warning("All model attempts failed for node ", node)
  return(NULL)
}


scoreME3 =  function(node, parents, data, args = list()) {
  
  if (node == args$group)
    return(0)
  require(INLA)
  
  local.distribution =
    fitME3(node = node, parents = setdiff(parents, args$group),
           group = args$group, data = data)
  
  if (is.null(local.distribution))
    return(-Inf)
  
  return(local.distribution$mlik[1])
}

#== score function INLA, user defind LKJ prior, for comparison ===# 
fitME4 <- function(node, parents, group, data) {
  data[[group]] <- as.numeric(data[[group]])
  data$rintercept <- rep(1, each = ncluster * csize)
  
  eta_values <- c(8, 10, 20) 
  
  if (length(parents) == 0) {
    m <- 2
  } else {
    m <- length(parents) + 1
    if (length(parents) > 2) eta_values <- c(8, 12, 20)
  }
  
  cinla <- list(int.strategy = 'eb')
  cfam <- list(hyper = list())
  
  # Build the formula based on number of parents
  if (length(parents) == 0) {
    formula_text <- paste0(node, " ~ 1 + f(rintercept, model = cmodel, replicate = ", group, ")")
  } else {
    rhs <- paste("1", paste(parents, collapse = " + "), sep = " + ")
    slope_names <- paste0("rSlope", seq_along(parents))
    
    data[[slope_names[1]]] <- rep(2, each = ncluster * csize)
    if (length(parents) > 1) {
      for (i in 2:length(parents)) {
        data[[slope_names[i]]] <- rep(i + 1, each = ncluster * csize)
      }
    }
    slope_terms <- paste0("f(", slope_names, ", ", parents, ", copy = 'rintercept', replicate = ", group, ")", collapse = " + ")
    
    formula_text <- paste0(node, " ~ ", rhs,
                           " + f(rintercept, model = cmodel, replicate = ", group, ")",
                           " + ", slope_terms)
  }
  
  for (eta in eta_values) {
    #message("Trying eta = ", eta, " for node ", node)
    
    cmodel <- cgeneric(
      model = "LKJ",
      n = m,
      eta = eta
    )
    
    formula <- as.formula(formula_text)
    model <- try(
      {
        inla(formula, family = "gaussian", data = data,
             control.inla = cinla,
             control.family = cfam)
      },
      silent = TRUE
    )
    
    if (!inherits(model, "try-error")) {
      return(model)
    } else {
      #message("INLA crashed for node ", node, " with eta = ", eta, ". Trying next eta...")
    }
  }
  
  warning("All model attempts failed for node ", node)
  return(NULL)
}


scoreME4 =  function(node, parents, data, args = list()) {
  
  if (node == args$group)
    return(0)
  require(INLA)
  
  local.distribution =
    fitME4(node = node, parents = setdiff(parents, args$group),
           group = args$group, data = data)
  
  if (is.null(local.distribution))
    return(-Inf)
  
  return(local.distribution$mlik[1])
}


#=== Parameter learning function MLE ====#
ldistME_lme = function(node, parents, group, data) {
  
  require(lme4)
  
  environment(nloptwrap)$defaultControl$maxeval = 5 * 1e3
  environment(nloptwrap)$defaultControl$xtol_abs = 1e-6
  environment(nloptwrap)$defaultControl$ftol_abs = 1e-6
  
  model = fitME(node = node, parents = parents, group = group, data = data)
  ngroups = nlevels(data[[group]])  
  
  ldist = list(
    coef = array(0, dim = c(length(parents) + 1, ngroups), 
                 dimnames = list(c("(Intercept)", parents), NULL)),
    sd = rep(0, ngroups)
  )
  
  if (is.null(model)) {
    environment(nloptwrap)$defaultControl$maxeval = 1e5
    environment(nloptwrap)$defaultControl$xtol_abs = 1e-5
    environment(nloptwrap)$defaultControl$ftol_abs = 1e-5
    
    model = fitME(node = node, parents = parents, group = group, data = data)
    
    
    environment(nloptwrap)$defaultControl$maxeval = 5 * 1e3
    environment(nloptwrap)$defaultControl$xtol_abs = 1e-6
    environment(nloptwrap)$defaultControl$ftol_abs = 1e-6
    
    if (is.null(model)) {
      warning("Failed to fit node ", node, " with lme4.")
      return(ldist)  
    }
  }
  
  fixefs = fixef(model)
  for (i in seq_len(ncol(ldist$coef))) {
    ldist$coef[, i] = fixefs
  }
  
  ranefs = ranef(model)[[group]]
  group_levels = levels(data[[group]])
  for (i in seq_len(ngroups)) {
    ldist$coef[, i] = ldist$coef[, i] + as.numeric(ranefs[group_levels[i], ])
  }
  
  resids = resid(model)
  for (i in seq_len(ngroups)) {
    ldist$sd[i] = sd(resids[data[[group]] == group_levels[i]])
  }
  
  ldist$coef[is.na(ldist$coef)] = 0
  ldist$sd[is.na(ldist$sd)] = 0
  
  return(ldist)
}

#=== Parameter learning function INLA _ fitME1 ====#

ldistME_inla1 = function(node, parents, group, data) {
  require(INLA)
  
  model = fitME1(node = node, parents = parents, group = group, data = data)
  if (is.null(model)) {
    warning("Model is NULL for node ", node)
    return(NULL)
  }
  
  ngroups = nlevels(data[[group]])
  group_levels = levels(data[[group]])
  
  rand_vars = c("(Intercept)", parents)
  m = length(rand_vars)
  
  ldist = list(
    coef = array(0, dim = c(m, ngroups),
                 dimnames = list(rand_vars, NULL)),
    sd = rep(0, ngroups)
  )
  
  if (is.null(model)) {
    
    model = fitME1(node = node, parents = parents, group = group, data = data)
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
  
  for (i in seq_len(ngroups)) {
    for (name in coef_names) {
      if (name %in% rand_vars) {
        ldist$coef[name, i] = fixefs[name]
      }
    }
  }
  
  re_df = model$summary.random$numid
  if (nrow(re_df) != m * ngroups) {
    stop("Mismatch in number of random effects and expected size.")
  }
  
  re_df$effect = rep(rand_vars, each = ngroups)
  re_df$group = rep(group_levels, times = m)
  
  ranefs_wide = reshape(re_df[, c("group", "effect", "mode")],
                        timevar = "effect", idvar = "group", direction = "wide")
  colnames(ranefs_wide) = gsub("mode\\.", "", colnames(ranefs_wide))
  rownames(ranefs_wide) = ranefs_wide$group
  ranefs_wide$group = NULL
  
  for (i in seq_len(ngroups)) {
    for (eff in rand_vars) {
      if (eff %in% colnames(ranefs_wide)) {
        ldist$coef[eff, i] = ldist$coef[eff, i] + ranefs_wide[group_levels[i], eff]
      }
    }
  }
  
  pred = model$summary.fitted.values$mean
  resids = data[[node]] - pred
  for (i in seq_len(ngroups)) {
    ldist$sd[i] = sd(resids[data[[group]] == group_levels[i]])
  }
  
  return(ldist)
}

#=== Parameter learning function INLA _ fitME2 ====#

ldistME_inla2 = function(node, parents, group, data) {
  require(INLA)
  
  model = fitME2(node = node, parents = parents, group = group, data = data)
  if (is.null(model)) {
    warning("Model is NULL for node ", node)
    return(NULL)
  }
  
  ngroups = nlevels(data[[group]])
  group_levels = levels(data[[group]])
  
  rand_vars = c("(Intercept)", parents)
  m = length(rand_vars)
  
  ldist = list(
    coef = array(0, dim = c(m, ngroups),
                 dimnames = list(rand_vars, NULL)),
    sd = rep(0, ngroups)
  )
  
  if (is.null(model)) {
    
    model = fitME2(node = node, parents = parents, group = group, data = data)
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
  
  for (i in seq_len(ngroups)) {
    for (name in coef_names) {
      if (name %in% rand_vars) {
        ldist$coef[name, i] = fixefs[name]
      }
    }
  }
  
  re_df = model$summary.random$numid
  if (nrow(re_df) != m * ngroups) {
    stop("Mismatch in number of random effects and expected size.")
  }
  
  re_df$effect = rep(rand_vars, each = ngroups)
  re_df$group = rep(group_levels, times = m)
  
  ranefs_wide = reshape(re_df[, c("group", "effect", "mode")],
                        timevar = "effect", idvar = "group", direction = "wide")
  colnames(ranefs_wide) = gsub("mode\\.", "", colnames(ranefs_wide))
  rownames(ranefs_wide) = ranefs_wide$group
  ranefs_wide$group = NULL
  
  for (i in seq_len(ngroups)) {
    for (eff in rand_vars) {
      if (eff %in% colnames(ranefs_wide)) {
        ldist$coef[eff, i] = ldist$coef[eff, i] + ranefs_wide[group_levels[i], eff]
      }
    }
  }
  
  pred = model$summary.fitted.values$mean
  resids = data[[node]] - pred
  for (i in seq_len(ngroups)) {
    ldist$sd[i] = sd(resids[data[[group]] == group_levels[i]])
  }
  
  return(ldist)
}

#=== Parameter learning function INLA _ fitME3 ====#
ldistME_inla3 = function(node, parents, group, data) {
  require(INLA)
  
  model = fitME3(node = node, parents = parents, group = group, data = data)
  if (is.null(model)) {
    warning("Model is NULL for node ", node)
    return(NULL)
  }
  
  ngroups = nlevels(data[[group]])
  group_levels = levels(data[[group]])
  
  rand_vars = c("(Intercept)", parents)
  m = length(rand_vars)
  
  ldist = list(
    coef = array(0, dim = c(m, ngroups),
                 dimnames = list(rand_vars, NULL)),
    sd = rep(0, ngroups)
  )
  
  if (is.null(model)) {
    
    model = fitME3(node = node, parents = parents, group = group, data = data)
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
  
  for (i in seq_len(ngroups)) {
    for (name in coef_names) {
      if (name %in% rand_vars) {
        ldist$coef[name, i] = fixefs[name]
      }
    }
  }
  if (length(rand_vars)==1){
    re_df = matrix(model$summary.random$rintercept$mode,
                   byrow = TRUE, nrow = ncluster)
    re_df <- re_df[, 1, drop = FALSE]
    colnames(re_df)= rand_vars
    re_df <- as.data.frame(re_df)
  } else{
    re_df = matrix(model$summary.random$rintercept$mode,
                   byrow = TRUE, nrow = ncluster)
    colnames(re_df)= rand_vars
    re_df <- as.data.frame(re_df)
  }
  
  for (i in seq_len(ngroups)) {
    for (eff in rand_vars) {
      if (eff %in% colnames(re_df)) {
        ldist$coef[eff, i] = ldist$coef[eff, i] + re_df[group_levels[i], eff]
      }
    }
  }
  
  pred = model$summary.fitted.values$mean
  resids = data[[node]] - pred
  for (i in seq_len(ngroups)) {
    ldist$sd[i] = sd(resids[data[[group]] == group_levels[i]])
  }
  
  return(ldist)
}

#=== Parameter learning function INLA _ fitME4 ====#

ldistME_inla4 = function(node, parents, group, data) {
  require(INLA)
  
  model = fitME4(node = node, parents = parents, group = group, data = data)
  if (is.null(model)) {
    warning("Model is NULL for node ", node)
    return(NULL)
  }
  
  ngroups = nlevels(data[[group]])
  group_levels = levels(data[[group]])
  
  rand_vars = c("(Intercept)", parents)
  m = length(rand_vars)
  
  ldist = list(
    coef = array(0, dim = c(m, ngroups),
                 dimnames = list(rand_vars, NULL)),
    sd = rep(0, ngroups)
  )
  
  if (is.null(model)) {
    
    model = fitME4(node = node, parents = parents, group = group, data = data)
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
  
  for (i in seq_len(ngroups)) {
    for (name in coef_names) {
      if (name %in% rand_vars) {
        ldist$coef[name, i] = fixefs[name]
      }
    }
  }
  
  if (length(rand_vars)==1){
    re_df = matrix(model$summary.random$rintercept$mode,
                   byrow = TRUE, nrow = ncluster)
    re_df <- re_df[, 1, drop = FALSE]
    colnames(re_df)= rand_vars
    re_df <- as.data.frame(re_df)
  } else{
    re_df = matrix(model$summary.random$rintercept$mode,
                   byrow = TRUE, nrow = ncluster)
    colnames(re_df)= rand_vars
    re_df <- as.data.frame(re_df)
  }
  
  for (i in seq_len(ngroups)) {
    for (eff in rand_vars) {
      if (eff %in% colnames(re_df)) {
        ldist$coef[eff, i] = ldist$coef[eff, i] + re_df[group_levels[i], eff]
      }
    }
  }
  
  pred = model$summary.fitted.values$mean
  resids = data[[node]] - pred
  for (i in seq_len(ngroups)) {
    ldist$sd[i] = sd(resids[data[[group]] == group_levels[i]])
  }
  
  return(ldist)
}

#== KLD to compaire parametrs of mle and INLA with the true parameters ===#
KL = function(P, Q) {
  Pdag = remove.node(bn.net(P), "G")
  Qdag = remove.node(bn.net(Q), "G")
  klvalue = 0  
  
  for (group in seq_along(coef(P$G))) {  
    
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
  
  pred = predict(fitted, node = "G", data = data, method = "bayes-lw")
  obs = data$G
  
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


#====Simulation study for m by m correlated RE with different Wishart priors=====#

#== Random DAG used to simulate the data ===#
TMBN = model2network("[G][Y1|G][Y2|G][Y8|G][Y3|Y1:G][Y4|Y1:G][Y5|Y2:G][Y6|Y1:Y5:G][Y9|Y1:Y4:Y8:G][Y7|Y5:Y6:G][Y10|Y4:Y6:G]")

shd_lme = list() 
shd_inla = list() 
shd_inla1 = list()  
shd_inla2 = list() 
shd_inla3 = list() 

KL_lme <- list()
KL_inla <-list()
KL_inla1 <-list()
KL_inla2 <-list()
KL_inla3 <-list()

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
  ncluster = 30
  csize = 5
  N = ncluster*csize
  data = data.frame("ID" = rep(1:N,
                               each = 1),
                    "G" = rep(1:ncluster,
                              each = csize))
  #== Y1|G
  
  
  sigmau01 = 0.35 
  u0 = rnorm(n = ncluster,
              mean = 0,
              sd = sigmau01)
  
  data$u0 = rep(u0,
                 each = 5, times=1)
  
  sigmae = sqrt(0.5)
  data$e = rnorm(N, 0, sigmae)
  
  data$b0 = 0.72
  data$Y1 = data$b0 +  data$u0 + data$e
  
  ##
  data$b0i = data$b0 + data$u0
  
  b0 = list()
  for(i in 1:ncluster){
    b0[i] =  data$b0i[data$G == i][1]
  }
  b0 = unlist(b0)
  
  sde = list()
  for(i in 1:ncluster){
    sde[i] =  sd(data$e[data$G == i])
  }
  sde = unlist(sde)
  
  
  CPTY1 = list(
    coef = matrix(b0, ncol = ncluster, byrow = TRUE,
                  dimnames = list("(Intercept)", as.character(seq(0, ncluster-1)))),
    sd = sde
  )
  data <- data[, !(names(data) %in% c("u0","u1", "e", "b0","b0i"))]
  
  #== Y2|Y1:G
  
  m= 2
  rho <- -0.4
  Sigma <- matrix(NA, m, m)
  diag(Sigma) <- c(0.4,0.35)
  for(i in 1:m) {
    for (j in 1:m) {
      if (i != j) {
        Sigma[i, j] <- rho^abs(i-j) * sqrt(Sigma[i, i] * Sigma[j, j])
      }
    }
  }
  random_effects <- rmvnorm(n = ncluster,sigma = Sigma)
  
  data$u0 = rep(random_effects[, 1],
                each = 5, times=1)
  
  data$u1 = rep(random_effects[, 2],
                each = 5, times=1)
  
  sigmae = sqrt(0.5) # 0.6
  data$e = rnorm(N, 0, sigmae)
  
  data$b0 = 1.3
  data$b1 = 0.3
  
  data$Y2 = data$b0 + data$b1*data$Y1 + data$u1*data$Y1 + data$u0 + data$e
  
  ##
  data$b0i = data$b0 + data$u0
  
  b0 = list()
  for(i in 1:ncluster){
    b0[i] =  data$b0i[data$G == i][1]
  }
  b0 = unlist(b0)
  
  data$b1i = data$b1 + data$u1
  b1 = list()
  for(i in 1:ncluster){
    b1[i] =  data$b1i[data$G == i][1]
  }
  b1 = unlist(b1)
  
  sde = list()
  for(i in 1:ncluster){
    sde[i] =  sd(data$e[data$G == i])
  }
  sde = unlist(sde)
  
  
  CPTY2 = list(
    coef = matrix(c(b0, b1), nrow = 2, ncol = ncluster, byrow = TRUE,
                  dimnames = list(c("(Intercept)", "Y1"), as.character(seq(0, ncluster-1)))),
    sd = sde
  )
  data <- data[, !(names(data) %in% c("u0","u1", "e", "b0","b0i", "b1i", "b1"))]
  
  #== Y3|G
  sigmau0 = 0.2 # 0.22
  u0 = rnorm(n = ncluster,
              mean = 0,
              sd = sigmau0)
  
  data$u0 = rep(u0,
                 each = 5, times=1)
  
  sigmae = sqrt(0.52)
  data$e = rnorm(N, 0, sigmae)
  
  data$b0 = -1.03
  data$Y3 = data$b0 +  data$u0 + data$e
  
  ##
  data$b0i = data$b0 + data$u0
  
  b0 = list()
  for(i in 1:ncluster){
    b0[i] =  data$b0i[data$G == i][1]
  }
  b0 = unlist(b0)
  
  sde = list()
  for(i in 1:ncluster){
    sde[i] =  sd(data$e[data$G == i])
  }
  sde = unlist(sde)
  
  
  CPTY3 = list(
    coef = matrix(b0, ncol = ncluster, byrow = TRUE,
                  dimnames = list("(Intercept)", as.character(seq(0, ncluster-1)))),
    sd = sde
  )
  data <- data[, !(names(data) %in% c("u0","u1", "e", "b0","b0i"))]
  
  #===Y4|Y2:G
  m= 2
  rho <- 0.54
  Sigma <- matrix(NA, m, m)
  diag(Sigma) <- c(0.48, 0.4)
  for(i in 1:m) {
    for (j in 1:m) {
      if (i != j) {
        Sigma[i, j] <- rho^abs(i-j) * sqrt(Sigma[i, i] * Sigma[j, j])
      }
    }
  }
  random_effects <- rmvnorm(n = ncluster,  sigma = Sigma)
  
  
  data$u0 = rep(random_effects[, 1],
                each = 5, times=1)
  
  data$u1 = rep(random_effects[, 2],
                each = 5, times=1)
  
  sigmae = sqrt(0.54) # 0.35
  data$e = rnorm(N, 0, sigmae)
  
  data$b0 = -0.32
  data$b1 = 0.56
  
  data$Y4 = data$b0 + data$b1*data$Y2 + data$u1*data$Y2 + data$u0 + data$e
  
  ##
  data$b0i = data$b0 + data$u0
  
  b0 = list()
  for(i in 1:ncluster){
    b0[i] =  data$b0i[data$G == i][1]
  }
  b0 = unlist(b0)
  
  data$b1i = data$b1 + data$u1
  b1 = list()
  for(i in 1:ncluster){
    b1[i] =  data$b1i[data$G == i][1]
  }
  b1 = unlist(b1)
  
  sde = list()
  for(i in 1:ncluster){
    sde[i] =  sd(data$e[data$G == i])
  }
  sde = unlist(sde)
  
  
  CPTY4 = list(
    coef = matrix(c(b0, b1), nrow = 2, ncol = ncluster, byrow = TRUE,
                  dimnames = list(c("(Intercept)", "Y2"), as.character(seq(0, ncluster-1)))),
    sd = sde
  )
  
  data <- data[, !(names(data) %in% c("u0","u1", "e", "b0","b0i", "b1i", "b1"))]
  
  #== Y5|Y3:G 
  m= 2
  rho <- 0.4
  Sigma <- matrix(NA, m, m)
  diag(Sigma) <- c(0.52, 0.4)
  for(i in 1:m) {
    for (j in 1:m) {
      if (i != j) {
        Sigma[i, j] <- rho^abs(i-j) * sqrt(Sigma[i, i] * Sigma[j, j])
      }
    }
  }
  random_effects <- rmvnorm(n = ncluster,  sigma = Sigma)
  
  
  
  data$u0 = rep(random_effects[, 1],
                each = 5, times=1)
  
  data$u1 = rep(random_effects[, 2],
                each = 5, times=1)
  
  sigmae = sqrt(0.450) #0.95
  data$e = rnorm(N, 0, sigmae)
  
  data$b0 = 2.63
  data$b1 = 1.19
  
  data$Y5 = data$b0 + data$b1*data$Y3 + data$u1*data$Y3 + data$u0 + data$e
  
  ##
  data$b0i = data$b0 + data$u0
  
  b0 = list()
  for(i in 1:ncluster){
    b0[i] =  data$b0i[data$G == i][1]
  }
  b0 = unlist(b0)
  
  data$b1i = data$b1 + data$u1
  b1 = list()
  for(i in 1:ncluster){
    b1[i] =  data$b1i[data$G == i][1]
  }
  b1 = unlist(b1)
  
  sde = list()
  for(i in 1:ncluster){
    sde[i] =  sd(data$e[data$G == i])
  }
  sde = unlist(sde)
  
  
  CPTY5 = list(
    coef = matrix(c(b0, b1), nrow = 2, ncol = ncluster, byrow = TRUE,
                  dimnames = list(c("(Intercept)", "Y3"), as.character(seq(0, ncluster-1)))),
    sd = sde
  )
  
  data <- data[, !(names(data) %in% c("u0","u1", "e", "b0","b0i", "b1i", "b1"))]
  
  #== Y6|Y1:Y3:G
  m= 3
  rho <- 0.3
  Sigma <-  matrix(NA, m, m)
  diag(Sigma) <- c(0.54, 0.43, 0.34)
  for(i in 1:m) {
    for (j in 1:m) {
      if (i != j) {
        Sigma[i, j] <- rho^abs(i-j) * sqrt(Sigma[i, i] * Sigma[j, j])
      }
    }
  }
  random_effects <- rmvnorm(n = ncluster,  sigma = Sigma)
  
  data$u0 = rep(random_effects[, 1],
                each = 5, times=1)
  
  data$u1 = rep(random_effects[, 2],
                each = 5, times=1)
  data$u2 = rep(random_effects[, 3],
                each = 5, times=1)
  
  sigmae = sqrt(0.4) #0.4
  data$e = rnorm(N, 0, sigmae)
  
  data$b0 = 1.03
  data$b1 = -0.27
  data$b2 = 0.77
  
  data$Y6 = data$b0 + data$b1*data$Y1 +data$b2*data$Y3 + data$u1*data$Y1 + 
    data$u2*data$Y3 + data$u0 + data$e
  
  ##
  data$b0i = data$b0 + data$u0
  
  b0 = list()
  for(i in 1:ncluster){
    b0[i] =  data$b0i[data$G == i][1]
  }
  b0 = unlist(b0)
  
  data$b1i = data$b1 + data$u1
  b1 = list()
  for(i in 1:ncluster){
    b1[i] =  data$b1i[data$G == i][1]
  }
  b1 = unlist(b1)
  
  data$b2i = data$b2 + data$u2
  b2 = list()
  for(i in 1:ncluster){
    b2[i] =  data$b2i[data$G == i][1]
  }
  b2 = unlist(b2)
  
  sde = list()
  for(i in 1:ncluster){
    sde[i] =  sd(data$e[data$G == i])
  }
  sde = unlist(sde)
  
  CPTY6 = list(
    coef = matrix(c(b0, b1, b2), nrow = 3, ncol = ncluster, byrow = TRUE,
                  dimnames = list(c("(Intercept)", "Y1","Y3"), as.character(seq(0, ncluster-1)))),
    sd = sde
  )
  
  data <- data[, !(names(data) %in% c("u0","u1","u2", "e", "b0","b0i", "b1i", "b1", "b2i", "b2"))]
  
  #=== Y7|Y4:G
  
  m= 2
  rho <- 0.4
  Sigma <- matrix(NA, m, m)
  diag(Sigma) <- c(0.54, 0.41)
  for(i in 1:m) {
    for (j in 1:m) {
      if (i != j) {
        Sigma[i, j] <- rho^abs(i-j) * sqrt(Sigma[i, i] * Sigma[j, j])
      }
    }
  }
  random_effects <- rmvnorm(n = ncluster,  sigma = Sigma)
  
  data$u0 = rep(random_effects[, 1],
                each = 5, times=1)
  
  data$u1 = rep(random_effects[, 2],
                each = 5, times=1)
  
  sigmae = sqrt(0.45) # 0.6
  data$e = rnorm(N, 0, sigmae)
  
  data$b0 = 2.2
  data$b1 = 0.98
  
  data$Y7 = data$b0 + data$b1*data$Y4 + data$u1*data$Y4 + data$u0 + data$e
  
  ##
  data$b0i = data$b0 + data$u0
  
  b0 = list()
  for(i in 1:ncluster){
    b0[i] =  data$b0i[data$G == i][1]
  }
  b0 = unlist(b0)
  
  data$b1i = data$b1 + data$u1
  b1 = list()
  for(i in 1:ncluster){
    b1[i] =  data$b1i[data$G == i][1]
  }
  b1 = unlist(b1)
  
  sde = list()
  for(i in 1:ncluster){
    sde[i] =  sd(data$e[data$G == i])
  }
  sde = unlist(sde)
  
  
  CPTY7 = list(
    coef = matrix(c(b0, b1), nrow = 2, ncol = ncluster, byrow = TRUE,
                  dimnames = list(c("(Intercept)", "Y4"), as.character(seq(0, ncluster-1)))),
    sd = sde
  )
  data <- data[, !(names(data) %in% c("u0","u1", "e", "b0","b0i", "b1i", "b1"))]
  
  #== Y8|Y5:G
  
  m= 2
  rho <- 0.52
  Sigma <- matrix(NA, m, m)
  diag(Sigma) <- c(0.45, 0.38)
  for(i in 1:m) {
    for (j in 1:m) {
      if (i != j) {
        Sigma[i, j] <- rho^abs(i-j) * sqrt(Sigma[i, i] * Sigma[j, j])
      }
    }
  }
  random_effects <- rmvnorm(n = ncluster,  sigma = Sigma)
  
  data$u0 = rep(random_effects[, 1],
                each = 5, times=1)
  
  data$u1 = rep(random_effects[, 2],
                each = 5, times=1)
  
  sigmae = sqrt(0.45) # 0.6
  data$e = rnorm(N, 0, sigmae)
  
  data$b0 = -0.5
  data$b1 = 0.62
  
  data$Y8 = data$b0 + data$b1*data$Y5 + data$u1*data$Y5 + data$u0 + data$e
  
  ##
  data$b0i = data$b0 + data$u0
  
  b0 = list()
  for(i in 1:ncluster){
    b0[i] =  data$b0i[data$G == i][1]
  }
  b0 = unlist(b0)
  
  data$b1i = data$b1 + data$u1
  b1 = list()
  for(i in 1:ncluster){
    b1[i] =  data$b1i[data$G == i][1]
  }
  b1 = unlist(b1)
  
  sde = list()
  for(i in 1:ncluster){
    sde[i] =  sd(data$e[data$G == i])
  }
  sde = unlist(sde)
  
  
  CPTY8 = list(
    coef = matrix(c(b0, b1), nrow = 2, ncol = ncluster, byrow = TRUE,
                  dimnames = list(c("(Intercept)", "Y5"), as.character(seq(0, ncluster-1)))),
    sd = sde
  )
  data <- data[, !(names(data) %in% c("u0","u1", "e", "b0","b0i", "b1i", "b1"))]
  
  #== Y9|Y3:Y4:Y5:G
  m=4
  rho <- 0.42
  Sigma <- matrix(NA, m, m)
  diag(Sigma) <- c(0.6, 0.52, 0.42, 0.35)
  for(i in 1:m) {
    for (j in 1:m) {
      if (i != j) {
        Sigma[i, j] <- rho^abs(i-j) * sqrt(Sigma[i, i] * Sigma[j, j])
      }
    }
  }
  random_effects <- rmvnorm(n = ncluster,  sigma = Sigma)
  
  data$u0 = rep(random_effects[, 1],
                each = 5, times=1)
  
  data$u1 = rep(random_effects[, 2],
                each = 5, times=1)
  data$u2 = rep(random_effects[, 3],
                each = 5, times=1)
  
  data$u3 = rep(random_effects[, 4],
                each = 5, times=1)
  sigmae = sqrt(0.55) # 0.85
  data$e = rnorm(N, 0, sigmae)
  
  data$b0 = 0.97
  data$b1 = 0.27
  data$b2 = 1.31
  data$b3 = -0.34
  
  data$Y9 = data$b0 + data$b1*data$Y3 +data$b2*data$Y4 +data$b3*data$Y5 + 
    data$u1*data$Y3 + data$u2*data$Y4 + data$u3*data$Y5 + 
    data$u0 + data$e
  
  ##
  data$b0i = data$b0 + data$u0
  
  b0 = list()
  for(i in 1:ncluster){
    b0[i] =  data$b0i[data$G == i][1]
  }
  b0 = unlist(b0)
  
  data$b1i = data$b1 + data$u1
  b1 = list()
  for(i in 1:ncluster){
    b1[i] =  data$b1i[data$G == i][1]
  }
  b1 = unlist(b1)
  
  data$b2i = data$b2 + data$u2
  b2 = list()
  for(i in 1:ncluster){
    b2[i] =  data$b2i[data$G == i][1]
  }
  b2 = unlist(b2)
  
  data$b3i = data$b3 + data$u3
  b3 = list()
  for(i in 1:ncluster){
    b3[i] =  data$b3i[data$G == i][1]
  }
  b3 = unlist(b3)
  
  sde = list()
  for(i in 1:ncluster){
    sde[i] =  sd(data$e[data$G == i])
  }
  sde = unlist(sde)
  
  CPTY9 = list(
    coef = matrix(c(b0, b1, b2, b3), ncol = ncluster, nrow = 4, byrow = TRUE,
                  dimnames = list(c("(Intercept)", "Y3", "Y4","Y5"), as.character(seq(0, ncluster-1)))),
    sd = sde
  )
  data <- data[, !(names(data) %in% c("u0","u1","u2","u3", "e", "b0","b0i", "b1i", "b1", "b2i", "b2","b3i","b3"))]
  
  #== Y10|Y4:G
  m= 2
  rho <- 0.44
  Sigma <- matrix(NA, m, m)
  diag(Sigma) <- c(0.58, 0.42)
  for(i in 1:m) {
    for (j in 1:m) {
      if (i != j) {
        Sigma[i, j] <- rho^abs(i-j) * sqrt(Sigma[i, i] * Sigma[j, j])
      }
    }
  }
  random_effects <- rmvnorm(n = ncluster,  sigma = Sigma)
  
  data$u0 = rep(random_effects[, 1],
                each = 5, times=1)
  
  data$u1 = rep(random_effects[, 2],
                each = 5, times=1)
  
  sigmae = sqrt(0.52) # 0.6
  data$e = rnorm(N, 0, sigmae)
  
  data$b0 = 1.92
  data$b1 = 1.19
  
  data$Y10 = data$b0 + data$b1*data$Y4 + data$u1*data$Y4 + data$u0 + data$e
  
  ##
  data$b0i = data$b0 + data$u0
  
  b0 = list()
  for(i in 1:ncluster){
    b0[i] =  data$b0i[data$G == i][1]
  }
  b0 = unlist(b0)
  
  data$b1i = data$b1 + data$u1
  b1 = list()
  for(i in 1:ncluster){
    b1[i] =  data$b1i[data$G == i][1]
  }
  b1 = unlist(b1)
  
  sde = list()
  for(i in 1:ncluster){
    sde[i] =  sd(data$e[data$G == i])
  }
  sde = unlist(sde)
  
  
  CPTY10 = list(
    coef = matrix(c(b0, b1), nrow = 2, ncol = ncluster, byrow = TRUE,
                  dimnames = list(c("(Intercept)", "Y4"), as.character(seq(0, ncluster-1)))),
    sd = sde
  )
  data <- data[, !(names(data) %in% c("u0","u1", "e", "b0","b0i", "b1i", "b1", "ID"))]
  
  #== CPT for G
  CPTG = matrix(rep(1/ncluster, ncluster), ncol = ncluster,
                dimnames = list(NULL, as.character(seq(0, ncluster-1))))
  
  # === True DAG and True Parameters
  #True DAG
  TMBN = model2network("[G][Y1|G][Y3|G][Y2|Y1:G][Y5|Y3:G][Y6|Y1:Y3:G][Y4|Y2:G][Y8|Y5:G][Y7|Y4:G][Y9|Y3:Y4:Y5:G][Y10|Y4:G]")
  
  #True parameters
  cgfit = custom.fit(TMBN, dist = list(G =CPTG, Y1 = CPTY1, Y2=CPTY2, Y3 = CPTY3, Y4 =CPTY4, Y5 = CPTY5, Y6 = CPTY6, Y7 = CPTY7, Y8 = CPTY8, Y9=CPTY9, Y10 = CPTY10))
  
  #== Structure learning ==#
  data$G<- as.factor(data$G)
  inbound.arcs = data.frame(from = setdiff(names(data), "G"), to = "G")
  outbound.arcs = data.frame(from = "G", to = setdiff(names(data), "G"))
  
  
  dag.lme = hc(data, score = "custom-score", fun = scoreME, args = list(group = "G"),
               whitelist = outbound.arcs)
  dag.inla = hc(data, score = "custom-score", fun = scoreME1, args = list(group = "G"),
                whitelist = outbound.arcs)
  dag.inla1 = hc(data, score = "custom-score", fun = scoreME2, args = list(group = "G"),
                 whitelist = outbound.arcs)
  dag.inla2 = hc(data, score = "custom-score", fun = scoreME3, args = list(group = "G"),
                 whitelist = outbound.arcs)
  dag.inla3 = hc(data, score = "custom-score", fun = scoreME4, args = list(group = "G"),
                 whitelist = outbound.arcs)
  
  TMBN_cp <-cpdag(TMBN)
  dag.lme_cp<-cpdag(dag.lme)
  dag.inla_cp<-cpdag(dag.inla)
  dag.inla1_cp<-cpdag(dag.inla1)
  dag.inla2_cp<-cpdag(dag.inla2)
  dag.inla3_cp<-cpdag(dag.inla3)
  
  shd_lme[[k]] = shd(TMBN_cp, dag.lme_cp) 
  shd_inla[[k]] = shd(TMBN_cp, dag.inla_cp)
  shd_inla1[[k]] = shd(TMBN_cp, dag.inla1_cp) 
  shd_inla2[[k]] = shd(TMBN_cp, dag.inla2_cp) 
  shd_inla3[[k]] = shd(TMBN_cp, dag.inla3_cp) 
  
  ##==== parameter learning_lme====##
  local.distributions =
    structure(vector(nnodes(dag.lme), mode = "list"), names = nodes(dag.lme))
  for (node in nodes(dag.lme)) {
    
    if (node == "G")
      local.distributions[["G"]] = prop.table(table(data$G))
    else
      local.distributions[[node]] =
        ldistME_lme(node = node, parents = setdiff(parents(dag.lme, node), "G"),
                    group = "G", data = data)
    
  }#FOR
  bn.lme = custom.fit(dag.lme, local.distributions)
  
  ##==== parameter learning_ INLA====##
  #==  INLA default wishart
  local.distributions =
    structure(vector(nnodes(dag.inla), mode = "list"), names = nodes(dag.inla))
  for (node in nodes(dag.inla)) {
    
    if (node == "G")
      local.distributions[["G"]] = prop.table(table(data$G))
    else
      local.distributions[[node]] =
        ldistME_inla1(node = node, parents = setdiff(parents(dag.inla, node), "G"),
                      group = "G", data = data)
    
  }#FOR
  bn.inla = custom.fit(dag.inla, local.distributions)
  
  #==  INLA1 user defind wishart
  local.distributions =
    structure(vector(nnodes(dag.inla1), mode = "list"), names = nodes(dag.inla1))
  for (node in nodes(dag.inla1)) {
    
    if (node == "G")
      local.distributions[["G"]] = prop.table(table(data$G))
    else
      local.distributions[[node]] =
        ldistME_inla2(node = node, parents = setdiff(parents(dag.inla1, node), "G"),
                      group = "G", data = data)
    
  }#FOR
  bn.inla1 = custom.fit(dag.inla1, local.distributions)
  
  #==  INLA2 default lkj
  local.distributions =
    structure(vector(nnodes(dag.inla2), mode = "list"), names = nodes(dag.inla2))
  for (node in nodes(dag.inla2)) {
    
    if (node == "G")
      local.distributions[["G"]] = prop.table(table(data$G))
    else
      local.distributions[[node]] =
        ldistME_inla3(node = node, parents = setdiff(parents(dag.inla2, node), "G"),
                      group = "G", data = data)
    
  }#FOR
  bn.inla2 = custom.fit(dag.inla2, local.distributions)
  
  #==  INLA2 user defined lkj
  local.distributions =
    structure(vector(nnodes(dag.inla3), mode = "list"), names = nodes(dag.inla3))
  for (node in nodes(dag.inla3)) {
    
    if (node == "G")
      local.distributions[["G"]] = prop.table(table(data$G))
    else
      local.distributions[[node]] =
        ldistME_inla4(node = node, parents = setdiff(parents(dag.inla3, node), "G"),
                      group = "G", data = data)
    
  }#FOR
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

shd_lme = unlist(shd_lme) 
shd_inla = unlist(shd_inla) 
shd_inla1 = unlist(shd_inla1)  
shd_inla2 = unlist(shd_inla2)
shd_inla3 = unlist(shd_inla3)

KL_lme <- unlist(KL_lme)
KL_inla <-unlist(KL_inla)
KL_inla1 <-unlist(KL_inla1)
KL_inla2 <-unlist(KL_inla2)
KL_inla3 <-unlist(KL_inla3)

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

stats_1<-data.frame(shd_lme, shd_inla, shd_inla1, shd_inla2, shd_inla3,
                    KL_lme, KL_inla, KL_inla1, KL_inla2, KL_inla3, 
                    F1_lme, F1_inla, F1_inla1, F1_inla2, F1_inla3,
                    F1_Slme, F1_Sinla, F1_Sinla1, F1_Sinla2, F1_Sinla3,
                    F1_Slme_1, F1_Sinla_1, F1_Sinla1_1, F1_Sinla2_1, F1_Sinla3_1)

setwd("BN_INLA/Output/27April")
out.file2 <- sprintf(
  "Gauss_MbyM2_%s_%s.rds",
  Sys.Date(), format(Sys.time(), "%H%M%S")
)
saveRDS(stats_1, file = out.file2)

