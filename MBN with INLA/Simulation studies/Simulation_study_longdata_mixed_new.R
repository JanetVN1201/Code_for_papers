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
#=== Custom score function MLE ====#

fitME = function(node, parents, group, t, data) {
  
  rhs = paste(c("1", setdiff(parents, t), t), collapse = " + ")  
  formula = as.formula(paste(node, "~", rhs, "+ (1 +", t, "|", group, ")"))
  
  is_binary = all(data[[node]] %in% c(0, 1), na.rm = TRUE)
  
  if (is_binary) {
    model = try(glmer(formula, data = data, family = binomial, 
                      control = glmerControl(optimizer = "bobyqa")))
  } else {
    data$t = data$t+1
    model = try(lmer(formula, data = data, 
                     control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE)))
  }
  
  
  if (is(model, "try-error"))
    return(NULL)
  
  return(model)
  
} #FITME


scoreME = function(node, parents, data, args = list()) {
  
  if (node == args$group)
    return(0)
  if (node == args$t)
    return(0)

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

#== custom score function INLA2 default Wishart prior

fitME1 = function(node, parents, group, t, data) {
  
  data$Y1 <- as.numeric(as.character(data$Y1))
  data$Y2 <- as.numeric(as.character(data$Y2))
  data$Y3 <- as.numeric(as.character(data$Y3))
  data$Y4 <- as.numeric(as.character(data$Y4))
  nid <- length(unique(data[[group]]))  
  data$numid <- as.numeric(as.factor(data[[group]]))  
  data$slopeid <- data$numid +nid

  rhs = paste(c("1", setdiff(parents, t), t), collapse = " + ")
  
  formula = as.formula(paste0(node, "~", rhs, 
                              "+ f(numid, model = 'iidkd', n =", 2 * nid, ", order = 2)",
                              "+ f(slopeid, t, copy = 'numid')"))
  
  is_binary = all(data[[node]] %in% c(0, 1), na.rm = TRUE)
  if (is_binary) {
    
    model = tryCatch(
      inla(formula, family = "binomial", data = data,
           control.inla = list(int.strategy = "eb")),
      error = function(e) {
        message("Error fitting model for node ", node, ": ", e$message)
        return(NULL)
      }
    )
  } else {  
    data$t = data$t +1
    model = tryCatch(
      inla(formula, family = "gaussian", data = data,
           control.inla = list(int.strategy = "eb")),
      error = function(e) {
        message("Error fitting model for node ", node, ": ", e$message)
        return(NULL)
      }
    )
  }
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

#=== Custom score function INLA user defined Wishart prior ====#

fitME2 = function(node, parents, group, t, data) {
  data$Y1 <- as.numeric(as.character(data$Y1))
  data$Y2 <- as.numeric(as.character(data$Y2))
  data$Y3 <- as.numeric(as.character(data$Y3))
  data$Y4 <- as.numeric(as.character(data$Y4))
  nid <- length(unique(data[[group]]))  
  data$numid <- as.numeric(as.factor(data[[group]]))  
  data$slopeid <- data$numid +nid
  
  rhs = paste(c("1", setdiff(parents, t), t), collapse = " + ")
  
  formula = as.formula(paste0(node, "~", rhs, 
                              "+ f(numid, model = 'iidkd', n = ", 2 * nid, 
                              ", order = 2, constr = FALSE, hyper = list(theta1 = list(param = c(6,1,1,0))))", 
                              "+ f(slopeid, t, copy = 'numid')"))
  
  is_binary = all(data[[node]] %in% c(0, 1), na.rm = TRUE)
  if (is_binary) {
    model = tryCatch(
      inla(formula, family = "binomial", data = data,
           control.inla = list(int.strategy = "eb")),
      error = function(e) {
        message("Error fitting model for node ", node, ": ", e$message)
        return(NULL)
      }
    )
  } else {  
    data$t = data$t +1
    model = tryCatch(
      inla(formula, family = "gaussian", data = data,
           control.inla = list(int.strategy = "eb")),
      error = function(e) {
        message("Error fitting model for node ", node, ": ", e$message)
        return(NULL)
      }
    )
  }
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
#====Custom score function INLA default LKJ prior===#
fitME3 <- function(node, parents, group, t, data) {
  
  data$Y1 <- as.numeric(as.character(data$Y1))
  data$Y2 <- as.numeric(as.character(data$Y2))
  data$Y3 <- as.numeric(as.character(data$Y3))
  data$Y4 <- as.numeric(as.character(data$Y4))
  
  data$ID <- as.numeric(data$ID)
  data$rintercept <- rep(1, each = nrepeat * nindividual)
  data$rSlope <- rep(2, each = nrepeat * nindividual)
  n <- 2
  eta_values <- c(10, 30)  # Try these eta values in order
  cfam <- list(hyper = list())
  cinla <- list(int.strategy = 'eb')
  rhs <- paste(c("1", setdiff(parents, t), t), collapse = " + ")

  for (eta in eta_values) {
    #message("Trying eta = ", eta, " for node ", node)
    cmodel <- cgeneric(model = "LKJ", n = n, eta = eta)
    
    formula <- as.formula(paste0(node, " ~ ", rhs,
                                 " + f(rintercept, model = cmodel, replicate = ID)",
                                 " + f(rSlope, ", t, ", copy = 'rintercept', replicate = ID)"))
    is_binary = all(data[[node]] %in% c(0, 1), na.rm = TRUE)
    
    if (is_binary){
      model <- try(
        inla(formula, family = "binomial", data = data,
             control.inla = cinla,
             control.family = cfam),
        silent = TRUE
      )
      
      # If the model is not an error object, return it
      if (!inherits(model, "try-error")) {
        return(model)
      } else {
        #message("INLA crashed with eta = ", eta, ". Trying next eta...")
      }
    } else{
    
    data$t = data$t +1
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
#============================================#
#====Custom score function INLA user defined LKJ prior===#
fitME4 <- function(node, parents, group, t, data) {
  
  data$Y1 <- as.numeric(as.character(data$Y1))
  data$Y2 <- as.numeric(as.character(data$Y2))
  data$Y3 <- as.numeric(as.character(data$Y3))
  data$Y4 <- as.numeric(as.character(data$Y4))
  
  data$ID <- as.numeric(data$ID)
  data$rintercept <- rep(1, each = nrepeat * nindividual)
  data$rSlope <- rep(2, each = nrepeat * nindividual)
  
  n <- 2
  eta_values <- c(13, 30)  # Try these eta values in order
  
  cinla <- list(int.strategy = 'eb')
  cfam <- list(hyper = list())
  rhs <- paste(c("1", setdiff(parents, t), t), collapse = " + ")
  
  for (eta in eta_values) {
    #message("Trying eta = ", eta, " for node ", node)
    cmodel <- cgeneric(model = "LKJ", n = n, eta = eta)
    
    formula <- as.formula(paste0(node, " ~ ", rhs,
                                 " + f(rintercept, model = cmodel, replicate = ID)",
                                 " + f(rSlope, ", t, ", copy = 'rintercept', replicate = ID)"))
    is_binary = all(data[[node]] %in% c(0, 1), na.rm = TRUE)
    
    if (is_binary){
      model <- try(
        inla(formula, family = "binomial", data = data,
             control.inla = cinla,
             control.family = cfam),
        silent = TRUE
      )
      
      # If the model is not an error object, return it
      if (!inherits(model, "try-error")) {
        return(model)
      } else {
        #message("INLA crashed with eta = ", eta, ". Trying next eta...")
      }
    } else{
      
      data$t = data$t +1
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
#==== F1-score for the structure learning ===#
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

#== Simulation study ==#
# Random DAG used to simulate the data
TMBN = model2network("[ID][t][Y1|ID:t][Y2|ID:t][Y8|ID:t][Y3|Y1:ID:t][Y4|Y1:ID:t][Y5|Y2:ID:t][Y6|Y1:Y5:ID:t][Y9|Y1:Y4:Y8:ID:t][Y7|Y5:Y6:ID:t][Y10|Y4:Y6:t:ID]")

shd_lme = list() 
shd_inla = list()
shd_inla1 = list()
shd_inla2 = list()
shd_inla3 = list()

F1_Slme <- list()
F1_Sinla<-list()
F1_Sinla1<- list()
F1_Sinla2<- list()
F1_Sinla3<- list()

F1_Slme_1 <- list()
F1_Sinla_1<-list()
F1_Sinla1_1<- list()
F1_Sinla2_1<- list()
F1_Sinla3_1<- list()

for(k in 1:2) {
nindividual = 50
nrepeat = 5
N = nrepeat*nindividual

data = data.frame(
  "ID" = rep(1:nindividual, each = nrepeat),
  "t" = rep(1:nrepeat, nindividual)
)
data$t <- data$t - 1
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

data$b0 = 0.5
data$b1 = -0.3

data$logodds = data$b0 + data$b1*data$t + data$u1*data$t + data$u0

data$probability = logistic(data$logodds)

data$Y1 = rbinom(n = 1:nrow(data),
                 size = 1,
                 prob = data$probability)

data <- data[, !(names(data) %in% c("u0","u1", "b0","b1","probability", "logodds"))]


#== Y2|ID:t
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

data$b0 = 0.7
data$b1 = -0.3

data$logodds = data$b0 + data$b1*data$t + data$u1*data$t + data$u0

data$probability = logistic(data$logodds)

data$Y2 = rbinom(n = 1:nrow(data),
                 size = 1,
                 prob = data$probability)

data <- data[, !(names(data) %in% c("u0","u1", "b0","b1", "probability", "logodds"))]

#=== Y3|Y1:ID:t
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

data$b0 = -0.5
data$b1 = 0.3
data$b2 = -0.9

data$logodds = data$b0 + data$b1*data$t +
  data$b2*data$Y1 + data$u1*data$t + 
  data$u0
data$probability = logistic(data$logodds)

data$Y3 = rbinom(n = 1:nrow(data),
                 size = 1,
                 prob = data$probability)


data <- data[, !(names(data) %in% c("u0","u1", "b0","b1", "b2", "probability", "logodds"))]

#=== Y4|Y1:ID:t
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

data$b0 = 0.32
data$b1 = -0.56
data$b2 = 0.48 

data$logodds =  data$b0 + data$b1*data$t + data$b2*data$Y1 + 
  data$u1*data$t + data$u0 

data$probability = logistic(data$logodds)

data$Y4 = rbinom(n = 1:nrow(data),
                 size = 1,
                 prob = data$probability)

data <- data[, !(names(data) %in% c("u0","u1", "b0","b1", "b2","logodds","probability"))]

#===Y5|Y2:ID:t
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

sigmae1 = sqrt(0.55) #0.95
data$e = rnorm(N, 0, sigmae1)

data$b0 = 2.63
data$b1 = 1.19
data$b2 = -0.5

data$Y5 = data$b0 + data$b1*data$t + data$b2*data$Y2 + 
  data$u1*data$t + data$u0 + data$e

data <- data[, !(names(data) %in% c("u0","u1", "e", "b0","b1", "b2"))]


#===Y6|Y1,Y5:t:ID
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
  data$b3*data$Y5 + data$u1*data$t + 
  data$u0 + data$e

data <- data[, !(names(data) %in% c("u0","u1","e", "b0","b1", "b2","b3"))]


#===Y7|Y5,Y6:ID:t
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

sigmae1 = sqrt(0.57) # 87
data$e = rnorm(N, 0, sigmae1)

data$b0 = 1.28
data$b1 = -0.43
data$b2 = 0.53
data$b3 =  -0.8

data$Y7 = data$b0 + data$b1*data$t + data$b2*data$Y5 + 
  data$b3*data$Y6 + data$u1*data$t + data$u0 + data$e


data <- data[, !(names(data) %in% c("u0","u1", "e", "b0","b1", "b2","b3"))]


#== Y8|ID:t
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

sigmae1 = sqrt(0.52) # 0.65
data$e = rnorm(N, 0, sigmae1)

data$b0 = 1
data$b1 = 0.5
data$Y8 = data$b0 + data$b1*data$t + data$u1*data$t + data$u0 + data$e

data <- data[, !(names(data) %in% c("u0","u1", "e", "b0","b1"))]

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

data$Y9 = data$b0 +data$b1*data$t + data$b2*data$Y1 +data$b3*data$Y4 +data$b4*data$Y8 + 
  data$u1*data$t + data$u0 + data$e

data <- data[, !(names(data) %in% c("u0","u1", "e", "b0","b1", "b2","b3","b4"))]


#===Y10|Y4,Y6:ID:time
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
data$b3 = 1.82

data$Y10 = data$b0 + data$b1*data$t + data$b2*data$Y4 +data$b3*data$Y6 + 
  data$u1*data$t + data$u0 + data$e

data <- data[, !(names(data) %in% c("u0","u1", "e", "b0","b1", "b2","b3"))]

#== True MBN ==#
# True DAG
TMBN = model2network("[ID][t][Y1|ID:t][Y2|ID:t][Y8|ID:t][Y3|Y1:ID:t][Y4|Y1:ID:t][Y5|Y2:ID:t][Y6|Y1:Y5:ID:t][Y9|Y1:Y4:Y8:ID:t][Y7|Y5:Y6:ID:t][Y10|Y4:Y6:t:ID]")


#=== Structure learning ===#
data$ID<- as.factor(data$ID)
data$t<- as.numeric(data$t)
data$Y1<- as.factor(data$Y1)
data$Y2<- as.factor(data$Y2)
data$Y3<- as.factor(data$Y3)
data$Y4<- as.factor(data$Y4)

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
dag.lme_cp <- cpdag(dag.lme)
dag.inla_cp <-cpdag(dag.inla)
dag.inla1_cp<-cpdag(dag.inla1)
dag.inla2_cp<-cpdag(dag.inla2)
dag.inla3_cp<-cpdag(dag.inla3)

shd_lme[[k]] = shd(TMBN_cp, dag.lme_cp) 
shd_inla[[k]] = shd(TMBN_cp, dag.inla_cp)
shd_inla1[[k]] = shd(TMBN_cp, dag.inla1_cp) 
shd_inla2[[k]] = shd(TMBN_cp, dag.inla2_cp) 
shd_inla3[[k]] = shd(TMBN_cp, dag.inla3_cp) 


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

F1_Slme <- unlist(F1_Slme)
F1_Sinla <-unlist(F1_Sinla)
F1_Sinla1<- unlist(F1_Sinla1)
F1_Sinla2<- unlist(F1_Sinla2)
F1_Sinla3<- unlist(F1_Sinla3)

F1_Slme_1 <- unlist(F1_Slme_1)
F1_Sinla_1 <-unlist(F1_Sinla_1)
F1_Sinla1_1<- unlist(F1_Sinla1_1)
F1_Sinla2_1<- unlist(F1_Sinla2_1)
F1_Sinla3_1<- unlist(F1_Sinla3_1)


stats<-data.frame(shd_lme, shd_inla, shd_inla1, shd_inla2, shd_inla3, 
                  F1_Slme, F1_Sinla, F1_Sinla1, F1_Sinla2, F1_Sinla2,
                  F1_Slme_1, F1_Sinla_1, F1_Sinla1_1, F1_Sinla2_1, 
                  F1_Sinla3_1)

out.file3 = "D:/King_Abdullah/MBN_INLA2/output_file3.rds"
saveRDS(stats, file = out.file3)
#readRDS("D:/King_Abdullah/MBN_INLA2/output_file3.rds")

