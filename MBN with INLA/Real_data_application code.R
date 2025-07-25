library(haven)
library(openxlsx)
library(ggplot2)
library(glmmTMB)
library(dplyr)
library(readxl)
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
library(tidyr)
library(ggplot2)

data<- read_excel("path/Sub_data.xlsx")
names(data)[names(data) == "gasterointertile"] <- "AGE"
names(data)[names(data) == "No_RTUF"] <- "RUTF"
names(data)[names(data) == "Diarrhea_cat"] <- "Diarrhea"
names(data)[names(data) == "whz06"] <- "WHZ"
names(data)[names(data) == "pnemonia"] <- "pneumonia"

#=== structure learning used for the analysis =========#
#===== MLE
fitME = function(node, parents, group, t, data) {
  parents_n = setdiff(parents, c(t, "RUTF"))
  rhs = paste(c("1", parents_n, paste0(t, "*", "RUTF")), collapse = " + ")
  
  formula = as.formula(paste(node, "~", rhs, "+ (1 +", t, "|", group, ")"))
  
  is_binary = all(data[[node]] %in% c(0, 1), na.rm = TRUE)
  
  if (is_binary) {
    model = try(glmer(formula, data = data, family = binomial, 
                      control = glmerControl(optimizer = "bobyqa")))
  } else {
    data$time = data$time+1
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
  if (node == "RUTF")
    return(0)
  if (node == "Gender")
    return(0)
  if (node == "Age")
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

#=== INLA default wishart 
fitME1 = function(node, parents, group, t, data) {
  
  data$Diarrhea <- as.numeric(as.character(data$Diarrhea))
  data$AGE <- as.numeric(as.character(data$AGE))
  data$pneumonia <- as.numeric(as.character(data$pneumonia))
  data$Dehayderation <- as.numeric(as.character(data$Dehayderation))
  
  nid <- length(unique(data[[group]]))  
  data$numid <- as.numeric(as.factor(data[[group]]))  
  data$slopeid <- data$numid +nid
  
  parents_n = setdiff(parents, c(t, "RUTF"))
  rhs = paste(c("1", parents_n, paste0(t, "*", "RUTF")), collapse = " + ")
  rhs = paste(c("1", parents_n, t, "RUTF"), collapse = " + ")
  
  
  formula = as.formula(paste0(node, "~", rhs, 
                              "+ f(numid, model = 'iidkd', n = ", 2 * nid, 
                              ", order = 2, constr = FALSE, hyper = list(theta1 = list(param = c(100,1,1,0))))", 
                              "+ f(slopeid, time, copy = 'numid')"))
  
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
    data$time = data$time +1
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
  if (node == "RUTF")
    return(0)
  if (node == "Gender")
    return(0)
  if (node == "Age")
    return(0)
  
  require(INLA)
  
  local.distribution =
    fitME1(node = node, parents = setdiff(parents, c(args$group, args$t)),  # Ensure both are removed
           group = args$group, t = args$t, data = data)
  
  if (is.null(local.distribution))
    return(-Inf)
  
  return(local.distribution$mlik[1])
  
}


#==== INLA wishart
fitME2 = function(node, parents, group, t, data) {
  
  data$Diarrhea <- as.numeric(as.character(data$Diarrhea))
  data$AGE <- as.numeric(as.character(data$AGE))
  data$pneumonia <- as.numeric(as.character(data$pneumonia))
  data$Dehayderation <- as.numeric(as.character(data$Dehayderation))
  
  nid <- length(unique(data[[group]]))  
  data$numid <- as.numeric(as.factor(data[[group]]))  
  data$slopeid <- data$numid +nid
  
  parents_n = setdiff(parents, c(t, "RUTF"))
  rhs = paste(c("1", parents_n, paste0(t, "*", "RUTF")), collapse = " + ")
  rhs = paste(c("1", parents_n, t, "RUTF"), collapse = " + ")
  
  
  formula = as.formula(paste0(node, "~", rhs, 
                              "+ f(numid, model = 'iidkd', n = ", 2 * nid, 
                              ", order = 2, constr = FALSE, hyper = list(theta1 = list(param = c(6,1,1,0))))", 
                              "+ f(slopeid, time, copy = 'numid')"))
  
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
    data$time = data$time +1
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
  if (node == "RUTF")
    return(0)
  if (node == "Gender")
    return(0)
  if (node == "Age")
    return(0)
  
  require(INLA)
  
  local.distribution =
    fitME2(node = node, parents = setdiff(parents, c(args$group, args$t)),  # Ensure both are removed
           group = args$group, t = args$t, data = data)
  
  if (is.null(local.distribution))
    return(-Inf)
  
  return(local.distribution$mlik[1])
  
}

#=== INLA LKJ
fitME3 <- function(node, parents, group, t, data) {
  
  data$Diarrhea <- as.numeric(as.character(data$Diarrhea))
  data$AGE <- as.numeric(as.character(data$AGE))
  data$pneumonia <- as.numeric(as.character(data$pneumonia))
  data$Dehayderation <- as.numeric(as.character(data$Dehayderation))
  
  data$ID <- as.numeric(data$ID)
  data$rintercept <- rep(1, each = 981)
  data$rSlope <- rep(2, each = 981)
  n <- 2
  eta_values <- c(10, 30)  # the value 30 is used for debugging purpose
  
  cinla <- list(int.strategy = 'eb')
  cfam <- list(hyper = list())
  
  parents_n = setdiff(parents, c(t, "RUTF"))
  rhs = paste(c("1", parents_n, paste0(t, "*", "RUTF")), collapse = " + ")
  
  for (eta in eta_values) {
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
      
      if (!inherits(model, "try-error")) {
        return(model)
      } else {
    
      }
    } else{
      
      data$time = data$time +1
      model <- try(
        inla(formula, family = "gaussian", data = data,
             control.mode = list(
               theta = c(4, 0),
               restart = TRUE),
             control.inla = cinla),
        silent = TRUE
      )
      
      
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

#=== structure learning

data$ID<- as.factor(data$ID)
data$time<- as.numeric(data$time)
data$RUTF<- as.numeric(data$RUTF)
data$Age<- as.numeric(data$Age)
data$WHZ<- as.numeric(data$WHZ)


data$Gender<- as.factor(data$Gender)
data$Diarrhea<- as.factor(data$Diarrhea)
data$Dehayderation<- as.factor(data$Dehayderation)
data$AGE<- as.factor(data$AGE)
data$pneumonia<- as.factor(data$pneumonia)


source_nodes = c("ID", "time", "RUTF")  
no_incoming = c("ID", "time", "RUTF", "Age", "Gender")  

all_nodes = names(data) 

whitelist.arcs <- expand.grid(
  from = source_nodes,
  to = setdiff(all_nodes, c(source_nodes, "Age", "Gender"))
)

blacklist.incoming <- expand.grid(
  from = all_nodes,
  to = no_incoming
)

dag.lme <- hc(data, score = "custom-score", fun = scoreME, 
              args = list(group = "ID", t= "time"), 
              whitelist = whitelist.arcs, 
              blacklist = blacklist.incoming)

dag.inla <- hc(data, score = "custom-score", fun = scoreME1, 
               args = list(group = "ID", t= "time"), 
               whitelist = whitelist.arcs, 
               blacklist = blacklist.incoming)

dag.inla1 <- hc(data, score = "custom-score", fun = scoreME2, 
                args = list(group = "ID", t= "time"), 
                whitelist = whitelist.arcs, 
                blacklist = blacklist.incoming)

dag.inla2 <- hc(data, score = "custom-score", fun = scoreME3, 
                args = list(group = "ID", t= "time"), 
                whitelist = whitelist.arcs, 
                blacklist = blacklist.incoming)
dag.inla3 <- hc(data, score = "custom-score", fun = scoreME4, 
                args = list(group = "ID", t= "time"), 
                whitelist = whitelist.arcs, 
                blacklist = blacklist.incoming)


highlighted_edges <- matrix(c(
  "ID", "Diarrhea",
  "ID", "Dehayderation",
  "ID", "WHZ",
  "ID", "pneumonia",
  "ID", "AGE"
), byrow = TRUE, ncol = 2,
dimnames = list(NULL, c("from", "to")))

par(mfrow = c(2, 2))
dag_lme<- graphviz.plot(dag.lme, 
                        highlight = list(arcs = highlighted_edges, col = "green", lwd =2),
                        layout = "dot", 
                        shape = "ellipse", 
                        fontsize = 14)
mtext("a) MBN using MLE", side = 1, line = 4.5, adj = 0, cex = 1)

dag_inla<- graphviz.plot(dag.inla, 
                         highlight = list(arcs = highlighted_edges, col = "green", lwd =2),
                         layout = "dot", 
                         shape = "ellipse", 
                         fontsize = 14)
mtext("b) MBN using Bayes 1", side = 1, line = 4.5, adj = 0, cex = 1)

dag_inla1<- graphviz.plot(dag.inla1, 
                          highlight = list(arcs = highlighted_edges, col = "green", lwd =2),
                          layout = "dot", 
                          shape = "ellipse", 
                          fontsize = 12)
mtext("c) MBN using Bayes 2", side = 1, line = 4, adj = 0, cex = 1)

dag_inla2<- graphviz.plot(dag.inla2, 
                          highlight = list(arcs = highlighted_edges, col = "green", lwd =2),
                          layout = "dot", 
                          shape = "ellipse", 
                          fontsize = 12)
mtext("d) MBN using Bayes 3", side = 1, line = 4, adj = 0, cex = 1)

