##NEW
# Running simulations
source("Quantile_joint/Sim_Aug_24/Sim.R")
NSim <- 200
SimParams <- list(
  nsujet = 500,   # Number of individuals
  psi = 0.8,      # Kumaraswamy scale parameter
  SD1 = 0.2,      # SD of random intercept
  SD2 = 0.2,      # SD of random slope
  rho = 0.1,      # Random intercept/slope correlation coefficient
  gapLongi = 0.1, # Gap between longitudinal measurements
  gap = 0.01,     # Used to generate many time points for permutation algorithm
  followup = 5    # Follow-up time
)

q0 = 0.5
cens0 = 0.006

GenerateParamSets <- function(q = q0) {
  # Varying parameters
  q_val <- q
  beta_0_val <- c(-0.5)
  beta_1_val <- c(-0.4, 0.4)
  beta_2_val <- c(0.4)
  beta_3_val <- c(0.5)
  gamma_1_val <- c(0.1)
  gamma_2_val <- c(0.3)
  phi_val <- c(-0.1, 0.1)
  
  ParamSets <- expand.grid(q = q_val,             # Quantile
                           beta_0 = beta_0_val,   # Longitudinal intercept
                           beta_1 = beta_1_val,   # Longitudinal slope
                           beta_2 = beta_2_val,   # Longitudinal binary covariate
                           beta_3 = beta_3_val,   # Longitudinal continuous covariate
                           gamma_1 = gamma_1_val, # Survival binary covariate [IGNORED FOR SIMPLICITY]
                           gamma_2 = gamma_2_val, # Survival continuous covariate [IGNORED FOR SIMPLICITY]
                           phi = phi_val)         # Association parameter
  return(ParamSets)
}

ParamSets <- GenerateParamSets()

inla.setOption(inla.timeout = 400)

AllSimRes <- RunSimulations(ParamSets, NSim, cens0, q0=q0)

# Combining all simulation results
CombRes <- do.call(rbind, AllSimRes)

# Sort combined results
CombRes <- CombRes %>%
  arrange(Set, Sort)

# Reordering columns
nonStaticParams <- IdentifyNonStaticParams(ParamSets)
CombRes <- CombRes[c("Set", "q", nonStaticParams, setdiff(names(CombRes), c("Set", "q", nonStaticParams)))]
CombRes$Q <- CombRes$q

# Identify non-static parameters
ParamSets <- GenerateParamSets()
nonStaticParams <- IdentifyNonStaticParams(ParamSets)

# Display non-static parameters for first row only
lastSet <- ""
for (i in 1:nrow(CombRes)) {
  if (CombRes$Set[i] != lastSet) {
    lastSet <- CombRes$Set[i]
  } else {
    for (param in nonStaticParams) {
      CombRes[[param]][i] <- ""
    }
  }
}

# Drop unwanted 'q.1' column
if("q.1" %in% names(CombRes)) {
  CombRes <- CombRes[ , !(names(CombRes) %in% c("q.1"))]
}

print(CombRes)

n <- 500000
eventRandom <- round(rexp(n, cens0) + 1, 0)
censorRandom <- round(runif(n, 1, 501), 0)
isCensored <- ifelse(eventRandom > censorRandom, 1, 0)
CRate <- round(sum(isCensored)/n,2)
print(paste("Censoring rate:", CRate))

# Generate LaTeX tables
UniqQs <- unique(CombRes$Q)
library(xtable)

for(QLevel in UniqQs) {
  SubDat <- subset(CombRes, Q == QLevel, select = -c(Set, Sort, RMSE, q, Q))
  
  SubDat$Value <- sprintf("%.2f", SubDat$Value)
  SubDat$Bias <- sprintf("%.4f", SubDat$Bias)
  SubDat$Coverage <- sprintf("%.3f", SubDat$Coverage)
  
  LaTeXTable <- xtable(SubDat, floating = FALSE, include.rownames = FALSE, include.colnames = FALSE, sanitize.text.function = identity)
  FileName <- paste0("Quantile_joint/Sim_Aug_24/ResultsLaTeX_q", gsub(" ", "_", QLevel), "_CR", gsub(" ", "_", CRate), "_N", gsub(" ", "_", SimParams$nsujet), ".txt")
  print(LaTeXTable, 
        file = FileName, 
        floating = FALSE, 
        include.rownames = FALSE, 
        include.colnames = FALSE,
        sanitize.text.function = identity,
        hline.after = c(-1, 0, nrow(SubDat)),
        only.contents = TRUE)
  write.csv(CombRes, paste0("Quantile_joint/Sim_Aug_24/ResultsLaTeX_q", gsub(" ", "_", QLevel), "_CR", gsub(" ", "_", CRate), "_N", gsub(" ", "_", SimParams$nsujet), ".csv"))
  cat("LaTeX file generated for Q =", QLevel, "at", FileName, "\n")
}
