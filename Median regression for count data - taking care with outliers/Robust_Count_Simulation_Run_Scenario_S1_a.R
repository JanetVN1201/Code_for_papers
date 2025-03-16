source("Robust_Count_Simulation_Master.R")

demoRes <- SaveModels(
  scenarioName = "S1",
  nsim = 10,
  N = 100,
  cVal = 1,
  beta0True = 2.0,
  beta1True = 0.3,
  alphaTrue = 0.6,
  etaTrue = 5,
  deltaTrue = 0.45,
  burnin = 4000,
  sample = 4000,
  thin = 2,
  adapt = 1000,
  nChains = 2
)