library(dplyr)
library(purrr)
library(stringr)

CombineScenarios <- function() {
  rdsFiles <- list.files(pattern = "^Scenario_.*\\.rds$")
  if(length(rdsFiles) == 0) {
    message("No scenario RDS files found.")
    return(NULL)
  }
  allObjs <- lapply(rdsFiles, readRDS)
  bigDF <- do.call(rbind, lapply(allObjs, function(x) x$resultsLong))
  bigDF
}

SummarizeSimRes <- function(allDataLong) {
  allDataLong2 <- allDataLong %>%
    mutate(true_value = dplyr::case_when(
      param == "beta0" ~ beta0True,
      param == "beta1" ~ beta1True,
      param == "alpha" ~ alphaTrue,
      param == "eta" ~ etaTrue,
      param == "delta" ~ deltaTrue,
      TRUE ~ NA_real_
    ))
  summaryTable <- allDataLong2 %>%
    group_by(
      scenario, N, cVal, beta0True, beta1True, alphaTrue, etaTrue, deltaTrue,
      model, param, true_value
    ) %>%
    summarize(
      n_reps = n(),
      Bias = mean(estimate - true_value, na.rm = TRUE),
      RMSE = sqrt(mean((estimate - true_value)^2, na.rm = TRUE)),
      Coverage = mean(cred_lower <= true_value & cred_upper >= true_value, 
                      na.rm = TRUE),
      AvgCIlen = mean(cred_upper - cred_lower, na.rm = TRUE),
      PowerExcl0 = mean((cred_upper < 0) | (cred_lower > 0), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(scenario, N, cVal, model, param)
  summaryTable
}

allData <- CombineScenarios()
summaryDF <- SummarizeSimRes(allData)
print(summaryDF)
View(summaryDF)