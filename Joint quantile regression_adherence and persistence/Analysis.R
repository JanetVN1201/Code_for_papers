###########
###SETUP###
###########

# Load packages
library(tidyverse)
library(ggplot2)
library(quantreg)
library(reshape2)
library(dplyr)
library(lme4)
library(glmmTMB)
library(parallel)
library(survival)
library(survminer)
library(INLA)
library(multinma)
library(sn)
library(zoo)
library(splines)
library(purrr)
library(INLAjoint)
library(xtable)
library(tidyr)
library(xtable)

set.seed(123)

# Working directory and read in datasets
setwd('C:/Users/a239866/OneDrive - Syneos Health/Divan @ Syneos Health - Linked Files/Research/Drug Adherence')
Cores <- ceiling(detectCores()/3) # Number of processing cores

TimePoints <- seq(0, 300, length.out = 100) # Time points for survival plots
TrtColors <- c("No Intervention" = "red", "Intervention" = "blue") # Map the treatments to colors

###############
###FUNCTIONS###
###############

RowNamesColumn <- function(df, NewColName = "ParameterName") {
  df[NewColName] <- rownames(df)
  return(df)
}

kable <- function(...) {
  knitr::kable(..., format = "html", table.attr = "style = \"color: white;\"")
}
shortestinterval <- function(postsims, width=0.95) { # Coded by Sean van der Merwe, UFS
  postsims |> sort() -> sorted.postsims
  round(length(postsims)*width) -> gap
  sorted.postsims |> diff(gap) |> which.min() -> pos
  sorted.postsims[c(pos, pos + gap)] }

theme_set(theme_bw())

# Function to simulate from qkumar distribution
rkum <- function(n, Kappa, psi, q) {
  beta <- log(1 - q)/log(1 - exp(-psi)) # Solve beta of qkumar
  alpha <- log(1 - (1 - q)^(1/beta))/log(Kappa) # Solve alpha of qkumar
  (1 - (1 - runif(n))^(1/beta))^(1/alpha) # Transform using inverse CDF
}

# Function to impute NA values with the mean for a specific column
ImputeNAMean <- function(data, column_name) {
  if (!column_name %in% names(data)) {
    stop("Column not found in the dataframe.")
  }
  column_mean <- mean(data[[column_name]], na.rm = TRUE)
  data[[column_name]][is.na(data[[column_name]])] <- column_mean
  return(data)
}

##########################################################
#########DATA MANAGEMENT AND PRELIMINARY ANALYSIS#########
##########################################################

# Read the datasets
Compliance <- read.csv('Datasets/AllCompliance.csv')
Persist <- read.csv('Datasets/persistence.csv')
ITT <- read.csv('Datasets/itt.csv')
DMData <- read.csv('Datasets/BaseChar.csv')
BaseChar <- read.csv('Datasets/gp.csv')
PatInfo <- read.csv('Datasets/patinfo.csv')

BaseChar <- ImputeNAMean(BaseChar, 'CARDRISK')
BaseChar <- ImputeNAMean(BaseChar, 'HDLCHOL')

# Match the ID fields
ITT$id <- ITT$RANDNR

# Data processing and merging
Compliance <- Compliance %>%
  filter(phase %in% 0:30) %>%
  mutate(
    phasen = as.numeric(phase),
  ) %>%
  arrange(id, phasen) %>%
  group_by(id) %>%
  mutate(
    CumulativeDur = cumsum(durdint),
    StudyDay = CumulativeDur
  ) %>%
  ungroup() %>%
  left_join(Persist, by = "id") %>%
  full_join(ITT, by = "id") %>%
  filter(itt == '1')
Compliance$ID = as.integer(as.factor(Compliance$id)) #Renumber subject IDs chronologically
Compliance <- Compliance %>% # Dataset in long format (interval-censored)
  arrange(ID, StudyDay) %>%
  group_by(ID) %>%
  mutate(
    time1 = lag(StudyDay, order_by = StudyDay, default = 0),
    time2 = StudyDay,
    is_LastRecord = row_number() == n()
  ) %>%
  mutate(EventTime = if_else(survstat == 1, pers, NA_real_)) %>%
  fill(EventTime, .direction = "downup") %>%
  mutate(
    time2 = case_when(
      !is.na(EventTime) & time2 >= EventTime ~ EventTime,
      is_LastRecord & time2 > pers ~ pers,
      TRUE ~ time2
    )
  ) %>%
  filter(is.na(EventTime) | time1 <= EventTime) %>%
  mutate(
    EventPers = if_else(!is.na(EventTime) & time2 == EventTime, 1, 0),
    EventPers = if_else(is_LastRecord & !is.na(EventTime) & time2 < EventTime, 1, EventPers)
  ) %>%
  ungroup()
Compliance$StudyDay <- Compliance$time2
Compliance$time2 <- pmin(Compliance$time2, Compliance$pers)
Compliance <- merge(Compliance, BaseChar, by = "RANDNR")

Compliance <- merge(Compliance, DMData[, c("RANDNR", "Age")], by = "RANDNR")
Compliance <- merge(Compliance, PatInfo[, c("RANDNR", "Group")], by = "RANDNR")
Compliance$int <- ifelse(is.na(Compliance$int), 
                         ifelse(Compliance$Group == 2, 0, Compliance$Group), 
                         Compliance$int)

Compliance <- Compliance %>% # Impute variables where MEMS data are not available
  mutate(
    phasen = ifelse(base_imputed == 1, 0, phasen),
    cod = ifelse(base_imputed == 1, 0.01, cod),
    StudyDay = ifelse(base_imputed == 1, 1, StudyDay),
    pers = ifelse(base_imputed == 1, 1, pers),
    EventPers = ifelse(base_imputed == 1, 1, EventPers),
    survstat = ifelse(base_imputed == 1, 1, survstat),
    time2 = ifelse(base_imputed == 1, 1, time2),
    EventTime = ifelse(base_imputed == 1, 1, EventTime),
  )

Compliance <- Compliance %>%
  mutate(
    StudyMonth = StudyDay/30,
    nummed = cod,
    cod = ifelse(cod > 0.99, 0.99, ifelse(cod < 0.01, 0.01, cod)),
    LogitCOD = qlogis(cod),
    MonthLabel = ifelse(phasen == 0, "Baseline", as.character(phasen)),
    IntLabel = recode(factor(int), '0' = 'No Intervention', '1' = 'Intervention'),
    Treat = ifelse(int == 0, 2, 1)
  )

Compliance <- Compliance %>%
  group_by(ID) %>%
  mutate(Q50PerPat = quantile(cod, probs = 0.5, na.rm = TRUE)) %>%
  ungroup()
OverallQ50 <- quantile(Compliance$cod, probs = 0.5, na.rm = TRUE)
Compliance <- Compliance %>%
  mutate(BelowQ50 = ifelse(Q50PerPat < OverallQ50, 0, 1),
         OverallQ50 = OverallQ50,
         BelowQ50Label = recode(factor(BelowQ50), 
                                '1' = 'Above Median', 
                                '0' = 'Below Median'))

write.csv(Compliance, "Manuscript/Output/AdhereData.csv")

Persistence <- Compliance %>%
  group_by(ID) %>%
  slice_min(order_by = as.numeric(phasen)) %>%
  ungroup()
Persistence <- Persistence %>%
  mutate(IntLabelNum = ifelse(IntLabel == "Intervention", 1, 0))
Persistence$PersMonth <- Persistence$pers/30

# Sort MonthLabel
NumericLevels <- sort(unique(as.numeric(Compliance$MonthLabel[Compliance$MonthLabel != "Baseline"])))
SortedLevels <- c("Baseline", as.character(NumericLevels))
Compliance$MonthLabel <- factor(Compliance$MonthLabel, levels = SortedLevels)

# Calculate the number of observations for each phasen, IntLabel, and survstat
CountData <- Compliance %>% 
  group_by(phasen, IntLabel, survstat) %>% 
  summarise(n = n()) %>%
  mutate(MonthLabel = ifelse(phasen == 0, "Baseline", as.character(phasen)))

# Function to relabel the survstat facet
SurvStatLabelChange <- function(x) {
  ifelse(x == "1", "Not Persistent", "Persistent")
}

# Boxplot of adherence by intervention group
FilteredBoxplotDataOverall <- Compliance %>% 
  semi_join(CountData %>% group_by(phasen, IntLabel) %>% summarise(n = sum(n)) %>% filter(n > 5), by = c("phasen", "IntLabel"))
FilteredJitterDataOverall <- Compliance %>% 
  semi_join(CountData %>% group_by(phasen, IntLabel) %>% summarise(n = sum(n)) %>% filter(n <= 5), by = c("phasen", "IntLabel"))
OverallBoxplot <- ggplot(Compliance, aes(x = MonthLabel, y = cod, fill = IntLabel)) +
  geom_boxplot(data = FilteredBoxplotDataOverall) +
  geom_jitter(data = FilteredJitterDataOverall, width = 0.2) +
  labs(x = "Month", 
       y = "Adherence", 
       fill = "Intervention Group") +
  scale_fill_manual(values = c("No Intervention" = "pink", "Intervention" = "lightblue")) +
  theme_light() +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    legend.position = 'bottom', 
    legend.direction = 'horizontal', 
    legend.justification = 'center',
    legend.text = element_text(size = 12),
    panel.grid.major = element_line(colour = "lightgray", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 12)
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  scale_x_discrete(limits = SortedLevels)
pdf("Manuscript/Output/AdhereBoxplotOverall.pdf", width = 7, height = 6)
print(OverallBoxplot)
dev.off()

# Boxplot of adherence by intervention group and persistence outcome (censor flag)
FilteredBoxplotData <- Compliance %>% 
  semi_join(CountData %>% filter(n > 5), by = c("phasen", "IntLabel", "survstat"))
FilteredJitterData <- Compliance %>% 
  semi_join(CountData %>% filter(n <= 5), by = c("phasen", "IntLabel", "survstat"))
histogram <- ggplot(Compliance, aes(x = MonthLabel, y = cod, fill = IntLabel)) +
  geom_boxplot(data = FilteredBoxplotData) +
  geom_jitter(data = FilteredJitterData, width = 0.2) +
  labs(x = "Month", 
       y = "Adherence", 
       fill = "Intervention Group") + 
  scale_fill_manual(values = c("No Intervention" = "pink", "Intervention" = "lightblue")) +
  theme_light() +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    legend.position = 'bottom', 
    legend.direction = 'horizontal', 
    legend.justification = 'center',
    legend.text = element_text(size = 12),
    panel.grid.major = element_line(colour = "lightgray", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 12)
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  scale_x_discrete(limits = SortedLevels) +
  facet_grid(rows = vars(survstat), labeller = labeller(survstat = SurvStatLabelChange))
pdf("Manuscript/Output/AdhereBoxplotPersist.pdf", width = 7, height = 7)
print(histogram)
dev.off()

# Quantile regression of adherence
taus <- c(0.1, 0.25, 0.5, 0.75, 0.9)
models <- lapply(taus, function(tau) rq(LogitCOD ~ int + StudyDay, data = Compliance, tau = tau))
coefs <- sapply(models, coef)
CoefsDF <- as.data.frame(t(coefs))
CoefsDF$quantile <- taus
CoefsLong <- melt(CoefsDF, id.vars = "quantile")
CoefsLong$variable <- recode(CoefsLong$variable, 
                             '(Intercept)' = "Intercept",
                             'int' = "Intervention Effect",
                             'StudyDay' = "Slope (Days)")
qregplot <-
  ggplot(CoefsLong, aes(x = quantile, y = value, color = variable)) +
  geom_line(linewidth = 0.8, linetype = "dashed") +
  geom_point(shape = 17, size = 3) +
  labs(x = "Quantile",
       y = "Coefficient Estimate") + 
  theme_minimal() +
  facet_wrap(~ variable, scales = "free_y", ncol = 2) +
  theme(legend.position = "none",
        axis.line.x = element_line(color = "black", linewidth = 0.5), 
        axis.line.y = element_line(color = "black", linewidth = 0.5)) +
  scale_x_continuous(breaks = taus, minor_breaks = FALSE) +
  guides(color = FALSE)
pdf("Manuscript/Output/AdhereQregplot.pdf", width = 7, height = 5)
qregplot
dev.off()

# Spaghetti plot of adherence over time by intervention group and censor flag
MovingWindow <- 30
Compliance %>% 
  arrange(StudyDay, IntLabel, survstat) %>% 
  group_by(StudyDay, IntLabel, survstat) %>%
  summarise(OverallMedian = median(cod, na.rm = TRUE),
            Q25 = quantile(cod, 0.25, na.rm = TRUE)) %>%
  group_by(IntLabel, survstat) %>% 
  arrange(StudyDay) %>%
  mutate(RollingMedian = rollapply(OverallMedian, width = MovingWindow, FUN = median, align = 'right', fill = NA),
         RollingQ25 = rollapply(Q25, width = MovingWindow, FUN = median, align = 'right', fill = NA)) %>%
  ungroup() -> ComplianceRolling
ProfilePlot <- ggplot(Compliance, aes(x = StudyDay, y = cod, group = ID)) +
  geom_line(aes(colour = IntLabel, linetype = IntLabel), alpha = 0.5) +
  geom_line(data = ComplianceRolling, aes(x = StudyDay, y = RollingMedian, group = 1), color = "black", linetype = "solid", size = 1.3) +
  geom_line(data = ComplianceRolling, aes(x = StudyDay, y = RollingQ25, group = 1), color = "black", linetype = "twodash", size = 1.3) +
  scale_color_manual(values = c("No Intervention" = "red", "Intervention" = "blue")) +
  scale_linetype_manual(values = c("No Intervention" = "dotted", "Intervention" = "dotted")) +
  scale_x_continuous(breaks = seq(0, max(Compliance$StudyDay, na.rm = TRUE), by = 30)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  facet_grid(rows = vars(survstat), cols = vars(IntLabel), labeller = labeller(survstat = SurvStatLabelChange)) + 
  labs(x = "Day", y = "Adherence") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.line = element_line(color = "black", linewidth = 0.5),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
ggsave("Manuscript/Output/AdhereProfilePlot.pdf", plot = ProfilePlot, width = 13, height = 8)

# Plot adherence data for a random sample of subjects
RandomSubjects <- sample(unique(Compliance$ID), 6)
ComplianceSubset <- subset(Compliance, ID %in% RandomSubjects)
FirstObservation <- ComplianceSubset %>% group_by(ID) %>% slice_head(n = 1)
DummyDay0 <- FirstObservation
DummyDay0$StudyDay <- 0
CombinedDummyAndFirst <- bind_rows(DummyDay0, FirstObservation)
IndividualPlot <- ggplot(ComplianceSubset, aes(x = StudyDay, y = cod, group = ID)) +
  geom_line(color = "blue", size = 0.9, linetype = 'dashed') +
  geom_point(color = "black", size = 3.5, shape = 19) +
  geom_line(data = CombinedDummyAndFirst, aes(x = StudyDay, y = cod, group = ID), 
            color = "blue", linetype = "dotted", linewidth = 0.6) +
  scale_x_continuous(breaks = seq(0, max(ComplianceSubset$StudyDay, na.rm = TRUE), by = 30)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(x = "Day", y = "Adherence") +
  theme_minimal() +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  facet_wrap(~ ID, ncol = 2)
ggsave("Manuscript/Output/AdhereIndividualPlot.pdf", plot = IndividualPlot, width = 8, height = 8)

# Fit random intercept-slope mixed model to adherence data
MixedModelWithSlope <- lmer(LogitCOD ~ StudyDay*IntLabel + (1 + StudyDay | ID), data = Compliance) # Mixed model with random slope term dropped
MixedModelNoSlope <- lmer(LogitCOD ~ StudyDay*IntLabel + (1 | ID), data = Compliance)
ModelComparison <- anova(MixedModelWithSlope, MixedModelNoSlope) # Perform a likelihood ratio test to compare the models
print(ModelComparison)

# Plot the Kaplan-Meier survival curve along with the number at risk table
IntLabelLevels <- levels(Persistence$IntLabel)
MyPalette <- setNames(c("red", "blue"), IntLabelLevels)
KMFit <- survfit(Surv(pers, survstat) ~ IntLabel, data = Persistence)
KMCurve <- data.frame(
  Time = KMFit$time,
  KM_est = KMFit$surv,
  KM_lower = KMFit$lower,
  KM_upper = KMFit$upper,
  Treatment = rep(names(KMFit$strata), KMFit$strata),
  Component = "Kaplan-Meier"
)
save(KMCurve, file = "Manuscript/Output/PersistSurvivalKM.RData")
SurvivalKM <- ggsurvplot(
  KMFit,
  data = Persistence,
  risk.table = TRUE,
  risk.table.title = "Remaining in Study",
  risk.table.title.size = 1,
  risk.table.height = 0.2,
  risk.table.y.text = FALSE,
  risk.table.y.text.col = TRUE,
  risk.table.y.text.size = 1,
  tables.theme = theme_cleantable(),
  table.y.text = FALSE,
  pval = FALSE,
  conf.int = FALSE,
  xlab = "Study Day",
  ylab = "Cumulative Proportion Persistent",
  break.time.by = 30,
  legend = "bottom",
  palette = MyPalette,
  legend.labs = IntLabelLevels,
  censor = FALSE
)
SurvivalKM$plot <- SurvivalKM$plot + theme(legend.title = element_blank())
SurvivalKM$plot <- SurvivalKM$plot + 
  theme(
    legend.text = element_text(size = 15),
    panel.grid.major = element_line(color = "lightgray", linetype = "solid", linewidth = 0.2),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(limits = c(0.7, 1), breaks = seq(0.7, 1, by = 0.1))
FinalPlot <- ggarrange(SurvivalKM$plot, SurvivalKM$table, 
                       ncol = 1, 
                       nrow = 2, 
                       heights = c(4, 0.8))
ggsave("Manuscript/Output/PersistSurvivalKM.pdf", FinalPlot, width = 8, height = 7)

# Plot the Kaplan-Meier survival curve by adherence level maintained during study (below/above median)
BelowQ50Levels <- levels(Persistence$BelowQ50Label)
MyPalette <- setNames(c("red", "blue"), BelowQ50Levels)
KMFit <- survfit(Surv(pers, survstat) ~ BelowQ50Label, data = Persistence)
KMCurve <- data.frame(
  Time = KMFit$time,
  KM_est = KMFit$surv,
  KM_lower = KMFit$lower,
  KM_upper = KMFit$upper,
  Treatment = rep(names(KMFit$strata), KMFit$strata),
  Component = "Kaplan-Meier"
)
save(KMCurve, file = "Manuscript/Output/PersistKMLevelMaintain.RData")
SurvivalKM <- ggsurvplot(
  KMFit,
  data = Persistence,
  risk.table = TRUE,
  risk.table.title = "Remaining in Study",
  risk.table.title.size = 1,
  risk.table.height = 0.2,
  risk.table.y.text = FALSE,
  risk.table.y.text.col = TRUE,
  risk.table.y.text.size = 1,
  tables.theme = theme_cleantable(),
  table.y.text = FALSE,
  pval = FALSE,
  conf.int = FALSE,
  xlab = "Study Day",
  ylab = "Cumulative Proportion Persistent",
  break.time.by = 30,
  legend = "bottom",
  palette = MyPalette,
  legend.labs = BelowQ50Levels,
  legend.title = "Adherence",
  censor = FALSE
)
SurvivalKM$plot <- SurvivalKM$plot + theme(legend.text = element_text(size = 15),
                                           legend.title = element_text(size = 15))
SurvivalKM$plot <- SurvivalKM$plot + 
  theme(
    legend.text = element_text(size = 15),
    panel.grid.major = element_line(color = "lightgray", linetype = "solid", linewidth = 0.2),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2))
FinalPlot <- ggarrange(SurvivalKM$plot, SurvivalKM$table, 
                       ncol = 1, 
                       nrow = 2, 
                       heights = c(4, 0.8))
ggsave("Manuscript/Output/PersistKMLevelMaintain.pdf", FinalPlot, width = 8, height = 7)

# Create a log-minus-log survival plot
ggsurvplot(
  KMFit, 
  data = Persistence, 
  fun = "cloglog", 
  risk.table = TRUE, 
  ggtheme = theme_minimal()
)

# Fit a Cox proportional hazards model
CoxModel <- coxph(Surv(pers, survstat) ~ IntLabel + Age + CARDRISK + HDLCHOL, data = Persistence)
TestRes <- cox.zph(CoxModel) # Test the proportional hazards assumption using Schoenfeld residuals
print(TestRes)
plot(TestRes) # Plot the Schoenfeld residuals using base R

# Association between adherence and persistence: fit time varying covariates Cox regression model
CoxTimeVarying <- coxph(formula = Surv(time1, time2, EventPers) ~ cod, data = Compliance)
summary(CoxTimeVarying)

# Number and percentage of subjects who had an event and who are censored, per treatment group
EventSummary <- Persistence %>%
  group_by(IntLabelNum, survstat) %>%
  summarise(n = n()) %>%
  mutate(percentage = (n/sum(n))*100) %>%
  arrange(IntLabelNum, survstat)
print(EventSummary)

# Data for INLA models
data.long <- Compliance[order(Compliance$ID), ]
data.surv <- Persistence[order(Persistence$ID), ]

ParamMap <- data.frame(
  Parameter = c("InteY", 
                "TIMEY",
                "TRT",
                "INTX",
                "AGE",
                "CRISK",
                "HCHOL",
                "precision for qkumar observations",
                "SD of IDY (1)", 
                "SD of IDY (2)",
                "rho of IDY", 
                "Beta for eta", 
                "int", 
                "Age",
                "CARDRISK",
                "HDLCHOL",
                "Precision for baseline.hazard"),
  TblLabel = c("$\\beta_0$", 
               "$\\beta_{\\text{time}}$", 
               "$\\beta_{\\text{group}}$", 
               "$\\beta_{\\text{tx}}$", 
               "$\\beta_{\\text{age}}$", 
               "$\\beta_{\\text{crisk}}$", 
               "$\\beta_{\\text{HDL}}$", 
               "$\\psi$", 
               "$\\sigma_{b_0}$", 
               "$\\sigma_{b_1}$",
               "$\\rho_{{b_0},{b_1}}$", 
               "$\\phi$", 
               "$\\gamma_{\\text{group}}$",
               "$\\gamma_{\\text{age}}$", 
               "$\\gamma_{\\text{crisk}}$", 
               "$\\gamma_{\\text{HDL}}$",
               "$\\tau_{h}$"),
  FigLabel = c("beta[0]", 
               "beta[time]", 
               "beta[group]", 
               "beta[tx]", 
               "beta[age]", 
               "beta[crisk]", 
               "beta[HDL]", 
               "psi",
               "sigma[b0]", 
               "sigma[b1]", 
               "rho[b0, b1]", 
               "phi",
               "gamma[group]", 
               "gamma[age]", 
               "gamma[crisk]",
               "gamma[HDL]",
               "tau[h]"),
  Sort = 1:17
)

######################################################################################################
###FIT THE LINEAR KUMARASWAMY MODEL FOR ADHERENCE JOINTLY WITH COX MODEL FOR PERSISTENCE USING INLA###
######################################################################################################

###############
###MODEL FIT###
###############

gapLongi <- 1/8 # Gap between longitudinal measurements
gap <- 1/80 # Used to generate many time points for permutation algorithm
followup <- 13 # Follow-up time
endpoint <- 16 # Prediction endpoint

# Model fit
NL <- nrow(data.long)
NS <- nrow(data.surv)

data.long$time <- data.long$StudyMonth
data.long$Y <- data.long$cod

data.long <- data.long[, c("ID", "time", "int", "Age", "CARDRISK", "HDLCHOL", "Y")]

data.surv$eventTimes <- data.surv$PersMonth
data.surv$Event <- data.surv$survstat

# Baseline characteristics stats
CalcGrpStats <- function(data, VarName, VarLabel) {
  stats <- data %>%
    group_by(int) %>%
    summarise(
      Mean = mean({{VarName}}, na.rm = TRUE),
      SD = sd({{VarName}}, na.rm = TRUE),
      Min = min({{VarName}}, na.rm = TRUE),
      Median = median({{VarName}}, na.rm = TRUE),
      Max = max({{VarName}}, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(int = as.character(int), Variable = VarLabel) %>%
    pivot_longer(cols = -c(int, Variable), names_to = "Statistic", values_to = "Value") %>%
    arrange(Variable, match(Statistic, c("Mean", "SD", "Min", "Median", "Max")))
  TotStats <- data %>%
    summarise(
      Mean = mean({{VarName}}, na.rm = TRUE),
      SD = sd({{VarName}}, na.rm = TRUE),
      Min = min({{VarName}}, na.rm = TRUE),
      Median = median({{VarName}}, na.rm = TRUE),
      Max = max({{VarName}}, na.rm = TRUE)
    ) %>%
    pivot_longer(cols = everything(), names_to = "Statistic", values_to = "Value") %>%
    mutate(int = "Total", Variable = VarLabel)
  bind_rows(stats, TotStats)
}
AgeStats <- CalcGrpStats(data.surv, Age, "Age (years)")
CardriskStats <- CalcGrpStats(data.surv, CARDRISK, "Cardiac risk score")
HDLStats <- CalcGrpStats(data.surv, HDLCHOL, "HDL cholesterol (mg/dL)")
AllStats <- bind_rows(AgeStats, CardriskStats, HDLStats)
FinalStats <- AllStats %>%
  pivot_wider(
    names_from = int,
    values_from = Value,
    names_prefix = "Group_"
  ) %>%
  select(Variable, Statistic, starts_with("Group_"))
FmtDec <- function(x) {
  sapply(x, function(y) {
    if (floor(y) == y) {
      return(as.character(as.integer(y)))
    } else {
      formatted <- formatC(y, format = "f", digits = 1)
      return(sub("0+$", "", sub("\\.0+$", "", formatted)))
    }
  })
}
FinalStats <- FinalStats %>%
  group_by(Variable) %>%
  mutate(Variable = if_else(Statistic == "Mean", Variable, "")) %>%
  ungroup()
FinalStats <- FinalStats %>%
  mutate(across(starts_with("Group_"), ~FmtDec(as.numeric(.x))))
LaTeXTable <- xtable(FinalStats, 
                     caption = "Descriptive Statistics of Continuous Variables", 
                     label = "tab:descriptiveStats", 
                     align = c("l", "l", "l", "r", "r", "r"),
                     digits = c(0, 0, 0, 3, 3, 3))
print(LaTeXTable, 
      file = "Manuscript/Output/DescriptiveStats.tex",
      include.rownames = FALSE, 
      floating = TRUE, 
      booktabs = TRUE)

# Prepare data for INLA
# Cox model structure for survival (with Bayesian smooth splines for the baseline hazard, i.e., "rw2")
CoxExt <- inla.coxph(YS ~ -1 + Intercept + int + Age + CARDRISK + HDLCHOL, 
                     control.hazard = list(model = "rw2", scale.model = TRUE, diagonal = 1e-2,
                                           constr = TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))),
                     data = list(YS = inla.surv(time = c(data.surv$eventTimes), event = c(data.surv$Event)),
                                 Intercept = rep(1, NS), 
                                 eta = data.surv$ID, 
                                 int = data.surv$int, 
                                 Age = data.surv$Age, 
                                 CARDRISK = data.surv$CARDRISK, 
                                 HDLCHOL = data.surv$HDLCHOL))

# Time weight for time-dependent covariates in survival (i.e., fixed and random slope)
t.weight <- CoxExt$data$baseline.hazard.time + 0.5*CoxExt$data$baseline.hazard.length
nsCox <- nrow(CoxExt$data) # Number of intervals for survival
NY <- NL + nsCox
IDcox <- CoxExt$data$expand..coxph # Random effects ID for survival
CoxExt$data$eta <- 1:nsCox # Unique ID for shared linear predictor Y
CoxExt$data$intX <- CoxExt$data$int*t.weight

# Match binary and continuous covariates with survival time intervals
IDBC <- merge(CoxExt$data, unique(cbind("expand..coxph" = data.long$ID)), by = "expand..coxph")

covariates <- list(
  InteY = rep(1, NY), # Intercept Y
  TIMEY = c(data.long$time, t.weight), # Time Y
  TRT = c(data.long$int, IDBC$int),
  INTX = c(data.long$time*data.long$int, IDBC$intX), # Treatment by time interaction
  AGE = c(data.long$Age, IDBC$Age),
  CRISK = c(data.long$CARDRISK, IDBC$CARDRISK),
  HCHOL = c(data.long$HDLCHOL, IDBC$HDLCHOL),
  IDY = c(data.long$ID, IDcox), # Random intercept Y
  IDY_s = c(NS + data.long$ID, NS + IDcox), # Random slope Y
  WY_s = c(data.long$time, t.weight), # Weight random slope Y
  u1 = c(rep(NA, NL), CoxExt$data$eta), # Y association with survival
  w1 = c(rep(NA, NL), rep(-1, nsCox)) # Y weight association with survival
)

Y.joint <- list(
  Y = c(data.long$Y, rep(NA, nsCox)),
  Y.eta = c(rep(NA, NL), rep(0, nsCox)) # Y association with survival
)

jointdf <- data.frame(covariates, Y.joint)
joint.DataCox <- c(as.list(inla.rbind.data.frames(jointdf, CoxExt$data)), CoxExt$data.list)
Yjoint <- joint.DataCox[c("Y", "Y.eta", "y..coxph")] # Outcomes (longitudinal and survival)
joint.DataCox$Y <- Yjoint

# Update formula from the Cox model structure to add longitudinal part
formulaJ <- update(CoxExt$formula, Yjoint ~ . -1 + InteY + TIMEY + INTX + TRT + AGE + CRISK + HCHOL +
                     f(IDY, model = "iidkd", order = 2, n = NS*2, constr = FALSE,
                       hyper = list(theta = list(param = c(10, 1, 1, 0)))) +
                     f(IDY_s, WY_s, copy ="IDY") +
                     f(u1, w1, model = "iid", hyper = list(prec = list(initial = -6, fixed = TRUE)), constr = FALSE) + # Remove for separate model
                     f(eta, copy = "u1", hyper = list(beta = list(fixed = FALSE, param = c(0, 0.16), initial = 1)))) # Remove for separate model

fitINLAModel <- function(q, qkumarModel = TRUE) {
  
  # INLA function call
  inla.setOption(inla.mode = "experimental")
  if (qkumarModel) {
    JMinla <- inla(verbose = FALSE,
                   formula = formulaJ,
                   family = c("qkumar", "Gaussian", CoxExt$family),
                   data = joint.DataCox,
                   control.fixed = list(mean = 0, prec = 0.16, mean.intercept = 0, prec.intercept = 0.16),
                   control.family = list(
                     list(control.link = list(quantile = q)),
                     list(hyper = list(prec = list(initial = 12, fixed = TRUE))),
                     list()),
                   E = joint.DataCox$E..coxph,
                   control.inla = list(int.strategy = "eb"),
                   control.compute = list(config = TRUE))
  } else {
    JMinla <- inla(verbose = FALSE,
                   formula = formulaJ,
                   family = c("beta", "Gaussian", CoxExt$family),
                   data = joint.DataCox,
                   control.fixed = list(mean = 0, prec = 0.16, mean.intercept = 0, prec.intercept = 0.16),
                   control.family = list(
                     list(control.link = list()),
                     list(hyper = list(prec = list(initial = 12, fixed = TRUE))),
                     list()),
                   E = joint.DataCox$E..coxph,
                   control.inla = list(int.strategy = "eb"),
                   control.compute = list(config = TRUE))
  }
  
  print(summary(JMinla))
  
  fixed <- JMinla$summary.fixed
  hyperpar <- JMinla$summary.hyperpar
  
  fixed <- cbind(Parameter = rownames(fixed), fixed)
  rownames(fixed) <- NULL
  hyperpar <- cbind(Parameter = rownames(hyperpar), hyperpar)
  rownames(hyperpar) <- NULL
  
  hyperpar <- hyperpar %>%
    mutate(Parameter = ifelse(Parameter == "precision parameter for the beta observations", 
                              "precision for qkumar observations", 
                              Parameter))
  
  combined <- bind_rows(fixed, hyperpar)
  
  IIDKD <- inla.iidkd.sample(10^4, JMinla, "IDY")
  
  VarCov <- matrix(unlist(IIDKD), nrow = 2^2)
  VarCovSD <- matrix(apply(VarCov, 1, sd), 2, 2)
  VarCov025 <- matrix(apply(VarCov, 1, function(x) quantile(x,  0.025)), 2, 2)
  VarCov975 <- matrix(apply(VarCov, 1, function(x) quantile(x,  0.975)), 2, 2)
  
  qq <- matrix(rowMeans(matrix(unlist(IIDKD), nrow = 2^2)), 2, 2)
  
  combined[combined$Parameter == "Theta1 for IDY", "Parameter"] <- "SD of IDY (1)"
  combined[combined$Parameter == "Theta2 for IDY", "Parameter"] <- "SD of IDY (2)"
  combined[combined$Parameter == "Theta3 for IDY", "Parameter"] <- "rho of IDY"
  
  combined[combined$Parameter == "SD of IDY (1)", "mean"] <- sqrt(diag(qq)[1])
  combined[combined$Parameter == "SD of IDY (2)", "mean"] <- sqrt(diag(qq)[2])
  combined[combined$Parameter == "rho of IDY", "mean"] <- qq[lower.tri(qq)]
  
  combined[combined$Parameter == "SD of IDY (1)", "sd"] <- VarCovSD[1, 1]
  combined[combined$Parameter == "SD of IDY (2)", "sd"] <- VarCovSD[2, 2]
  combined[combined$Parameter == "rho of IDY", "sd"] <- VarCovSD[1, 2]
  
  combined[combined$Parameter == "SD of IDY (1)", "0.025quant"] <- VarCov025[1, 1]
  combined[combined$Parameter == "SD of IDY (2)", "0.025quant"] <- VarCov025[2, 2]
  combined[combined$Parameter == "rho of IDY", "0.025quant"] <- VarCov025[1, 2]
  
  combined[combined$Parameter == "SD of IDY (1)", "0.975quant"] <- VarCov975[1, 1]
  combined[combined$Parameter == "SD of IDY (2)", "0.975quant"] <- VarCov975[2, 2]
  combined[combined$Parameter == "rho of IDY", "0.975quant"] <- VarCov975[1, 2]
  
  combined <- merge(combined, ParamMap, by = "Parameter", all.x = TRUE)
  combined$q <- q
  
  n_post_sims <- 2000
  psims <- inla.posterior.sample(n_post_sims, JMinla)
  
  save(combined, psims, file = paste0("Manuscript/Output/AdhereKUMJointRObject_", q, ".RData"))
  save(IIDKD, file = paste0("Manuscript/Output/AdhereKUMJointIIDKD_", q, ".RData"))
  
  #######################
  ###MODEL DIAGNOSTICS###
  #######################
  
  observed <- data.long$Y
  n_obs_long <- length(observed)
  pred_raw_prec <- psims |> sapply(\(sim) {
    names(sim$hyperpar)[names(sim$hyperpar) == "precision parameter for the beta observations"] <- "precision for qkumar observations"
    sim$hyperpar["precision for qkumar observations"]
  })
  
  predicted_raw <- psims |> sapply(\(sim) {
    sim$latent[seq_len(n_obs_long)]
  })
  predicted_INLA <- seq_len(n_post_sims) |> sapply(\(sim) {
    rkum(n_obs_long, plogis(predicted_raw[,sim]), pred_raw_prec[sim], q)
  })
  standardised_residual_INLA <- sapply(seq_len(n_obs_long), \(i) {
    mean(predicted_INLA[i, ] <= observed[i])
  })
  yhat_INLA <- rowMeans(predicted_INLA)
  
  n_latent <- nrow(psims[[1]]$latent)
  latent_names <- psims[[1]]$latent |> rownames()
  latent_pos_InteY <- which(latent_names |> startsWith("InteY"))
  latent_pos_Intercept <- which(latent_names |> startsWith("Intercept"))
  # Number of predictors after intercept, time, and binary treatment indicator:
  n_other_preditors <- n_latent - latent_pos_InteY - 2
  
  beta_effects <- psims |> sapply(\(sim) {
    sim$latent[seq(latent_pos_InteY, n_latent),]
  })
  beta_pred <- cbind(rep(1, n_obs_long), 
                     data.long$time,
                     data.long$time*data.long$int,
                     as.matrix(data.long[,3:(ncol(data.long)-1)])) %*% beta_effects
  gamma_effects <- psims |> sapply(\(sim) {
    sim$latent[seq(latent_pos_Intercept, latent_pos_InteY - 1),]
  })
  random_effects_level <- psims |> sapply(\(sim) {
    sim$latent[which(latent_names |> startsWith("IDY"))[seq_len(NS)],]
  })
  random_effects_slope <- psims |> sapply(\(sim) {
    sim$latent[which(latent_names |> startsWith("IDY"))[seq_len(NS) + NS],]
  })
  eta_pred <- beta_pred + random_effects_level[data.long$ID,] + random_effects_slope[data.long$ID,]*data.long$time
  predicted_manual <- seq_len(n_post_sims) |> sapply(\(sim) {
    rkum(n_obs_long, plogis(eta_pred[,sim]), pred_raw_prec[sim], q)
  })
  standardised_residual_manual <- sapply(seq_len(n_obs_long), \(i) {
    mean(predicted_manual[i, ] <= observed[i])
  })
  yhat_manual <- rowMeans(predicted_manual)
  
  QQPlot <- function(std_res, NObs = length(std_res)) {
    par(mar = c(4.5, 5.0, 3.0, 0.5))
    qqplot(seq_len(NObs)/(NObs + 1), std_res, main = "Residual QQ Plot", xlab = "Expected Residual", ylab = "Observed Residual", col = 'darkblue', bty = 'n', cex = 0.5, pch = 16, cex.axis = 1.5, cex.lab = 2.0, cex.main = 2.0)
    lines(c(0, 1), c(0, 1), lwd = 2)
  }
  
  ResPlot <- function(std_res, yhat, NObs = length(std_res)) {
    x <- seq_len(NObs)/(NObs + 1)
    y <- std_res[order(yhat)]
    ylvls <- c(0.25, 0.5, 0.75)
    par(mar = c(4.5, 5.0, 3.0, 0.5))
    plot(x, y, main = "Residuals vs Fits", type = 'n', xlab = "Fitted Value", ylab = "Observed Residual", yaxt = 'n', bty = 'n', cex.axis = 1.5, cex.lab = 2.0, cex.main = 2.0)
    axis(2, at = c(0, ylvls, 1), cex.axis = 1.5)
    abline(h = ylvls, lty = 2, lwd = 2, col = "darkgreen")
    points(x, y, cex = 0.5, pch = 16)
    for (j in seq_along(ylvls)) {
      abline(quantreg::rq(y ~ x, ylvls[j]), col = j + 1, lwd = 4)
    }
  }
  
  standardised_residual_INLA |> QQPlot()
  standardised_residual_INLA |> ResPlot(yhat_INLA)
  
  standardised_residual_manual |> QQPlot()
  standardised_residual_manual |> ResPlot(yhat_manual)
  
  ######################
  ###BY-PATIENT PLOTS###
  ######################
  
  time_seq <- seq(0, endpoint, gap)
  n_time_seq <- length(time_seq)
  int_seq <- rep(1, n_time_seq)
  
  area_centers <- (time_seq[-1] + time_seq[-n_time_seq])/2
  
  rolling_mean <- function(x, k = 15) {
    n <- length(x)
    y <- x
    x <- c(x, rep(x[n], k))
    for (i in 1:n) {
      y[i] <- x[i:(i+k-1)] |> mean(na.rm = TRUE)
    }
    y
  }
  
  if (q == 0.5) {
    
    plot_random_patients <- function(n = 20, q) {
      patient_IDs <- sample(unique(data.surv[data.surv$survstat == 0, ]$ID), size = n, replace = FALSE)
      
      for (sbj in patient_IDs) {
        surv_sbj <- data.surv |> filter(ID == sbj)
        rows_sbj <- which(data.long$ID == sbj)
        nrows_sbj <- length(rows_sbj)
        long_sbj <- data.long |> filter(ID == sbj)
        predicted_sbj <- predicted_INLA[rows_sbj,]
        fitted_sbj <- predicted_raw[rows_sbj,]
        random_effects_level_sbj <- random_effects_level[sbj,]
        random_effects_slope_sbj <- random_effects_slope[sbj,]
        
        sbj_intervals <- fitted_sbj |> plogis() |> apply(1, shortestinterval)
        sbj_pred_intervals <- predicted_sbj |> apply(1, shortestinterval)
        sbj_shade_data <- data.frame(
          Time = long_sbj$time,
          Low = sbj_intervals[1,], 
          High = sbj_intervals[2,]
        )
        plot_line_names <- c('Observation', paste0('Fitted q=', q), 'Lower Limit', 'Upper Limit', 'Counter-factual')
        sbj_plot_data <- data.frame(
          Item = rep(plot_line_names, each = nrows_sbj), 
          Time = rep(long_sbj$time, times = length(plot_line_names)), 
          Value = c(long_sbj$Y, 
                    fitted_sbj |> plogis() |> rowMeans(), 
                    sbj_pred_intervals[1,], 
                    sbj_pred_intervals[2,], 
                    (fitted_sbj +
                       ((1-2*long_sbj$int)*long_sbj$time) %*% beta_effects[3,,drop=FALSE] +
                       (1-2*long_sbj$int) %*% beta_effects[4,,drop=FALSE]
                    ) |> plogis() |> rowMeans()
          )
        )
        
        sbj_plot_data |> ggplot(aes(x = Time)) + 
          geom_ribbon(aes(ymin = Low, ymax = High), data = sbj_shade_data, alpha = 0.2) + 
          geom_line(aes(y = Value, group = Item, colour = Item), data = \(x) filter(x, Item != 'Observation')) + 
          geom_point(aes(y = Value, group = Item, colour = Item), data = \(x) filter(x, Item == 'Observation')) + 
          scale_colour_discrete(breaks = plot_line_names[c(4,1,2,5,3)]) + 
          theme(legend.title = element_blank()) + 
          ylim(0,1)
        
        plot_line_names <- plot_line_names[-5]
        sbj_plot_data <- data.frame(
          Item = rep(plot_line_names, each = nrows_sbj), 
          Time = rep(long_sbj$time, times = length(plot_line_names)), 
          Value = c(long_sbj$Y, 
                    fitted_sbj |> plogis() |> rowMeans(), 
                    sbj_pred_intervals[1,] |> rolling_mean(5), 
                    sbj_pred_intervals[2,] |> rolling_mean(5)
          )
        )
        
        observation_data <- filter(sbj_plot_data, Item == 'Observation')
        upper_limit_data <- filter(sbj_plot_data, Item == 'Upper Limit')
        lower_limit_data <- filter(sbj_plot_data, Item == 'Lower Limit')
        fitted_data <- filter(sbj_plot_data, Item == paste0('Fitted q=', q))
        
        p <- ggplot() + 
          geom_ribbon(data = sbj_shade_data, aes(x = Time, ymin = Low, ymax = High), fill = "grey", alpha = 0.4) +
          geom_line(data = upper_limit_data, aes(x = Time, y = Value), size = 1.3, colour = "black", linetype = "dashed") +
          geom_line(data = lower_limit_data, aes(x = Time, y = Value), size = 1.3, colour = "black", linetype = "dashed") +
          geom_line(data = fitted_data, aes(x = Time, y = Value), size = 1.3, colour = "red") +
          geom_point(data = observation_data, aes(x = Time, y = Value), size = 4, colour = "blue") +
          scale_x_continuous(name = "Day") +
          scale_y_continuous(name = "Adherence", limits = c(0, 1)) +
          theme_minimal() +
          theme(legend.position = "none",
                axis.text = element_text(size = 14),
                axis.title = element_text(size = 16),
                axis.line = element_line(colour = "black"),
                panel.border = element_blank(),
                panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.2),
                panel.grid.minor.y = element_line(colour = "grey95", linewidth = 0.1),
                panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.2),
                panel.grid.minor.x = element_line(colour = "grey95", linewidth = 0.1))
        
        pdf(paste0("Manuscript/Output/AdhereKUMJointProfile_", q, "_", sbj, ".pdf"), width = 8, height = 6)
        print(p)
        dev.off()
        
        sbj_obs_data <- data.frame(
          Item = rep(plot_line_names[1], nrows_sbj), 
          Time = long_sbj$time, 
          Value = long_sbj$Y
        )
        sbj_explanatory_values <- long_sbj[1,3:(ncol(long_sbj)-1)] |> unlist()
        beta_pred_sbj <- (c(int_seq,
                            time_seq,
                            time_seq*sbj_explanatory_values[1],
                            rep(sbj_explanatory_values,
                                each = n_time_seq)
        ) |> matrix(n_time_seq)
        ) %*% beta_effects
        gamma_pred_sbj <- (c(int_seq, 
                             rep(sbj_explanatory_values, 
                                 each = n_time_seq)
        ) |> matrix(n_time_seq)
        ) %*% gamma_effects
        fitted_sbj <- beta_pred_sbj +
          int_seq %*% t(random_effects_level_sbj) +
          time_seq %*% t(random_effects_slope_sbj)
        predicted_sbj <- seq_len(n_post_sims) |> sapply(\(sim) {
          rkum(n_time_seq, plogis(fitted_sbj[,sim]), pred_raw_prec[sim], q)
        })
        sbj_intervals <- fitted_sbj |> plogis() |> apply(1, shortestinterval)
        sbj_pred_intervals <- predicted_sbj |> apply(1, shortestinterval)
        sbj_shade_data <- data.frame(
          Time = time_seq, 
          Low = sbj_intervals[1,], 
          High = sbj_intervals[2,]
        )
        sbj_line_data <- data.frame(
          Item = rep(plot_line_names[-1], each = n_time_seq), 
          Time = rep(time_seq, times = length(plot_line_names[-1])), 
          Value = c(fitted_sbj |> plogis() |> rowMeans(), 
                    sbj_pred_intervals[1,] |> rolling_mean(), 
                    sbj_pred_intervals[2,] |> rolling_mean()
          )
        )
        
        sbj_line_data |> ggplot(aes(x = Time)) + 
          geom_ribbon(aes(ymin = Low, ymax = High), data = sbj_shade_data, alpha = 0.2) + 
          geom_line(aes(y = Value, group = Item, colour = Item)) + 
          geom_point(aes(y = Value, group = Item, colour = Item), data = sbj_obs_data) + 
          scale_colour_discrete(breaks = plot_line_names[c(4,1,2,5,3)]) + 
          theme(legend.title = element_blank()) + 
          ylim(0,1)
        
      }
      
    }
    
    plot_random_patients(n = 20, q = q)
    
  }
  
  #######################################################
  ###PLOT: HAZARD FUNCTION, LONGITUDINAL, AND SURVIVAL###
  #######################################################
  
  bh_items <- latent_names |> startsWith("baseline.hazard") |> which()
  bh_times <- JMinla$summary.random[["baseline.hazard"]]$ID
  bh_sims_raw <- psims |> sapply(\(sim) {
    sim$latent[bh_items,]
  }) |> t()
  bh_gap <- diff(bh_times) |> mean()
  
  sigma <- psims |> sapply(\(sim) {
    sim$hyperpar["Precision for baseline.hazard"]^(-0.5)
  })
  current_time <- bh_times |> last()
  while (current_time < endpoint) {
    current_time <- current_time + bh_gap
    k <- ncol(bh_sims_raw)
    new_sims <- bh_sims_raw[,k]*2 - bh_sims_raw[,k-1] + rnorm(n_post_sims, 0 , sigma)
    bh_sims_raw <- cbind(bh_sims_raw, new_sims)
    bh_times <- c(bh_times, current_time)
  }
  rm(current_time, k, new_sims)
  
  surv_intercept_column <- latent_names |> startsWith("Intercept") |> which()
  bh_sims <- seq_len(n_time_seq) |> sapply(\(i) {
    bh_sims_raw[, match(TRUE, bh_times > time_seq[i]) - 1]
  })
  phi_sims <- psims |> sapply(\(sim) {
    sim$hyperpar["Beta for eta"]
  })
  
  averages <- data.long[,4:(ncol(data.long)-1), drop = FALSE] |> colMeans()
  n_treatments <- 2
  treatments <- seq_len(n_treatments)
  treatment_plot_lines <- paste0("Treatment=", treatments-1)
  beta_binary <- treatments |> lapply(\(tr) {
    (c(int_seq, 
       time_seq, 
       rep(tr - 1, n_time_seq)*time_seq,
       rep(tr - 1, n_time_seq), 
       rep(averages, each = n_time_seq)
    ) |> matrix(n_time_seq)
    ) %*% beta_effects
  })
  predicted <- treatments |> lapply(\(tr) {
    seq_len(n_post_sims) |> sapply(\(sim) {
      rkum(n_time_seq, plogis(beta_binary[[tr]][,sim]), pred_raw_prec[sim], q)
    })
  })
  beta_intervals <- treatments |> lapply(\(tr) {
    beta_binary[[tr]] |> plogis() |> apply(1, shortestinterval)
  })
  pred_intervals <- treatments |> lapply(\(tr) {
    predicted[[tr]] |> apply(1, shortestinterval)
  })
  gamma_pred <- treatments |> lapply(\(tr) {
    (c(int_seq, 
       rep(tr - 1, n_time_seq),
       rep(averages, each = n_time_seq)
    ) |> matrix(n_time_seq)
    ) %*% gamma_effects
  })
  hazard <- treatments |> lapply(\(tr) {
    exp(bh_sims + t(beta_binary[[tr]])*phi_sims + t(gamma_pred[[tr]]))
  })
  h_intervals <- treatments |> lapply(\(tr) {
    hazard[[tr]] |> apply(2, shortestinterval)
  })
  S <- treatments |> lapply(\(tr) {
    exp(-(((hazard[[tr]][,-1] + hazard[[tr]][,-n_time_seq])/2*gap) |> apply(1, cumsum)))
  })
  S_intervals <- treatments |> lapply(\(tr) {
    S[[tr]] |> apply(1, shortestinterval)
  })
  
  long_plot_data <- data.frame(
    Treatment = rep(treatment_plot_lines, each = n_time_seq), 
    Time = rep(time_seq, times = n_treatments),
    Estimate = c(beta_binary |> sapply(rowMeans)) |> plogis(),
    Cred_Lower = c(beta_intervals |> sapply(\(x) x[1,] |> rolling_mean())), 
    Cred_Upper = c(beta_intervals |> sapply(\(x) x[2,] |> rolling_mean())),
    Component = "Longitudinal"
  )
  
  long_plot_data |> ggplot(aes(x = Time, group = Treatment, colour = Treatment, fill = Treatment)) +
    geom_ribbon(aes(ymin = Cred_Lower, ymax = Cred_Upper), alpha = 0.2, linetype = 3) + 
    geom_line(aes(y = Estimate), linewidth = 1.2) + 
    scale_colour_discrete(breaks = treatment_plot_lines) + 
    theme(legend.title = element_blank()) + 
    ylab("Longitudinal Outcome")
  
  pdf(paste0("Manuscript/Output/AdhereKUMJointLong_", q, ".pdf"), width = 8, height = 7)
  plot <- ggplot(long_plot_data, aes(x = Time*30, group = Treatment, colour = factor(Treatment), fill = factor(Treatment))) + 
    geom_ribbon(aes(ymin = Cred_Lower, ymax = Cred_Upper), alpha = 0.2, linetype = "dotted") + 
    geom_line(aes(y = Estimate), linewidth = 1.2) + 
    scale_colour_manual(values = c("Treatment=0" = "red", "Treatment=1" = "blue"), labels = c("No Intervention", "Intervention")) +
    scale_fill_manual(values = c("Treatment=0" = "red", "Treatment=1" = "blue"), labels = c("No Intervention", "Intervention")) + 
    theme(legend.title = element_blank(), 
          legend.position = "bottom",
          legend.background = element_blank(),
          legend.text = element_text(size = 16),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          panel.border = element_blank(),
          panel.grid.major = element_line(colour = "lightgray", linetype = "solid", linewidth = 0.2),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = "white", colour = NA)) + 
    xlab("Day") + 
    ylab("Adherence") +
    scale_x_continuous(breaks = seq(0, 299, 30), limits = c(0, 299)) +
    scale_y_continuous(breaks = seq(0, 1, 0.05), limits = c(0.6, 1))
  
  print(plot)
  dev.off()
  
  h_plot_data <- data.frame(
    Treatment = rep(treatment_plot_lines, each = n_time_seq), 
    Time = rep(time_seq, times = n_treatments), 
    Estimate = c(hazard |> sapply(colMeans)),
    Cred_Lower = c(h_intervals |> sapply(\(x) x[1,] |> rolling_mean())), 
    Cred_Upper = c(h_intervals |> sapply(\(x) x[2,] |> rolling_mean())),
    Component = "Hazard"
  )
  
  h_plot_data |> ggplot(aes(x = Time, group = Treatment, colour = Treatment, fill = Treatment)) + 
    geom_ribbon(aes(ymin = Cred_Lower, ymax = Cred_Upper), alpha = 0.2, linetype = 3) + 
    geom_line(aes(y = Estimate), linewidth = 1.2) + 
    scale_colour_discrete(breaks = treatment_plot_lines) + 
    theme(legend.title = element_blank()) + 
    ylab("Hazard Function") + 
    scale_y_log10()
  
  pdf(paste0("Manuscript/Output/AdhereKUMJointHazard_", q, ".pdf"), width = 8, height = 7)
  plot <- ggplot(h_plot_data, aes(x = Time*30, group = Treatment, colour = factor(Treatment), fill = factor(Treatment))) + 
    geom_ribbon(aes(ymin = Cred_Lower, ymax = Cred_Upper), alpha = 0.2, linetype = "dotted") + 
    geom_line(aes(y = Estimate), linewidth = 1.2) + 
    scale_colour_manual(values = c("Treatment=0" = "red", "Treatment=1" = "blue"), labels = c("No Intervention", "Intervention")) +
    scale_fill_manual(values = c("Treatment=0" = "red", "Treatment=1" = "blue"), labels = c("No Intervention", "Intervention")) + 
    theme(legend.title = element_blank(), 
          legend.position = "bottom",
          legend.background = element_blank(),
          legend.text = element_text(size = 16),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          panel.border = element_blank(),
          panel.grid.major = element_line(colour = "lightgray", linetype = "solid", linewidth = 0.2),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = "white", colour = NA)) + 
    xlab("Day") + 
    ylab("Hazard Function") + 
    scale_y_log10() +
    scale_x_continuous(breaks = seq(0, 299, 30), limits = c(0, 299)) +
    scale_y_continuous(breaks = seq(0, 1, 0.009), limits = c(1e-4, 0.05))
  
  print(plot)
  dev.off()
  
  S_future_line_data <- data.frame(
    Treatment = rep(treatment_plot_lines, each = length(area_centers)), 
    Time = rep(area_centers, times = length(treatment_plot_lines)), 
    Estimate = c(S |> sapply(rowMeans)),
    Cred_Lower = c(S_intervals |> sapply(\(x) x[1,] |> rolling_mean())), 
    Cred_Upper = c(S_intervals |> sapply(\(x) x[2,] |> rolling_mean())),
    Component = "Survival"
  )
  
  S_future_line_data |> ggplot(aes(x = Time, group = Treatment, colour = Treatment, fill = Treatment)) + 
    geom_ribbon(aes(ymin = Cred_Lower, ymax = Cred_Upper), alpha = 0.2, linetype = 3) + 
    geom_line(aes(y = Estimate), linewidth = 1.2) + 
    scale_colour_discrete(breaks = treatment_plot_lines) + 
    theme(legend.title = element_blank()) + 
    ylim(0,1) + ylab("Survival Probability")
  
  pdf(paste0("Manuscript/Output/AdhereKUMJointSurv_", q, ".pdf"), width = 8, height = 7)
  plot <- ggplot(S_future_line_data, aes(x = Time*30, group = Treatment, colour = factor(Treatment), fill = factor(Treatment))) + 
    geom_ribbon(aes(ymin = Cred_Lower, ymax = Cred_Upper), alpha = 0.2, linetype = "dotted") + 
    geom_line(aes(y = Estimate), linewidth = 1.2) + 
    scale_colour_manual(values = c("Treatment=0" = "red", "Treatment=1" = "blue"), labels = c("No Intervention", "Intervention")) +
    scale_fill_manual(values = c("Treatment=0" = "red", "Treatment=1" = "blue"), labels = c("No Intervention", "Intervention")) + 
    theme(legend.title = element_blank(), 
          legend.position = "bottom",
          legend.background = element_blank(),
          legend.text = element_text(size = 16),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          panel.border = element_blank(),
          panel.grid.major = element_line(colour = "lightgray", linetype = "solid", linewidth = 0.2),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = "white", colour = NA)) + 
    xlab("Day") + 
    ylab("Proportion Persistent") +
    scale_x_continuous(breaks = seq(0, 299, 30), limits = c(0, 299)) +
    scale_y_continuous(breaks = seq(0, 1, 0.025), limits = c(0.85, 1))
  
  print(plot)
  dev.off()
  
  data.surv$Treatment <- data.long$int[match(data.surv$ID, data.long$ID)]
  km_fit <- survfit(Surv(eventTimes, Event) ~ Treatment, data = data.surv)
  km_curve_data <- data.frame(
    Time = km_fit$time, 
    KM_est = km_fit$surv, 
    KM_lower = km_fit$lower, 
    KM_upper = km_fit$upper, 
    Treatment = rep(names(km_fit$strata), km_fit$strata),
    Component = "Kaplan-Meier"
  )
  
  km_curve_data |> ggplot(aes(x = Time, group = Treatment, colour = Treatment, fill = Treatment)) + 
    geom_ribbon(aes(ymin = KM_lower, ymax = KM_upper), alpha = 0.2) + 
    geom_line(aes(y = KM_est)) + 
    theme(legend.title = element_blank()) + 
    ylab("Survival Probability")
  
  S_future_line_data |> ggplot(aes(x = Time, 
                                   group = Treatment, 
                                   colour = Treatment, 
                                   fill = Treatment, 
                                   linetype = Component)) + 
    geom_ribbon(aes(ymin = Cred_Lower, ymax = Cred_Upper), alpha = 0.2) + 
    geom_ribbon(aes(ymin = KM_lower, ymax = KM_upper), data = km_curve_data, alpha = 0.2) + 
    geom_line(aes(y = KM_est), data = km_curve_data, linewidth = 1) + 
    geom_line(aes(y = Estimate), linewidth = 1) + 
    scale_colour_discrete(breaks = treatment_plot_lines) + 
    theme(legend.title = element_blank()) + 
    ylim(0,1) + ylab("Survival Probability")
  
  ### Capturing treatment effect
  
  TreatEff <- plogis(beta_binary[[2]]) - plogis(beta_binary[[1]])
  calculate_stats <- \(x) { 
    diff_dens <- density(x)
    diff_interval <- shortestinterval(x)
    c(
      Mean = mean(x), 
      Median = median(x), 
      Mode = diff_dens$x[which.max(diff_dens$y)], 
      Lower = diff_interval[1], 
      Upper = diff_interval[2]
    )}
  TreatEff_stats <- TreatEff |> apply(1, calculate_stats) |> t()
  TreatEff <- data.frame(Time = time_seq, TreatEff_stats)
  S_TreatEff <- S[[2]] - S[[1]]
  S_TreatEff_stats <- S_TreatEff |> apply(1, calculate_stats) |> t()
  S_TreatEff <- data.frame(Time = area_centers, S_TreatEff_stats)
  list(TreatmentEffect = TreatEff) |> 
    openxlsx::write.xlsx(paste0("Manuscript/Output/AdhereKUMJointTrtEffect_", q, ".xlsx"))
  
  summary_points <- seq(1, followup/gap , length = 5)
  TreatEff[summary_points, ] |> knitr::kable(digits = 3)
  S_TreatEff[summary_points[-1], ] |> knitr::kable(digits = 3)
  
  if (q == 0.5) {
    
    #########################
    ###DYNAMIC PREDICTIONS###
    #########################
    
    pred_plot <- function(subjects = 1:2, draw_plots = TRUE) {
      # Get posterior simulations and process them
      psims <- inla.posterior.sample(n_post_sims, JMinla)
      observed <- data.long$Y
      n_obs_long <- length(observed)
      pred_raw_prec <- psims |> sapply(\(sim) {
        sim$hyperpar["precision for qkumar observations"]
      })
      predicted_raw <- psims |> sapply(\(sim) {
        sim$latent[seq_len(n_obs_long)]
      })
      predicted_INLA <- seq_len(n_post_sims) |> sapply(\(sim) {
        rkum(n_obs_long, plogis(predicted_raw[,sim]), pred_raw_prec[sim], q)
      })
      NS <- length(unique(data.surv$ID))
      n_latent <- nrow(psims[[1]]$latent)
      latent_names <- psims[[1]]$latent |> rownames()
      latent_pos_InteY <- which(latent_names |> startsWith("InteY"))
      latent_pos_Intercept <- which(latent_names |> startsWith("Intercept"))
      # Number of predictors after intercept, time, and binary treatment indicator:
      n_other_preditors <- n_latent - latent_pos_InteY - 2
      beta_effects <- psims |> sapply(\(sim) {
        sim$latent[seq(latent_pos_InteY, n_latent),]
      })
      gamma_effects <- psims |> sapply(\(sim) {
        sim$latent[seq(latent_pos_Intercept, latent_pos_InteY - 1),]
      })
      random_effects_level <- psims |> sapply(\(sim) {
        sim$latent[which(latent_names |> startsWith("IDY"))[seq_len(NS)],]
      })
      random_effects_slope <- psims |> sapply(\(sim) {
        sim$latent[which(latent_names |> startsWith("IDY"))[seq_len(NS) + NS],]
      })
      rolling_mean <- function(x, k = 15) {
        n <- length(x)
        y <- x
        x <- c(x, rep(x[n], k))
        for (i in 1:n) {
          y[i] <- x[i:(i+k-1)] |> mean(na.rm = TRUE)
        }
        y  }
      phi_sims <- psims |> sapply(\(sim) {
        sim$hyperpar["Beta for eta"]
      })
      plot_line_names <- c('Observation', paste0('Fitted q=', q), 'Lower Limit', 'Upper Limit')
      S_plot_lines <- c("Upper Limit", paste0('Fitted q=', q), "Lower Limit")
      # Extend timeline beyond fitted portion
      time_seq <- seq(0, endpoint, gap)
      n_time_seq <- length(time_seq)
      sigma <- psims |> sapply(\(sim) {
        sim$hyperpar["Precision for baseline.hazard"]^(-0.5)
      })
      bh_items <- latent_names |> startsWith("baseline.hazard") |> which()
      bh_times <- JMinla$summary.random[["baseline.hazard"]]$ID
      bh_sims_raw <- psims |> sapply(\(sim) {
        sim$latent[bh_items,]
      }) |> t()
      bh_gap <- diff(bh_times) |> mean()
      current_time <- bh_times |> last()
      while (current_time < endpoint) {
        current_time <- current_time + bh_gap
        k <- ncol(bh_sims_raw)
        new_sims <- bh_sims_raw[,k]*2 - bh_sims_raw[,k-1] + rnorm(n_post_sims, 0 , sigma)
        bh_sims_raw <- cbind(bh_sims_raw, new_sims)
        bh_times <- c(bh_times, current_time)
      }
      surv_intercept_column <- latent_names |> startsWith("Intercept") |> which()
      bh_sims <- seq_len(n_time_seq) |> sapply(\(i) {
        bh_sims_raw[, match(TRUE, bh_times > time_seq[i]) - 1]
      })
      area_centers <- (time_seq[-1] + time_seq[-n_time_seq])/2
      plots <- vector('list', length(subjects))
      # Draw plots for each subject
      for (sbj_i in seq_along(subjects)) {
        sbj <- subjects[sbj_i]
        surv_sbj <- data.surv |> filter(ID == sbj)
        rows_sbj <- which(data.long$ID == sbj)
        nrows_sbj <- length(rows_sbj)
        long_sbj <- data.long |> filter(ID == sbj)
        predicted_sbj <- predicted_INLA[rows_sbj,]
        fitted_sbj <- predicted_raw[rows_sbj,]
        random_effects_level_sbj <- random_effects_level[sbj,]
        random_effects_slope_sbj <- random_effects_slope[sbj,]
        int_seq <- rep(1, n_time_seq)
        sbj_obs_data <- data.frame(
          Item = rep(plot_line_names[1], nrows_sbj), 
          Time = long_sbj$time, 
          Value = long_sbj$Y, 
          Component = "Observation"
        )
        sbj_explanatory_values <- long_sbj[1,3:(ncol(long_sbj)-1)] |> unlist()
        beta_pred_sbj <- (c(int_seq, 
                            time_seq, 
                            time_seq*sbj_explanatory_values[1], 
                            rep(sbj_explanatory_values, 
                                each = n_time_seq)
        ) |> matrix(n_time_seq)
        ) %*% beta_effects
        gamma_pred_sbj <- (c(int_seq, 
                             rep(sbj_explanatory_values, 
                                 each = n_time_seq)
        ) |> matrix(n_time_seq)
        ) %*% gamma_effects
        fitted_sbj <- beta_pred_sbj +
          int_seq %*% t(random_effects_level_sbj) +
          time_seq %*% t(random_effects_slope_sbj)
        predicted_sbj <- seq_len(n_post_sims) |> sapply(\(sim) {
          rkum(n_time_seq, plogis(fitted_sbj[,sim]), pred_raw_prec[sim], q)
        })
        sbj_intervals <- fitted_sbj |> plogis() |> apply(1, shortestinterval)
        sbj_pred_intervals <- predicted_sbj |> apply(1, shortestinterval)
        sbj_intervals <- fitted_sbj |> plogis() |> apply(1, shortestinterval)
        sbj_pred_intervals <- predicted_sbj |> apply(1, shortestinterval)
        sbj_shade_data <- data.frame(
          Time = time_seq,
          Low = sbj_intervals[1,], 
          High = sbj_intervals[2,]
        )
        sbj_line_data <- data.frame(
          Item = rep(plot_line_names[-1], each = n_time_seq), 
          Time = rep(time_seq, times = length(plot_line_names[-1])), 
          Value = c(fitted_sbj |> plogis() |> rowMeans(), 
                    sbj_pred_intervals[1,] |> rolling_mean(21), 
                    sbj_pred_intervals[2,] |> rolling_mean(21)
          ), 
          Component = "Longitudinal"
        )
        sbj_line_data <- subset(sbj_line_data, Time <= 299/30)
        sbj_line_data <- sbj_line_data[!grepl("Lower Limit|Upper Limit", sbj_line_data$Item, ignore.case = TRUE), ]
        h_sbj <- exp(bh_sims + t(fitted_sbj)*phi_sims + t(gamma_pred_sbj))
        area_sbj <- (h_sbj[,-1] + h_sbj[,-n_time_seq])/2*gap
        H_sbj <- area_sbj[, (area_centers > surv_sbj$eventTimes)] |> apply(1, cumsum)
        S_sbj <- exp(-H_sbj)
        S_sbj_times <- area_centers[(area_centers > surv_sbj$eventTimes)]
        S_sbj_n_times <- length(S_sbj_times)
        S_sbj_intervals <- S_sbj |> apply(1, shortestinterval)
        S_sbj_line_data <- data.frame(
          Item = rep(S_plot_lines, each = S_sbj_n_times), 
          Time = rep(S_sbj_times, times = length(S_plot_lines)), 
          Value = c(S_sbj_intervals[2,],
                    S_sbj |> apply(1, mean), 
                    S_sbj_intervals[1,]
          ), 
          Component = "Survival"
        )
        LULimDat <- S_sbj_line_data %>% filter(Item %in% c("Lower Limit", "Upper Limit"))
        other_data <- S_sbj_line_data %>% filter(!Item %in% c("Lower Limit", "Upper Limit"))
        plots[[sbj_i]] <- ggplot() + 
          geom_line(data = LULimDat, aes(x = Time*30, y = Value, group = Item), colour = "black", linetype = "dashed", size = 1.0) +
          geom_line(data = other_data, aes(x = Time*30, y = Value, group = Item), colour = "black", linetype = "solid", size = 1.0) +
          geom_line(data = sbj_line_data, aes(x = Time*30, y = Value, group = Item), colour = "maroon", size = 1.0) +
          geom_point(data = sbj_obs_data, aes(x = Time*30, y = Value), colour = "purple", size = 4, shape = 19) +
          geom_vline(xintercept = surv_sbj$eventTimes*30, colour = "black", linetype = "dotted") +
          scale_x_continuous(breaks = seq(0, 400, 30), limits = c(0, 400), name = "Day") +
          scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0.3, 1), name = "Adherence", sec.axis = sec_axis(~., breaks = seq(0, 1, 0.1), name = "Proportion Persistent")) +
          theme(
            legend.position = "none",
            panel.background = element_blank(),
            panel.border = element_blank(),
            panel.grid.major = element_line(colour = "lightgrey", linewidth = 0.2),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 19),
            axis.title.x = element_text(size = 19),
            axis.title.y = element_text(size = 19),
            axis.title.y.right = element_text(size = 19),
            axis.ticks = element_line(colour = "black")
          )
        
        print(plots[[sbj_i]])
      }
      if (draw_plots == TRUE) {
        for (sbj_i in seq_along(subjects)) {
          plots[[sbj_i]] |> print()
          pdf(file = paste0("Manuscript/Output/AdhereKUMJointPred_", q, "_", sbj_i, ".pdf"), width = 8, height = 6)
          print(plots[[sbj_i]])
          dev.off()
        }
      }
      post_sims_list <- list(
        INLA_sims_list = psims,
        pred_raw_prec = pred_raw_prec,
        fitted_raw = predicted_raw, 
        predicted_INLA = predicted_INLA, 
        beta_effects = beta_effects, 
        gamma_effects = gamma_effects, 
        random_effects_level = random_effects_level,
        random_effects_slope = random_effects_slope, 
        phi_sims = phi_sims,
        sigma = sigma, 
        time_seq = time_seq, 
        bh_sims = bh_sims,
        plots = plots
      )
      post_sims_list |> invisible()
    }
    
    subjects_with_lots_of_data <- merge(data.surv, data.long, by = "ID") %>%
      filter(EventPers == 0) %>%
      .$ID %>%
      factor(levels = seq_len(NS)) %>%
      table() %>%
      order(decreasing = TRUE) %>%
      .[1:80]
    
    post_sims_list <- pred_plot(subjects = subjects_with_lots_of_data)
    
  }
  
  return(list(combined = combined, TreatEff = TreatEff, S_TreatEff = S_TreatEff))
  
}

#######################################################################
###PLOT: PARAMETER ESTIMATES AND CONFIDENCE INTERVALS OVER QUANTILES###
#######################################################################

result_q999 <- fitINLAModel(0.999, qkumarModel = FALSE)
combined_q999 <- result_q999$combined
TreatEff_q999 <- result_q999$TreatEff
S_TreatEff_q999 <- result_q999$S_TreatEff

result_q10 <- fitINLAModel(0.1)
combined_q10 <- result_q10$combined
TreatEff_q10 <- result_q10$TreatEff
S_TreatEff_q10 <- result_q10$S_TreatEff

result_q25 <- fitINLAModel(0.25)
combined_q25 <- result_q25$combined
TreatEff_q25 <- result_q25$TreatEff
S_TreatEff_q25 <- result_q25$S_TreatEff

result_q50 <- fitINLAModel(0.5)
combined_q50 <- result_q50$combined
TreatEff_q50 <- result_q50$TreatEff
S_TreatEff_q50 <- result_q50$S_TreatEff

result_q75 <- fitINLAModel(0.75)
combined_q75 <- result_q75$combined
TreatEff_q75 <- result_q75$TreatEff
S_TreatEff_q75 <- result_q75$S_TreatEff

result_q95 <- fitINLAModel(0.95)
combined_q95 <- result_q95$combined
TreatEff_q95 <- result_q95$TreatEff
S_TreatEff_q95 <- result_q95$S_TreatEff

combined <- bind_rows(combined_q10,
                      combined_q25,
                      combined_q50,
                      combined_q75,
                      combined_q95,
                      combined_q999)

combined <- combined %>% filter(Parameter != "Intercept")

combined <- combined %>% arrange(Sort)

LongData <- combined %>%
  pivot_longer(cols = c("mean", "0.025quant", "0.975quant"),
               names_to = "statistic",
               values_to = "value") %>%
  mutate(quant_stat = paste(statistic, "for_q", q, sep = "_")) %>%
  select(-q, -statistic)

WideData <- LongData %>%
  pivot_wider(names_from = quant_stat, values_from = value)

LongData <- WideData %>%
  group_by(Parameter, TblLabel, FigLabel, Sort) %>%
  summarize(across(starts_with("mean_for_q_"), mean, na.rm = TRUE),
            across(starts_with("0.025quant_for_q_"), mean, na.rm = TRUE),
            across(starts_with("0.975quant_for_q_"), mean, na.rm = TRUE), 
            .groups = 'drop')

print(LongData)

FormatBCI <- function(lower, upper) {
  LowerFormat <- sprintf("%0.2f", lower)
  UpperFormat <- sprintf("%0.2f", upper)
  return(c(LowerFormat, UpperFormat))
}

qValues <- sort(unique(combined$q))
for (qVal in qValues) {
  LowerCol <- paste0("0.025quant_for_q_", qVal)
  UpperCol <- paste0("0.975quant_for_q_", qVal)
  BCILowerCol <- paste0("BCI_lower_for_q_", qVal)
  BCIUpperCol <- paste0("BCI_upper_for_q_", qVal)
  
  BCIVal <- mapply(FormatBCI, LongData[[LowerCol]], LongData[[UpperCol]], SIMPLIFY = FALSE)
  
  LongData[[BCILowerCol]] <- paste0("[", sapply(BCIVal, `[`, 1), ";")
  LongData[[BCIUpperCol]] <- paste0(sapply(BCIVal, `[`, 2), "]")
}

for (qVal in qValues) {
  BlankCol <- paste0("blank_for_q_", qVal)
  LongData[[BlankCol]] <- ""
}

OrderCol <- c("Sort", unlist(lapply(qValues, function(q) {
  c(paste0("mean_for_q_", q), paste0("BCI_lower_for_q_", q), paste0("BCI_upper_for_q_", q), paste0("blank_for_q_", q))
})))
OrderCol <- OrderCol[1:(length(OrderCol) - 1)]

CombRes <- LongData %>%
  select(TblLabel, all_of(OrderCol)) %>%
  arrange(Sort) %>%
  select(-Sort) %>%
  rename(Parameter = TblLabel)

NCol <- grep("mean_for_q_", names(CombRes))
CombRes[NCol] <- lapply(CombRes[NCol], function(x) round(as.numeric(x), 2))

DecPlace <- c(0,
              Parameter = 0, 
              mean_for_q_0.10 = 2, BCI_lower_for_q_0.10 = 2, BCI_upper_for_q_0.10 = 2, blank_for_q_0.10 = 0, 
              mean_for_q_0.25 = 2, BCI_lower_for_q_0.25 = 2, BCI_upper_for_q_0.25 = 2, blank_for_q_0.25 = 0, 
              mean_for_q_0.5 = 2, BCI_lower_for_q_0.5 = 2, BCI_upper_for_q_0.5 = 2, blank_for_q_0.5 = 0, 
              mean_for_q_0.75 = 2, BCI_lower_for_q_0.75 = 2, BCI_upper_for_q_0.75 = 2, blank_for_q_0.75 = 0, 
              mean_for_q_0.95 = 2, BCI_lower_for_q_0.95 = 2, BCI_upper_for_q_0.95 = 2)

LaTeXTable <- xtable(select(CombRes, !c(blank_for_q_0.95, mean_for_q_0.999, BCI_lower_for_q_0.999, BCI_upper_for_q_0.999)), digits = DecPlace)

print.xtable(LaTeXTable, 
             file = "Manuscript/Output/AdhereKUMJointQuantiles.tex", 
             floating = FALSE, 
             include.rownames = FALSE, 
             include.colnames = FALSE,
             sanitize.text.function = identity,
             hline.after = c(-1, 0, nrow(LaTeXTable)),
             only.contents = TRUE)
print(LaTeXTable)

PrepMerge <- function(df, quantile) {
  df %>%
    select(Time, Mean, Lower, Upper) %>%
    mutate(
      Mean = ifelse(round(100*Mean, 2) == 0, "0.00", sprintf("%.2f", round(100*Mean, 2))),
      Lower = ifelse(round(100*Lower, 2) == 0, "0.00", sprintf("%.2f", round(100*Lower, 2))),
      Upper = ifelse(round(100*Upper, 2) == 0, "0.00", sprintf("%.2f", round(100*Upper, 2))),
      CI = paste0("[", Lower, "; & ", Upper, "] &")
    ) %>%
    rename(
      !!paste0("Mean_q", quantile) := Mean,
      !!paste0("CI_q", quantile) := CI
    ) %>%
    select(-Lower, -Upper)
}

# Treatment effect of longitudinal model

TreatEff_q10Prep <- PrepMerge(TreatEff_q10, "10")
TreatEff_q25Prep <- PrepMerge(TreatEff_q25, "25")
TreatEff_q50Prep <- PrepMerge(TreatEff_q50, "50")
TreatEff_q75Prep <- PrepMerge(TreatEff_q75, "75")
TreatEff_q95Prep <- PrepMerge(TreatEff_q95, "95")

TrtCombMerged <- reduce(list(TreatEff_q10Prep,
                             TreatEff_q25Prep,
                             TreatEff_q50Prep,
                             TreatEff_q75Prep,
                             TreatEff_q95Prep),
                        full_join, by = "Time")

FiltTreatEff <- TrtCombMerged %>%
  filter(Time %in% 0:10) %>%
  mutate(Time = as.integer(round(Time*30)))

LaTeXTable <- xtable(FiltTreatEff)

print.xtable(LaTeXTable, 
             include.rownames = FALSE, include.colnames = FALSE,
             file = "Manuscript/Output/AdhereKUMJointLongTrtEffect.tex", 
             floating = FALSE,
             sanitize.text.function = function(x) {x},
             hline.after = c(-1, 0, nrow(LaTeXTable)),
             only.contents = TRUE,
             booktabs = TRUE)

# Treatment effect of survival model

S_TreatEff_q10Prep <- PrepMerge(S_TreatEff_q10, "10")
S_TreatEff_q25Prep <- PrepMerge(S_TreatEff_q25, "25")
S_TreatEff_q50Prep <- PrepMerge(S_TreatEff_q50, "50")
S_TreatEff_q75Prep <- PrepMerge(S_TreatEff_q75, "75")
S_TreatEff_q95Prep <- PrepMerge(S_TreatEff_q95, "95")

S_TrtCombMerged <- reduce(list(S_TreatEff_q10Prep,
                               S_TreatEff_q25Prep,
                               S_TreatEff_q50Prep,
                               S_TreatEff_q75Prep,
                               S_TreatEff_q95Prep),
                          full_join, by = "Time")

S_FiltTreatEff <- S_TrtCombMerged %>%
  mutate(Time = round((Time - 0.00625), 5)) %>%
  mutate(IntCheck = as.integer(Time)) %>%
  filter(Time == IntCheck) %>%
  mutate(Time = Time*30) %>%
  filter(Time <= 300) %>%
  select(-IntCheck)

S_LaTeXTable <- xtable(S_FiltTreatEff)

print.xtable(S_LaTeXTable, 
             include.rownames = FALSE, include.colnames = FALSE,
             file = "Manuscript/Output/AdhereKUMJointSurvTrtEffect.tex", 
             floating = FALSE,
             sanitize.text.function = function(x) {x},
             hline.after = c(-1, 0, nrow(LaTeXTable)),
             only.contents = TRUE,
             booktabs = TRUE)

###FIGURES###

SelectParam <- c("InteY", "TIMEY", "TRT", "INTX", "AGE", "CRISK", "HCHOL", "Beta for eta", "int", "Age", "CARDRISK", "HDLCHOL")

FiltCombined <- combined %>%
  filter(Parameter %in% SelectParam)

MeanEstimates <- combined %>%
  filter(q == 0.999 & Parameter %in% SelectParam) %>%
  mutate(q = 0.5)

OrderFigLab <- ParamMap$FigLabel[order(ParamMap$Sort)]
FiltCombined$FigLabel <- factor(FiltCombined$FigLabel, levels = OrderFigLab)
MeanEstimates$FigLabel <- factor(MeanEstimates$FigLabel, levels = OrderFigLab)

plot <- ggplot(subset(FiltCombined, q != 0.999), aes(x = q, y = mean, group = FigLabel)) +
  geom_line(size = 1.3, linetype = 1, color = "blue") +
  geom_point(shape = 17, size = 5, color = "blue") +
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.32, fill = "steelblue") +
  geom_line(aes(y = `0.025quant`), color = "blue", linetype = "dashed", size = 1.3) +
  geom_line(aes(y = `0.975quant`), color = "blue", linetype = "dashed", size = 1.3) +
  geom_point(data = MeanEstimates, aes(x = q, y = mean), shape = 15, size = 5, color = "darkred") +
  geom_errorbar(data = MeanEstimates, aes(x = q, ymin = `0.025quant`, ymax = `0.975quant`), width = 0.1, color = "darkred", size = 1.2) +
  facet_wrap(~ FigLabel, ncol = 3, scales = 'free', labeller = label_parsed, strip.position = "left") +
  ylab("Parameter") +
  scale_x_continuous(name = "Quantile", breaks = unique(FiltCombined$q)) +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_blank(),
    strip.placement = "outside",
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = 16),
    strip.background = element_blank(),
    panel.grid.major = element_line(size = 0.2, linetype = 'dotted', colour = "grey"),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 14),
    legend.position = 'bottom',
    panel.spacing = unit(0.5, "lines"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  )

pdf("Manuscript/Output/AdhereKUMJointQuantiles.pdf", width = 12, height = 11)
print(plot)
dev.off()

svg("Manuscript/Output/AdhereKUMJointQuantiles.svg", width = 12, height = 11)
print(plot)
dev.off()