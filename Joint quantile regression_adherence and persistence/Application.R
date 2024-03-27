# App to demonstrate dynamic predictions for joint model
# 
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

options(scipen = 12)

# Load packages
c('shiny', 'rhandsontable', 'tidyverse', 'plotly', 'INLA', 'PermAlgo', 'mvtnorm', 'survival') |>
  sapply(\(lib) { 
    library(package = lib, verbose = FALSE, quietly = TRUE, warn.conflicts = FALSE, character.only = TRUE)
  }) -> load_success

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Dynamically updating subject predictions"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          actionButton("predictbutton", "Fit and predict"), 
          radioButtons("groupbuttons", "Group", 
                       c("Treatment" = 1, "Control" = 0), 
                       inline = TRUE, width = '100%'), 
          sliderInput("predictor", "Standardised predictor", -5, 5, 0, 0.01, -2, width = '100%'), 
          rHandsontableOutput("adherencetable", height = 600),
          width = 3
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotlyOutput("jointplot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  showNotification(ui='Tip 1', action='Keep an eye out for messages popping up in this corner, they will let you know if something goes right or wrong.', duration = 8, type='message')
  showNotification(ui='Tip 2', action='Start by entering some adherence data and clicking on the Predict button.', duration = 30, type='message')
    
  set.seed(123)
  
  shortestinterval <- function(postsims, width=0.95) { # Coded by Sean van der Merwe, UFS
    postsims |> sort() -> sorted.postsims
    round(length(postsims)*width) -> gap
    sorted.postsims |> diff(gap) |> which.min() -> pos
    sorted.postsims[c(pos, pos + gap)] }
  
  theme_set(theme_bw())
  
  # Define simulation parameters
  ParamList <- list(
    nsujet = 500, # Number of individuals
    q = 0.5, # Quantile of longitudinal measurements
    # Continuous proportion data
    beta_0 = 2.5, # Intercept
    beta_1 = -0.3, # Slope
    beta_2 = -0.8, # Binary covariate
    beta_3 = 0.1, # Continuous covariate
    psi = 0.9, # Kumaraswamy scale parameter
    SD1 = 0.4, # SD of random intercept
    SD2 = 0.4, # SD of random slope
    rho = 0.2, # Random intercept/slope correlation coefficient
    # Survival data
    phi = -0.3, # Effect of Y's linear predictor on the risk of event
    gamma_1 = c(-0.5), # Binary covariate
    gamma_2 = c(0.1), # Continuous covariate
    gapLongi = 1/8, # Gap between longitudinal measurements
    gap = 1/80, # Used to generate many time points for permutation algorithm
    followup = 10, # Follow-up time
    endpoint = 15 # Prediction endpoint
  )
  
  # Data generation and model fitting functions

  # Function to simulate from qkumar distribution
  rkum <- function(n, Kappa, psi, q) {
    beta <- log(1 - q)/log(1 - exp(-psi)) # Solve beta of qkumar
    alpha <- log(1 - (1 - q)^(1/beta))/log(Kappa) # Solve alpha of qkumar
    (1 - (1 - runif(n))^(1/beta))^(1/alpha) # Transform using inverse CDF
  }

  # This is a function that takes the parameter list as input and produces data sets
  gendata <- function(p) {
    Covar <- p$rho*p$SD1*p$SD2
    Sigma <- matrix(c(p$SD1^2, Covar, Covar, p$SD2^2), 2) # Random effects covariance matrix
    mestime <- seq(0, p$followup, p$gap) # Measurement times
    timesLongi <- mestime[sapply(seq(0, p$followup, p$gapLongi), \(timepoint) 
                                 which.min(abs(mestime - timepoint)))] # Visit times
    time <- rep(mestime, p$nsujet) # Time column
    nmesindiv <- p$followup/p$gap + 1 # Max number of individual measurements
    nmesy <- nmesindiv*p$nsujet # Max total number of longitudinal measurements
    idY <- rep(1:p$nsujet, each = nmesindiv) # Individual ID
    MVnorm <- rmvnorm(p$nsujet, rep(0, 2), Sigma)
    b_int <- rep(MVnorm[, 1], each = nmesindiv) # Random intercept Y
    b_slo <- rep(MVnorm[, 2], each = nmesindiv) # Random slope Y
    binX <- rbinom(p$nsujet, 1, 0.5) # Binary covariate
    ctsX <- rnorm(p$nsujet, 1, 0.5) # Continuous covariate
    # Replicate binX and ctsX for each time point of each subject
    rep_binX <- rep(binX, each = nmesindiv)
    rep_ctsX <- rep(ctsX, each = nmesindiv)
    
    # Linear predictors
    linPredY <- (p$beta_0 + b_int) + (p$beta_1 + b_slo)*time + p$beta_2*binX[idY] + p$beta_3*ctsX[idY]
    
    # Continuous outcome Y
    Y <- rkum(nmesy, Kappa = plogis(linPredY), psi = p$psi, q = p$q) # Logit link usage
    
    # Permutation algorithm to generate survival times dependent on linear predictors
    DatTmp <- permalgorithm(p$nsujet, nmesindiv, 
                            Xmat = matrix(c(linPredY, rep_binX, rep_ctsX), nrow = p$nsujet*nmesindiv),
                            eventRandom = round(rexp(p$nsujet, 0.003) + 1, 0), # ~40% events
                            censorRandom = runif(p$nsujet, 1, nmesindiv), # Uniform random censoring
                            XmatNames = c("linPredY", "binPred", "otherPred"), # Association
                            betas = c(p$phi, p$gamma_1, p$gamma_2) # Association and hazard regression parameters
    )
    # Extract last line for each ID (= event/censoring time)
    survDat <- DatTmp |> 
      group_by(Id) |> 
      summarise(
        Id = last(Id), 
        eventTimes = mestime[last(Stop) + 1], 
        Event = last(Event)
      ) |> as.data.frame()
    survDat$binX <- binX[survDat$Id]  # Align binX with survival data
    survDat$ctsX <- ctsX[survDat$Id]  # Align ctsX with survival data
    DatTmp <- DatTmp |> mutate(
      time = mestime[Start + 1], # Measurement time of the longitudinal outcomes
      Uid = paste(Id, time) # Unique identifier to match covariates and longitudinal outcomes
    )
    longDat3 <- merge(DatTmp[, c("Uid", "Id", "time")], 
                      data.frame(Uid = paste(idY, time), 
                                 binX = rep_binX, 
                                 ctsX = rep_ctsX, 
                                 Y = Y), by = "Uid")
    longDat <- sapply(longDat3[longDat3$time %in% timesLongi, -1], as.numeric)
    longDat <- as.data.frame(longDat[order(longDat[, "Id"], longDat[, "time"]), ])
    list(survDat = survDat, longDat = longDat)
  }
  
  fitINLA <- function(d, q = 0.5) {
    survDat <- d$survDat
    longDat <- d$longDat
    
    # Model fit
    NL <- nrow(longDat)
    NS <- nrow(survDat)
    
    # Prepare data for INLA
    # Cox model structure for survival (with Bayesian smooth splines for the baseline hazard, i.e., "rw2")
    CoxExt <- inla.coxph(YS ~ -1 + Intercept + binX + ctsX, 
                         control.hazard = list(model = "rw2", scale.model = TRUE, 
                                               diagonal = 1e-2, constr = TRUE, 
                                               hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))),
                         data = list(YS = inla.surv(time = c(survDat$eventTimes), event = c(survDat$Event)),
                                     Intercept = rep(1, NS), 
                                     eta = survDat$Id,
                                     binX = survDat$binX, 
                                     ctsX = survDat$ctsX
                         )
    )
    
    # Time weight for time-dependent covariates in survival (i.e., fixed and random slope)
    t.weight <- CoxExt$data$baseline.hazard.time + 0.5*CoxExt$data$baseline.hazard.length
    nsCox <- nrow(CoxExt$data) # Number of intervals for survival
    NY <- NL + nsCox
    IDcox <- CoxExt$data$expand..coxph # Random effects ID for survival
    CoxExt$data$eta <- 1:nsCox # Unique ID for shared linear predictor Y
    
    # Replicate binX and ctsX according to structure of CoxExt$data
    Rep_binX <- survDat$binX[CoxExt$data$expand..coxph]
    Rep_ctsX <- survDat$ctsX[CoxExt$data$expand..coxph]
    
    # Merge replicated covariates with CoxExt$data
    IDBC <- data.frame(CoxExt$data, binX = Rep_binX, ctsX = Rep_ctsX)
    
    covariates <- data.frame(
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
      Y = c(longDat$Y, rep(NA, nsCox)), # Y
      Y.eta = c(rep(NA, NL), rep(0, nsCox)) # Y association with survival
    )
    
    jointdf <- data.frame(covariates, Y.joint)
    joint.DataCox <- c(as.list(inla.rbind.data.frames(jointdf, CoxExt$data)), CoxExt$data.list)
    Yjoint <- joint.DataCox[c("Y", "Y.eta", "y..coxph")] # Outcomes (longitudinal and survival)
    joint.DataCox$Y <- Yjoint
    
    # Update formula from the Cox model structure to add longitudinal part
    formulaJ <- update(CoxExt$formula, Yjoint ~ . -1 + InteY + TIMEY + binXY + ctsXY +
                         f(IDY, model = "iid2d", n = NS*2, constr = FALSE,
                           hyper = list(theta = list(param = c(10, 1, 1, 0)))) +
                         f(IDY_s, WY_s, copy ="IDY") +
                         f(u1, w1, model = "iid", 
                           hyper = list(prec = list(initial = -6, fixed = TRUE)), 
                           constr = FALSE) +
                         f(eta, copy = "u1", 
                           hyper = list(beta = list(fixed = FALSE, param = c(0, 0.16), 
                                                    initial = 1))))
    
    # INLA function call
    inla.setOption(inla.mode = "experimental")
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
                   control.compute=list(config = TRUE)
    )
    JMinla
  }
  
  d <- gendata(ParamList)
  
  pred_plot <- function(d, INLAfit, parameters, sbj = 501, n_post_sims = 1000) {
    # Get posterior simulations and process them
    psims <- inla.posterior.sample(n_post_sims, INLAfit)
    observed <- d$longDat$Y
    n_obs_long <- length(observed)
    pred_raw_prec <- psims |> sapply(\(sim) {
      sim$hyperpar["precision for qkumar observations"]
    })
    predicted_raw <- psims |> sapply(\(sim) {
      sim$latent[seq_len(n_obs_long)]
    })
    predicted_INLA <- seq_len(n_post_sims) |> sapply(\(sim) {
      rkum(n_obs_long, plogis(predicted_raw[,sim]), pred_raw_prec[sim], ParamList$q)
    })
    NS <- d$survDat$Id |> unique() |> length()
    n_latent <- nrow(psims[[1]]$latent)
    latent_names <- psims[[1]]$latent |> rownames()
    latent_pos_InteY <- which(latent_names |> startsWith("InteY"))
    latent_pos_Intercept <- which(latent_names |> startsWith("Intercept"))
    # Number of predictors after intercept, time, and binary treatment indicator:
    n_other_preditors <- n_latent - latent_pos_InteY - 2
    beta_effects <- psims |> sapply(\(sim) {
      sim$latent[seq(latent_pos_InteY, n_latent),]
    })
    beta_pred <- cbind(rep(1, n_obs_long), as.matrix(d$longDat[,2:(ncol(d$longDat)-1)])) %*% beta_effects
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
    plot_line_names <- c('Observation', paste0('Fitted q=', ParamList$q), 'Lower Limit', 'Upper Limit')
    S_plot_lines <- c("Upper Limit", "Fitted q=0.5", "Lower Limit")
    # Extent timeline beyond fitted portion
    time_seq <- seq(0, ParamList$endpoint, ParamList$gap)
    n_time_seq <- length(time_seq)
    sigma <- psims |> sapply(\(sim) {
      sim$hyperpar["Precision for baseline.hazard"]^(-0.5)
    })
    bh_items <- latent_names |> startsWith("baseline.hazard") |> which()
    bh_times <- INLAfit$summary.random[["baseline.hazard"]]$ID
    bh_sims_raw <- psims |> sapply(\(sim) {
      sim$latent[bh_items,]
    }) |> t()
    bh_gap <- diff(bh_times) |> mean()
    current_time <- bh_times |> last()
    while (current_time < ParamList$endpoint) {
      current_time <- current_time + bh_gap
      k <- ncol(bh_sims_raw)
      new_sims <- bh_sims_raw[,k]*2 - bh_sims_raw[,k-1] + rnorm(n_post_sims, 0 , sigma)
      bh_sims_raw <- cbind(bh_sims_raw, new_sims)
      bh_times <- c(bh_times, current_time)
    }
    surv_intercept_column <- latent_names |> startsWith("Intercept") |> which()
    bh_intercept <- psims |> sapply(\(sim) {
      sim$latent[surv_intercept_column,]
    })
    bh_sims <- seq_len(n_time_seq) |> sapply(\(i) {
      bh_sims_raw[, match(TRUE, bh_times > time_seq[i]) - 1] + bh_intercept
    })
    area_centers <- (time_seq[-1] + time_seq[-n_time_seq])/2
    # Draw plot
    surv_sbj <- d$survDat |> filter(Id == sbj)
    rows_sbj <- which(d$longDat$Id == sbj)
    nrows_sbj <- length(rows_sbj)
    long_sbj <- d$longDat |> filter(Id == sbj)
    predicted_sbj <- predicted_INLA[rows_sbj,]
    fitted_sbj <- predicted_raw[rows_sbj,]
    random_effects_level_sbj <- random_effects_level[sbj,]
    random_effects_slope_sbj <- random_effects_slope[sbj,]
    int_seq <- rep(1, n_time_seq)
    sbj_obs_data <- data.frame(
      Item = rep(plot_line_names[1], nrows_sbj), 
      Time = long_sbj$time, 
      Value = long_sbj$Y, 
      Component = "Observations"
    )
    sbj_explanatory_values <- long_sbj[1,3:(ncol(long_sbj)-1)] |> unlist()
    beta_pred_sbj <- (c(int_seq, 
                        time_seq, 
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
      rkum(n_time_seq, plogis(fitted_sbj[,sim]), pred_raw_prec[sim], ParamList$q)
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
    h_sbj <- exp(bh_sims + t(fitted_sbj)*phi_sims + t(gamma_pred_sbj))
    area_sbj <- (h_sbj[,-1] + h_sbj[,-n_time_seq])/2*ParamList$gap
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
    joint_plot <- sbj_line_data |> ggplot(aes(x = Time)) + 
      geom_ribbon(aes(ymin = Low, ymax = High), 
                  data = sbj_shade_data, 
                  alpha = 0.2, fill = '#BB2200FF') + 
      geom_line(aes(y = Value, group = Item, colour = Item, linetype = Component)) + 
      geom_point(aes(y = Value, group = Item, colour = Item), data = sbj_obs_data) + 
      scale_colour_discrete(breaks = plot_line_names[c(4,1,2,5,3)]) + 
      theme(legend.title = element_blank()) + 
      ylim(0,1) + 
      geom_line(aes(y = Value, group = Item, colour = Item, linetype = Component), 
                data = S_sbj_line_data) + 
      geom_vline(xintercept = surv_sbj$eventTimes)

    joint_plot |> ggplotly()
  }
  
  num_periods_total <- 120
  new_num_times <- ceiling(runif(1,0,40))
  new_times <- seq(0, (new_num_times-1))*ParamList$gapLongi
  RVs <- reactiveValues(DF = data.frame(
    Period = seq_len(num_periods_total), 
    Proportion = c(rkum(new_num_times, 0.9 - (new_times)*0.07, ParamList$psi, ParamList$q), 
                   rep(NA_real_, num_periods_total - new_num_times))
  ))
                        
  plot_object <- eventReactive(input$predictbutton, {
    showNotification(ui='Calculating', action='Please allow time to refit model.', duration = 10, type='message')
    new_props <- as.numeric(as.data.frame(hot_to_r(input$adherencetable))[[2]])
    new_props[new_props <= 0] <- NA_real_; new_props[new_props >= 1] <- NA_real_
    new_num_times <- max(which(!is.na(new_props)))
    new_times <- seq(0, (new_num_times-1))*ParamList$gapLongi
    dnew <- d
    newId <- max(d$survDat$Id) + 1
    selected_treatment <- as.numeric(input$groupbuttons)
    selected_predictor <- as.numeric(input$predictor)
    dnew$longDat <- dnew$longDat |> rbind(data.frame(
      Id = newId, 
      time = new_times, 
      binX = selected_treatment, 
      ctsX = selected_predictor, 
      Y = new_props[seq_len(new_num_times)]
    ))
    dnew$survDat <- dnew$survDat |> rbind(data.frame(
      Id = newId, 
      eventTimes = max(new_times) + (ParamList$gapLongi/2), 
      Event = 0,
      binX = selected_treatment, 
      ctsX = selected_predictor
    ))
    new_inla <- fitINLA(dnew, q = ParamList$q)
    showNotification(ui='Calculating', action='Predicting and plotting.', duration = 8, type='message')
    dnew |> pred_plot(new_inla, ParamList, sbj = newId)
  })
  
  output$jointplot <- renderPlotly({
    plot_object()
  })
  
  output$adherencetable <- renderRHandsontable({
    mytable <- rhandsontable(RVs$DF, rowHeaders = NULL, width = "100%", height = 580, colHeaders = c('Time Period', 'Proportion'))
    mytable <- hot_col(hot=mytable, col=1, type='numeric', format='0', readOnly = TRUE)
    mytable <- hot_col(hot=mytable, col=2, type='numeric', format='0.000', readOnly = FALSE)
    mytable
  })
  
  # Manual input validation:
  # (currently not working, relying on automatic validation which ensures numeric only)
  # observeEvent(input$adherencetable, {
  #   newtable <- as.data.frame(hot_to_r(input$adherencetable))
  #   oldtable <- RVs$DF
  #   changes <- which(newtable!=oldtable)
  #   if (length(changes) > 0) { if (changes[1]==1) { changes <- changes[-1] } }
  #   if (length(changes) > 0) {
  #     i <- changes[1] %% num_periods_total
  #     v <- newtable[i,2] |> as.numeric()
  #     valid_value <- FALSE
  #     if (is.finite(v)) {
  #       if ((v > 0) && (v < 1)) {
  #         valid_value <- TRUE
  #       }
  #     }
  #     if (valid_value) {
  #       oldtable[i, 2] <- round(v, 4)
  #     } else {
  #       showNotification(ui='Not a proportion', action='Please enter a proportion in the form 0.####, with no spaces or commas.', duration = 10, type='error')
  #       if (oldtable[num_periods_total,1] == num_periods_total) {
  #         oldtable[num_periods_total,1] <- num_periods_total + 0.1
  #       } else {
  #         oldtable[num_periods_total,1] <- num_periods_total
  #       }
  #     }
  #   }
  #   RVs$DF <<- oldtable
  # })
  # 
  
}

# Run the application 
shinyApp(ui = ui, server = server)
