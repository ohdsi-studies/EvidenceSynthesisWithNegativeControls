library(dplyr)
library(survival)

# Simulation settings ------------------------------------------------------------------------------
#' Create simulation settings
#'
#' @param esSettings            A settings object created either by `EvidenceSynthesis::createSimulationSettings()`
#'                              or `EvidenceSynthesis::createSccsSimulationSettings()`. The effect 
#'                              size settings (e.g. `hazardRatio` and `randomEffectSd`) will not 
#'                              apply to negative controls.
#' @param nNegativeControls     Number of negative controls.
#' @param nOutcomesOfInterest   Number of outcomes of interest.
#' @param nCaptureProcessChars  Number of data capture process characteristics.
#' @param cpcConsistency        How much each database's capture process characteristics will
#'                              differ. Set to 0.1 for low consistency, or 25 for high consistency.
#' @param bias0                 Baseline bias.
#' @param biasCpcSd             SD of the bias associated with each data capture characteristic.
#' @param minDatabasePopulation Minimum number of patients per database.
#' @param maxDatabasePopulation Maximum number of patients per database.
#' 
#' @details
#' Data capture characteristics are sampled in two stages:
#' 1. Sample the prevalence of each characteristic from a beta distribution with alpha = 1/cpcConsistency
#' and beta = 1/cpcConsistency. Higher cpcConsistency means the prevalences will be closer to 0 and 1.
#' 2. Sample the binary characteristics for each database from a binomial with 1 trial and probability
#' equal to the prevalences sampled in step 1.
#'
#' @returns
#' A settings object to be used in `simulateOne()`.
#'
createSimulationSettings <- function(
    esSettings = EvidenceSynthesis::createSimulationSettings(
      nSites = 5,
      n = 10000,
      treatedFraction = 0.2,
      nStrata = 10,
      minBackgroundHazard = 2e-07,
      maxBackgroundHazard = 2e-05,
      hazardRatio = 2,
      randomEffectSd = 0.25
    ),
    nNegativeControls = 50,
    nOutcomesOfInterest = 10,
    nCaptureProcessChars = 20,
    cpcConsistency = 25,
    bias0 = 0.1,
    biasCpcSd = 0.05
) {
  # Note: currently only supporting per-database n, but no other per-DB parameters
  args <- list()
  for (name in names(formals())) {
    args[[name]] <- get(name)
  }
  return(args)
}

# Simulation functions =----------------------------------------------------------------------------
# settings = createSimulationSettings()
# settings = createSimulationSettings(esSettings = EvidenceSynthesis::createSccsSimulationSettings())
simulateData <- function(seed, settings) {
  set.seed(seed)
  
  simulateOutcome <- function(settings, databaseCpChars, nPatients, outcomeId) {
    
    COHORT_METHOD <- "Cohort method" 
    SCCS <- "SCCS"
    
    isNegativeControl <- outcomeId <= settings$nNegativeControls
    design <- if (is(settings$esSettings, "simulationSettings")) COHORT_METHOD else SCCS
    nDatabases <- settings$esSettings$nSites
    
    biasCpc <- rnorm(settings$nCaptureProcessChars, 0, settings$biasCpcSd) # Bias associated with each data capture process characteristic.
    databaseBias <- settings$bias0 + databaseCpChars %*% biasCpc
    
    # Simulate 1 population at a time so we can add the bias to the effect:
    populations <- vector("list", nDatabases)
    for (i in seq_len(nDatabases)) {
      esSettings <- settings$esSettings
      esSettings$nSites <- 1
      esSettings$n <- if (length(settings$esSettings$n) == 1) settings$esSettings$n else settings$esSettings$n[i]
      if (isNegativeControl) {
        if (design == COHORT_METHOD) {
          esSettings$hazardRatio <- exp(databaseBias[i]) 
        } else {
          esSettings$rateRatio <- exp(databaseBias[i])  
        }
        esSettings$randomEffectSd <- 0
      } else {
        if (design == COHORT_METHOD) {
          esSettings$hazardRatio <- exp(databaseBias[i]) * settings$esSettings$hazardRatio
        } else {
          esSettings$rateRatio <- exp(databaseBias[i]) * settings$esSettings$rateRatio
        }
        esSettings$randomEffectSd <- settings$esSettings$randomEffectSd
      } 
      populations[[i]] <- EvidenceSynthesis::simulatePopulations(esSettings)[[1]]
    }
    
    normalApproximations <- vector("list", nDatabases)
    nonNormalApproximations <- vector("list", nDatabases)
    for (i in seq_len(nDatabases)) {
      if (design == COHORT_METHOD) {
        cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
                                                  data = populations[[i]],
                                                  modelType = "cox")
        treatmentVariable <- "x"
      } else {
        xColnames <- colnames(populations[[i]])[grep("x[0-9]+", colnames(populations[[i]]))]
        formula <- as.formula(paste("y ~ a +", paste0(xColnames, collapse = " + "), " + strata(stratumId) + offset(log(time))"))
        cyclopsData <- Cyclops::createCyclopsData(formula, data = populations[[i]], modelType = "cpr")
        treatmentVariable <- "a"
      }
      cyclopsFit <- suppressWarnings(Cyclops::fitCyclopsModel(cyclopsData))
      normalApproximations[[i]] <- EvidenceSynthesis::approximateLikelihood(
        cyclopsFit = cyclopsFit, 
        parameter = treatmentVariable, 
        approximation = "normal"
      ) |>
        mutate(databaseId = i,
               outcomeId = !!outcomeId)
      
      nonNormalApproximations[[i]] <- EvidenceSynthesis::approximateLikelihood(
        cyclopsFit = cyclopsFit, 
        parameter = treatmentVariable, 
        approximation = "adaptive grid"
      ) |>
        mutate(databaseId = i,
               outcomeId = !!outcomeId)
    }
    normalApproximations <- bind_rows(normalApproximations)
    nonNormalApproximations <- bind_rows(nonNormalApproximations)
    return(list(nonNormalApproximations = nonNormalApproximations,
                normalApproximations = normalApproximations))  
  }
  
  # Data capture characteristics per database (binary):
  cpcPrevalences <- rbeta(settings$nCaptureProcessChars, shape1 = 1 / settings$cpcConsistency, shape2 = 1 / settings$cpcConsistency)
  databaseCpChars = matrix(rbinom(settings$esSettings$nSites * settings$nCaptureProcessChars, 1, cpcPrevalences),
                           byrow = TRUE,
                           nrow = settings$esSettings$nSites ,
                           ncol = settings$nCaptureProcessChars)
  
  nOutcomes <- settings$nNegativeControls + settings$nOutcomesOfInterest
  normalApproximations <- vector("list", nOutcomes)
  nonNormalApproximations <- vector("list", nOutcomes)
  for (outcomeId in seq_len(nOutcomes)) {
    approximations <- simulateOutcome(settings = settings,
                                      databaseCpChars = databaseCpChars,
                                      outcomeId = outcomeId)
    normalApproximations[[outcomeId]] <- approximations$normalApproximations
    nonNormalApproximations[[outcomeId]] <- approximations$nonNormalApproximations
  }
  normalApproximations <- bind_rows(normalApproximations)
  nonNormalApproximations <- bind_rows(nonNormalApproximations)
  return(list(nonNormalApproximations = nonNormalApproximations,
              normalApproximations = normalApproximations))  
}

#' Run one simulation study
#'
#' @param seed           The random seed to use when generating the data.
#' @param settings       An object created with the `createSimulationSettings()` function.
#' @param methodFunction A function that applies the method to the simulated data. See details 
#'                       below.
#' @param ...            Additional parameters to be passed to the `methodFunction`
#'
#' @details
#' The `methodFunction` should have at least two input parameters: 
#' - `data`, the simulated data. This is a list of two dataframes:
#'     - `normalApproximations`, with fields `logRr`, `seLogRr`, `outcomeId`, and `databaseId`.
#'     - `nonNormalApproximations`. This currently contains an adaptive grid. with fields `point`,
#'       `value`, `outcomeId`, and `databaseId`
#' - `settings` The same settings object passed to `simulateOne`.
#' 
#' `databaseId` will be a sequential integer from 1 to `settings$esSettings$nSites`.
#' `outcomeId` is also a sequential integer. 1 ... settings$nNegativeControls are the negative 
#' controls, and settings$nNegativeControls+1 ... settings$nNegativeControls+settings$nOutcomesOfInterest
#' are the outcomes of interest for which estimates need to be produced (see below).
#' 
#' The `methodFunction` is expected to return a data frame with one row per outcome of interest, and
#' the following fields:
#' - `logRr`: Log of the estimated mean effect size.
#' - `logLb`: Log of the lower bound of the 95% confidence or credible interval of the mean effect
#'   size.
#' - `logUb`: Log of the upper bound of the 95% confidence or credible interval of the mean effect
#'   size.
#' - `seLogRr`: Standard error around the (log of the) mean effect size. 
#' - `logPiLb`: Log of the lower bound of the prediction interval.
#' - `logPiUb`: Log of the upper bound of the prediction interval.
#' - `seLogPi`: Standard error of the prediction interval.
#' - `tau`: Random effect: standard deviation around the log of the mean effect size.
#' - `tauLb`: Random effect: Lower bound of the 95% confidence or credible interval around the 
#'            standard deviation around the log of the mean effect size.
#' - `tauUb`: Random effect: Upper bound of the 95% confidence or credible interval around the 
#'            standard deviation around the log of the mean effect size.
#' - `nDatabases`: Number of databasese with valid estimates that contributed to the meta-analytic 
#'            estimate. Can be used when analyzing the results.
#' 
#' @returns
#' Metrics to be used as input for the `evaluateResults()` function.
#' 
simulateOne <- function(seed, settings, methodFunction, ...) {
  cacheFolder <- paste("Simulations/cache", rlang::hash(settings), sep = "_")
  if (!dir.exists(cacheFolder)) {
    dir.create(cacheFolder)
  }
  dummy <- list(...)
  dummy$methodFunction <- methodFunction
  estimatesCacheFilename <- file.path(cacheFolder, sprintf("Estimates_%s_%d.rds", rlang::hash(dummy), seed))
  if (file.exists(estimatesCacheFilename)) {
    estimates <- readRDS(estimatesCacheFilename)
  } else {
    dataCacheFilename <- file.path(cacheFolder, sprintf("Data_%d.rds", seed))
    if (file.exists(dataCacheFilename)) {
      data <- readRDS(dataCacheFilename)
    } else {
      data <- simulateData(seed, settings)
      saveRDS(data, dataCacheFilename)
    }
    args <- list(...)
    args$data <- data
    args$settings <- settings
    estimates <- do.call(methodFunction, args) 
    estimates <- estimates |>
      mutate(seed = seed)
    saveRDS(estimates, estimatesCacheFilename)
  }
  return(estimates)
}

# Code for verifying simulation --------------------------------------------------------------------
# settings <- createSimulationSettings()
plotSystematicErrorDistributions <- function(settings, seed = NULL) {
  # Compute the systematic error distribution within each database and plot it. Used to check if
  # distributions look similar to what is observed in real-world
  
  data <- simulateData(seed, settings)
  
  x <- seq(log(0.1), log(10), length.out = 100)
  compute <- function(x, mcmc) {
    yMcmc <- dnorm(rep(x, nrow(mcmc$chain)), mean = mcmc$chain[, 1], sd = 1/sqrt(mcmc$chain[, 2]))
    return(quantile(yMcmc, c(0.025, 0.5, 0.975)))
  }
  plotData <- list()
  for (i in seq_len(settings$esSettings$nSites)) {
    ncApproximations <- data$normalApproximations |>
      filter(databaseId == i, outcomeId <= settings$nNegativeControls)
    null <- EmpiricalCalibration::fitMcmcNull(
      logRr = ncApproximations$logRr,
      seLogRr = ncApproximations$seLogRr
    )
    ys <- sapply(x, compute, mcmc = attr(null, "mcmc"))
    y <- ys[2, ]
    yMaxLb <- ys[1, ]
    yMaxUb <- ys[3, ]
    normFactor <- max(ys[2, ])
    y <- y / normFactor
    yMaxLb <- yMaxLb / normFactor
    yMaxUb <- yMaxUb / normFactor
    plotData[[i]] <- tibble(
      databaseId = sprintf("Database %s", i),
      x = x,
      yMax = y,
      yMaxLb = yMaxLb,
      yMaxUb =  yMaxUb,
      yMin = 0
    )
    # Are we calibrated (using leave-one-out)?
    # EmpiricalCalibration::plotCalibration(
    #   logRr = ncApproximations$logRr,
    #   seLogRr = ncApproximations$seLogRr,
    #   useMcmc = F
    # )
  }
  plotData <- bind_rows(plotData)
  breaks <- c(0.25, 1, 4, 8)
  ggplot2::ggplot(plotData) +
    ggplot2::geom_vline(xintercept = log(breaks), colour = "#AAAAAA", lty = 1, size = 0.4) +
    ggplot2::geom_vline(xintercept = 0, size = 0.8) +
    ggplot2::geom_ribbon(ggplot2::aes(x = .data$x, ymax = .data$yMax, ymin = .data$yMin), fill = "#FF2700", alpha = 0.6, data = plotData) +
    ggplot2::geom_ribbon(ggplot2::aes(x = .data$x, ymax = .data$yMaxUb, ymin = .data$yMax), fill = "#FF2700", alpha = 0.3, data = plotData) +
    ggplot2::coord_cartesian(xlim = log(c(0.1, 10)), ylim = c(0, 2)) +
    ggplot2::scale_x_continuous("Systematic Error", breaks = log(breaks), labels = breaks) +
    ggplot2::facet_grid(databaseId ~ .) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_text(size = 14),
                   axis.title.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 14),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   strip.text.x = ggplot2::element_text(size = 14),
                   strip.text.y.left = ggplot2::element_text(size = 14, angle = 0, hjust = 0),
                   strip.background = ggplot2::element_blank(),
                   panel.spacing.y = ggplot2::unit(0, "lines"))
}
 

# Baseline approaches ------------------------------------------------------------------------------
# data = simulateData(1, settings)
applyCurrentApproach <- function(data, settings, bayesian = TRUE, approximation = "normal") {
  # Apply the current approach, where we first meta-analyse all outcomes, then use the meta-analytic
  # estimates to fit the empirical null and calibrate.
  
  if (!bayesian && approximation != "normal") {
    stop("Cannot perform non-Bayesian meta-analysis on non-normal approximations")
  }
  
  # group = groups[[2]]
  performSingleMetaAnalysis <- function(group) {
    outcomeId <- group$outcomeId[[1]]
    if (approximation != "normal") {
      group <- group |>
        group_by(databaseId) |>
        group_split()
    }
    if (bayesian) {
      estimate <- suppressMessages(EvidenceSynthesis::computeBayesianMetaAnalysis(group, showProgressBar = FALSE))
      estimate <- tibble(logRr = estimate$mu,
                         seLogRr = estimate$muSe,
                         logLb = estimate$mu95Lb,
                         logUb = estimate$mu95Ub,
                         logPiLb = estimate$predictionInterval95Lb,
                         logPiUb = estimate$predictionInterval95Ub,
                         seLogPi = (estimate$predictionInterval95Ub - estimate$predictionInterval95Lb) / (2 * qnorm(0.975)),
                         nDatabases = sum(!is.na(group$seLogRr)),
                         outcomeId = outcomeId)
    } else {
      meta <- meta::metagen(
        TE = group$logRr, 
        seTE = group$seLogRr,
        sm = "RR", 
        prediction = TRUE,
        control = list(maxiter=1000, stepadj=0.5))
      s <- summary(meta)
      rnd <- s$random
      estimate <- tibble(logRr = rnd$TE,
                         seLogRr = rnd$seTE,
                         logLb = rnd$lower,
                         logUb = rnd$upper,
                         logPiLb = s$predict$lower,
                         logPiUb = s$predict$upper,
                         seLogPi = (s$predict$upper - s$predict$lower) / (2 * qnorm(0.975)),
                         nDatabases = sum(!is.na(group$seLogRr)),
                         outcomeId = outcomeId)
    }
    return(estimate)
  }
  
  if (approximation == "normal") {
    groups <- data$normalApproximations |>
      group_by(outcomeId) |>
      group_split()
  } else {
    groups <- data$nonNormalApproximations |>
      group_by(outcomeId) |>
      group_split()
  }
  
  estimates <- lapply(groups, performSingleMetaAnalysis)
  estimates <- bind_rows(estimates)
  
  ncEstimates <- estimates |>
    filter(outcomeId <= settings$nNegativeControls)
  ooiEstimates <- estimates |>
    filter(outcomeId > settings$nNegativeControls)
  
  null <- EmpiricalCalibration::fitMcmcNull(
    logRr = ncEstimates$logRr,
    seLogRr = ncEstimates$seLogRr
  )
  calibratedCi <- EmpiricalCalibration::calibrateConfidenceInterval(
    logRr = ooiEstimates$logRr,
    seLogRr = ooiEstimates$seLogRr,
    model = EmpiricalCalibration::convertNullToErrorModel(null)
  ) |>
    transmute(logRr = logRr, logLb = logLb95Rr, logUb = logUb95Rr, seLogRr = seLogRr)
  nullPi <- EmpiricalCalibration::fitMcmcNull(
    logRr = (ncEstimates$logPiLb + ncEstimates$logPiUb) / 2,
    seLogRr = ncEstimates$seLogPi
  )
  calibratedPi <- EmpiricalCalibration::calibrateConfidenceInterval(
    logRr = (ooiEstimates$logPiLb + ooiEstimates$logPiUb) / 2,
    seLogRr = ooiEstimates$seLogPi,
    model = EmpiricalCalibration::convertNullToErrorModel(nullPi)
  ) |>
    transmute(logPiLb = logLb95Rr, logPiUb = logUb95Rr, seLogPi = seLogRr)
  estimatesOois <- bind_cols(calibratedCi, calibratedPi) |>
    mutate(tau = as.numeric(NA), tauLb = as.numeric(NA), tauUb = as.numeric(NA)) |>
    mutate(nDatabases = ooiEstimates$nDatabases)
  return(estimatesOois)
}

applyNaiveApproach <- function(data, settings, bayesian = TRUE) {
  # Apply the naive approach: calibrate estimates within each database, then meta-analyse calibrated
  # estimates. This should in theory be bad because it does not take into account that systematic
  # error is correlated between databases.
  
  calibratedEstimates <- vector("list", settings$esSettings$nSites)
  for (i in seq_len(settings$esSettings$nSites)) {
    ncApproximations <- data$normalApproximations |>
      filter(outcomeId <= settings$nNegativeControls,
             databaseId == i)
    ooiApproximations <- data$normalApproximations |>
      filter(outcomeId > settings$nNegativeControls,
             databaseId == i)
    if (nrow(ooiApproximations) == 0) {
      next
    }
    null <- suppressWarnings(EmpiricalCalibration::fitMcmcNull(
      logRr = ncApproximations$logRr,
      seLogRr = ncApproximations$seLogRr
    ))
    calibratedEstimates[[i]] <- EmpiricalCalibration::calibrateConfidenceInterval(
      logRr = ooiApproximations$logRr,
      seLogRr = ooiApproximations$seLogRr,
      model = EmpiricalCalibration::convertNullToErrorModel(null)
    ) |>
      mutate(databaseId = i,
             outcomeId = ooiApproximations$outcomeId)
  }
  calibratedEstimates <- bind_rows(calibratedEstimates)
  groups <- calibratedEstimates |>
    group_by(outcomeId) |>
    group_split()
  estimatesOois <- vector("list", settings$nOutcomesOfInterest)
  # group = groups[[1]]
  for (group in groups) {
    if (bayesian) {
      estimate <- suppressMessages(EvidenceSynthesis::computeBayesianMetaAnalysis(group, showProgressBar = FALSE))
      estimatesOois[[group$outcomeId[1]]] <- tibble(logRr = estimate$mu,
                                                    seLogRr = estimate$muSe,
                                                    logLb = estimate$mu95Lb,
                                                    logUb = estimate$mu95Ub,
                                                    logPiLb = estimate$predictionInterval95Lb,
                                                    logPiUb = estimate$predictionInterval95Ub,
                                                    seLogPi = (estimate$predictionInterval95Ub - estimate$predictionInterval95Lb) / (2 * qnorm(0.975)),
                                                    tau = estimate$tau,
                                                    tauLb = estimate$tau95Lb,
                                                    tauUb = estimate$tau95Ub,
                                                    nDatabases = sum(!is.na(group$seLogRr)))
    } else {
      meta <- meta::metagen(
        TE = group$logRr, 
        seTE = group$seLogRr, 
        sm = "RR", 
        prediction = TRUE,
        control = list(maxiter=1000, stepadj=0.5))
      s <- summary(meta)
      rnd <- s$random
      estimatesOois[[group$outcomeId[1]]] <- tibble(logRr = rnd$TE,
                                                    seLogRr = rnd$seTE,
                                                    logLb = rnd$lower,
                                                    logUb = rnd$upper,
                                                    logPiLb = s$predict$lower,
                                                    logPiUb = s$predict$upper,
                                                    seLogPi = s$predict$seTE,
                                                    tau = s$tau,
                                                    tauLb = s$lower.tau,
                                                    tauUb = s$upper.tau,
                                                    nDatabases = sum(!is.na(group$seLogRr)))
    }
  }
  estimatesOois <- bind_rows(estimatesOois)
  return(estimatesOois)
}

# Evaluation function ------------------------------------------------------------------------------
computeMeanPrecision <- function(lb, ub) {
  se <- (ub-lb) / (2 * qnorm(0.975))
  se <- pmax(0.001, se)
  precision <- 1/se^2
  return(exp(mean(log(precision))))
}

evaluateResults <- function(results, settings) {
  results <- bind_rows(results)
  
  if (is(settings$esSettings, "simulationSettings")) {
    trueLogRr <- log(settings$esSettings$hazardRatio)
  } else {
    trueLogRr <- log(settings$esSettings$rateRatio)
  }
  trueTau <- settings$esSettings$randomEffectSd
  # Deliberately computing precision based on width of intervals and not per-method SE estimates for 
  # consistency:
  metrics <- results |>
    summarise(coverageCi = mean(logLb <= trueLogRr & logUb >= trueLogRr),
              precisionCi = computeMeanPrecision(logLb, logUb),
              coveragePi = mean(pnorm(logPiUb, trueLogRr, trueTau) - pnorm(logPiLb, trueLogRr, trueTau)),
              precisionPi = computeMeanPrecision(logPiLb, logPiUb),
              coverageTau = mean(tauLb <= trueTau & tauUb >= trueTau),
              precisionTau = computeMeanPrecision(tauLb, tauUb))
  writeLines(paste(paste(names(metrics), format(as.vector(metrics), digits = 3), sep = " = "), collapse = ", "))
  invisible(metrics)
}
