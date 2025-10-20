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
    biasCpcSd = 0.1
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
                         seLogPi = s$predict$seTE,
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
    mutate(tau = as.numeric(NA), tauLb = as.numeric(NA), tauUb = as.numeric(NA))
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
                                                    tauUb = estimate$tau95Ub)
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
                                                    tauUb = s$upper.tau)
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
