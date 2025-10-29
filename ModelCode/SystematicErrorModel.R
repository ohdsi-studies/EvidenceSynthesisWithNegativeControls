# Use multivariate normal distribution to generalize old systematic error model to multiple
# databases.

source("ModelCode/PredictionInterval.R")

library(mvtnorm)
library(dplyr)
library(tidyr)

# Fitting the systematic errmor model ------------------------------------------

fitSystematicErrorModel <- function(data) {
  
  # Non-Bayesian model fitting. Outputs MLE of mean vector and covariance matrix. Currently using
  # full-rank covariance matrix.

  constructCovarianceMatrix <- function(params, nDatabases) {
    choleskyParams <- params[(nDatabases + 1):length(params)]
    choleskyMatrix <- matrix(0, nrow = nDatabases, ncol = nDatabases)
    choleskyMatrix[lower.tri(choleskyMatrix, diag = TRUE)] <- choleskyParams
    covarianceMatrix <- choleskyMatrix %*% t(choleskyMatrix)
    return(covarianceMatrix)
  }

  computeNegativeLogLikelihood <- function(params, logRrMatrix, seLogRrMatrix, nDatabases) {
    meanVector <- params[1:nDatabases]
    covarianceMatrix <- constructCovarianceMatrix(params, nDatabases)
    totalLogLikelihood <- 0
    for (i in 1:nrow(logRrMatrix)) {
      logRr <- logRrMatrix[i, ]
      seLogRr <- seLogRrMatrix[i, ]

      # Check for missing data for this outcome
      validIndices <- !is.na(seLogRr)

      if (sum(validIndices) < 2) {
        next 
      }
      # Subset matrices to non-missing data
      logRrSubset <- logRr[validIndices]
      seLogRrSubset <- seLogRr[validIndices]
      meanVectorSubset <- meanVector[validIndices]
      covarianceMatrixSubset <- covarianceMatrix[validIndices, validIndices]

      # Within-database variance matrix (measurement error)
      varianceMatrixSubset <- diag(seLogRrSubset^2)

      # Total covariance for this outcome
      totalCovarianceSubset <- covarianceMatrixSubset + varianceMatrixSubset

      logLikelihood <- tryCatch({
        mvtnorm::dmvnorm(logRrSubset, mean = meanVectorSubset, sigma = totalCovarianceSubset, log = TRUE)
      }, error = function(e) {
        return(-1e10)
      })
      totalLogLikelihood <- totalLogLikelihood + logLikelihood
    }
    return(-totalLogLikelihood)
  }

  # Get database IDs and count
  databaseIds <- unique(data$databaseId)
  nDatabases <- length(databaseIds)

  dataWide <- data %>%
    pivot_wider(
      id_cols = outcomeId,
      names_from = databaseId,
      values_from = c(logRr, seLogRr),
      names_sort = TRUE
    )

  databaseIds <- sort(databaseIds)
  logRrCols <- paste0("logRr_", databaseIds)
  seLogRrCols <- paste0("seLogRr_", databaseIds)
  logRrMatrix <- as.matrix(dataWide[, logRrCols])
  seLogRrMatrix <- as.matrix(dataWide[, seLogRrCols])

  initialMean <- rep(0, nDatabases)
  initialCholesky <- t(chol(diag(nDatabases) * 0.1))
  initialCholeskyParams <- initialCholesky[lower.tri(initialCholesky, diag = TRUE)]
  initialParams <- c(initialMean, initialCholeskyParams)
  optimResult <- optim(
    par = initialParams,
    fn = computeNegativeLogLikelihood,
    logRrMatrix = logRrMatrix,
    seLogRrMatrix = seLogRrMatrix,
    nDatabases = nDatabases,
    method = "L-BFGS-B",
    control = list(maxit = 2000)
  )
  if (optimResult$convergence != 0) {
    warning(sprintf("Optim failed to converge with covergence = %s", optimResult$convergence))
  }
  finalParams <- optimResult$par
  finalMean <- finalParams[1:nDatabases]
  names(finalMean) <- databaseIds

  finalCovarianceMatrix <- constructCovarianceMatrix(finalParams, nDatabases)
  rownames(finalCovarianceMatrix) <- databaseIds
  colnames(finalCovarianceMatrix) <- databaseIds
  return(list(
    mean = finalMean,
    covarianceMatrix = finalCovarianceMatrix
  ))
}

# Calibration functions --------------------------------------------------------
calibrateCiRandomEffects <- function(model, newData) {
  modelDbIds <- names(model$mean)
  dataDbIds <- newData$databaseId[!is.na(newData$seLogRr)]
  availableDbIds <- intersect(modelDbIds, dataDbIds)

  if (length(availableDbIds) == 0) {
    stop("None of the databases in newData are present in the systematic error model.")
  }

  # Subset the model parameters and new data to the common set of databases
  meanSys <- model$mean[availableDbIds]
  covSys <- model$covarianceMatrix[availableDbIds, availableDbIds, drop = FALSE]
  newData <- newData[newData$databaseId %in% availableDbIds, ]
  newData <- newData[match(availableDbIds, newData$databaseId), ] # Ensure order

  nDatabases <- length(availableDbIds)

  if (nDatabases == 1) {
    muHatFinal <- newData$logRr - meanSys
    muVarFinal <- covSys[1, 1] + newData$seLogRr^2
    ciLowerLog <- muHatFinal - 1.96 * sqrt(muVarFinal)
    ciUpperLog <- muHatFinal + 1.96 * sqrt(muVarFinal)
    return(data.frame(
      estimate = exp(as.numeric(muHatFinal)),
      ciLower = exp(as.numeric(ciLowerLog)),
      ciUpper = exp(as.numeric(ciUpperLog)),
      tau2 = NA,
      nDatabases = nDatabases,
      logRr = as.numeric(muHatFinal),
      seLogRr = sqrt(muVarFinal)
    ))

  } else {
    yPrime <- newData$logRr - meanSys
    sigmaTotal <- covSys + diag(newData$seLogRr^2)
    ones <- rep(1, nDatabases)
    negLogLike <- function(logTau2) {
      tau2 <- exp(logTau2)
      if (is.infinite(tau2)) return(1e10)
      V <- sigmaTotal + diag(tau2, nDatabases)
      VInv <- try(solve(V), silent = TRUE)
      if (inherits(VInv, "try-error")) return(1e10)
      muHat <- (t(ones) %*% VInv %*% yPrime) / (t(ones) %*% VInv %*% ones)
      ll <- mvtnorm::dmvnorm(yPrime, mean = rep(as.numeric(muHat), nDatabases), sigma = V, log = TRUE)
      return(-ll)
    }
    opt <- optim(par = log(0.01), fn = negLogLike, method = "BFGS", hessian = TRUE)
    estimatedTau2 <- exp(opt$par)
    fisherInfo <- opt$hessian[1, 1]
    seLogTau2 <- NA
    tau2CiLower <- NA
    tau2CiUpper <- NA
    if (fisherInfo > 0) { 
      seLogTau2 <- sqrt(1 / fisherInfo)
      ciLowerLogTau2 <- opt$par - 1.96 * seLogTau2
      ciUpperLogTau2 <- opt$par + 1.96 * seLogTau2
      tau2CiLower <- exp(ciLowerLogTau2)
      tau2CiUpper <- exp(ciUpperLogTau2)
    }
    VFinal <- sigmaTotal + diag(estimatedTau2, nDatabases)
    VFinalInv <- solve(VFinal)
    muHatFinal <- (t(ones) %*% VFinalInv %*% yPrime) / (t(ones) %*% VFinalInv %*% ones)
    muVarFinal <- 1 / (t(ones) %*% VFinalInv %*% ones)
    ciLowerLog <- muHatFinal - 1.96 * sqrt(muVarFinal)
    ciUpperLog <- muHatFinal + 1.96 * sqrt(muVarFinal)

    return(data.frame(
      logRr = as.numeric(muHatFinal),
      seLogRr = sqrt(muVarFinal),
      logLb = as.numeric(ciLowerLog),
      logUb = as.numeric(ciUpperLog),
      tau2 = estimatedTau2,
      seLogTau2 = seLogTau2,
      tau2CiLower = tau2CiLower,
      tau2CiUpper = tau2CiUpper
    ))
  }
}


calibrateCiBayesian <- function(model, newData) {

  logLikelihood <- function(params, yPrime, sigmaTotal) {
    mu <- params[1]
    tau <- params[2]
    if (tau < 0)
      return(-Inf)
    nDatabases <- length(yPrime)
    V <- sigmaTotal + diag(tau ^ 2, nDatabases)
    meanVec <- rep(mu, nDatabases)
    ll <- tryCatch({
      mvtnorm::dmvnorm(yPrime, mean = meanVec, sigma = V, log = TRUE)
    }, error = function(e) {
      return(-1e10)
    })
    return(ll)
  }

  logPrior <- function(params) {
    mu <- params[1]
    tau <- params[2]
    muPrior <- dnorm(mu, mean = 0, sd = 2, log = TRUE)
    tauPrior <- dnorm(tau, mean = 0, sd = 0.5, log = TRUE)
    return(muPrior + tauPrior)
  }

  minLogPosterior <- function(params, yPrime, sigmaTotal) {
    if (any(is.nan(params)) || any(is.infinite(params))) {
      return(-Inf)
    }

    ll <- logLikelihood(params, yPrime, sigmaTotal)
    lp <- logPrior(params)

    return(- ll - lp)
  }

  proposalFunction <- function(param, scale) {
    dim <- length(param)
    draw <- rnorm(dim, mean = param, sd = scale)

    # Tau cannot be negative:
    draw[2] <- abs(draw[2])
    return(draw)
  }

  runMetropolisMcmc <- function(startValue, iterations, scale, yPrime, sigmaTotal) {
    dim <- length(startValue)
    chain <- array(dim = c(iterations + 1, dim))
    logLik <- array(dim = c(iterations + 1, 1))
    acc <- array(dim = c(iterations + 1, 1))

    logLik[1] <- -minLogPosterior(startValue, yPrime, sigmaTotal)
    chain[1, ] <- c(startValue)
    acc[1] <- 1

    for (i in 1:iterations) {
      # print(paste('itr =', i))
      proposal <- proposalFunction(chain[i, ], scale = scale)
      newLogLik <- tryCatch(-minLogPosterior(proposal, yPrime, sigmaTotal), error = function(e) {
        -1e+10
      })

      # print(paste(paste(proposal, collapse = ","), newLogLik))
      prob <- exp(newLogLik - logLik[i])
      if (runif(1) < prob) {
        chain[i + 1, ] <- proposal
        logLik[i + 1] <- newLogLik
        acc[i + 1] <- 1
      } else {
        chain[i + 1, ] <- chain[i, ]
        logLik[i + 1] <- logLik[i]
        acc[i + 1] <- 0
      }
    }
    result <- list(logLik = logLik, chain = chain, acc = acc)
    return(result)
  }

  binarySearchMu <- function(modeMu,
                             modeTau,
                             alpha = 0.1,
                             yPrime,
                             sigmaTotal,
                             precision = 1e-07) {
    q <- qchisq(1 - alpha, 1) / 2
    L <- modeMu
    H <- 10
    llMode <- -minLogPosterior(c(modeMu, modeTau), yPrime = yPrime, sigmaTotal = sigmaTotal)
    while (H >= L) {
      M <- L + (H - L) / 2
      llM <- -minLogPosterior(c(M, modeTau), yPrime = yPrime, sigmaTotal = sigmaTotal)
      metric <- llMode - llM - q
      # writeLines(paste('M =', M, 'Metric = ',metric))
      if (metric > precision) {
        H <- M
      } else if (-metric > precision) {
        L <- M
      } else {
        return(abs(M - modeMu))
      }
      if (M == modeMu || M == 10) {
        return(0)
      }
    }
  }

  binarySearchTau <- function(modeMu,
                              modeTau,
                              alpha = 0.1,
                              yPrime,
                              sigmaTotal,
                              precision = 1e-07) {
    q <- qchisq(1 - alpha, 1) / 2
    llMode <- -minLogPosterior(c(modeMu, modeTau), yPrime = yPrime, sigmaTotal = sigmaTotal)
    L <- modeTau
    for (i in 1:10) {
      H <- modeTau + exp(i)
      llM <- -minLogPosterior(c(modeMu, H), yPrime = yPrime, sigmaTotal = sigmaTotal)
      metric <- llMode - llM - q
      if (metric > 0) {
        break
      }
    }
    if (i == 10) {
      return(0)
    }
    while (H >= L) {
      M <- L + (H - L) / 2
      llM <- -minLogPosterior(c(modeMu, M), yPrime = yPrime, sigmaTotal = sigmaTotal)
      metric <- llMode - llM - q
      # writeLines(paste('M =', M, 'Metric = ',metric))
      if (metric > precision) {
        H <- M
      } else if (-metric > precision) {
        L <- M
      } else {
        return(abs(M - modeTau))
      }
      if (M == modeTau) {
        return(0)
      }
    }
  }

  modelDbIds <- names(model$mean)
  newData <- newData |>
    filter(!is.na(seLogRr))
  availableDbIds <- intersect(modelDbIds, newData$databaseId)

  if (length(availableDbIds) == 0) {
    stop("None of the databases in newData are present in the systematic error model.")
  }

  meanSys <- model$mean[availableDbIds]
  covSys <- model$covarianceMatrix[availableDbIds, availableDbIds, drop = FALSE]
  newData <- newData[newData$databaseId %in% availableDbIds, ]
  newData <- newData[match(availableDbIds, newData$databaseId), ]

  nDatabases <- length(availableDbIds)

  # if (nDatabases == 1) {
  #   warning("Only one database estimate is available. Bayesian analysis for tau is not possible. Returning frequentist result.")
  #   null <- c(meanSys, sqrt(covSys))
  #   class(null) <- "null"
  #   calibratedCi <- EmpiricalCalibration::calibrateConfidenceInterval(
  #     logRr = newData$logRr,
  #     seLogRr = newData$seLogRr,
  #     model = EmpiricalCalibration::convertNullToErrorModel(null)
  #   )
  #   return(data.frame(
  #     mu = calibratedCi$logRr,
  #     muLb = calibratedCi$logLb95Rr,
  #     muUb = calibratedCi$logUb95Rr,
  #     tau = NA,
  #     tauLb = NA,
  #     tauUb = NA,
  #     piLb = NA,
  #     piUb = NA
  #   ))
  # }

  yPrime <- newData$logRr - meanSys
  if (nDatabases == 1) {
    sigmaTotal <- covSys + newData$seLogRr^2
  } else {
    sigmaTotal <- covSys + diag(newData$seLogRr^2)
  }
  fit <- optim(c(0, 1), minLogPosterior, yPrime = yPrime, sigmaTotal = sigmaTotal)

  # Profile likelihood for roughly correct scale:
  scale <- binarySearchMu(modeMu = fit$par[1],
                          modeTau = fit$par[2],
                          yPrime = yPrime,
                          sigmaTotal = sigmaTotal)
  scale <- c(scale, binarySearchTau(fit$par[1],
                                    fit$par[2],
                                    yPrime = yPrime,
                                    sigmaTotal = sigmaTotal))
  mcmc <- runMetropolisMcmc(fit$par, iterations = 100000, scale, yPrime = yPrime, sigmaTotal = sigmaTotal)
  mean(mcmc$acc)
  traces <- mcmc$chain
  mu = mean(traces[, 1])
  hdiMu <- HDInterval::hdi(traces[, 1], credMass = 0.95)
  hdiTau <- HDInterval::hdi(traces[, 2], credMass = 0.95)

  estimate <- list()
  attr(estimate, "traces") <- traces
  pi <- computePredictionInterval(estimate)

  return(data.frame(
    mu = mu,
    muLb = hdiMu[1],
    muUb = hdiMu[2],
    tau = median(traces[, 2]),
    tauLb = hdiTau[1],
    tauUb = hdiTau[2],
    piLb = pi[1],
    piUb = pi[2],
    seLogRr = sqrt(mean((traces[, 1] - mu)^2)),
    nDatabases = nDatabases
  ))
}


# Applying the systematic error model ------------------------------------------
applyGeneralizedModel <- function(data, settings, bayesian = TRUE) {
  # Apply Martijn's generalized calibration model. Currently only supports non-Bayesian approach for
  # the error model, but the meta-analysis is Bayesian.
  
  ncApproximations <- data$normalApproximations |>
    filter(outcomeId <= settings$nNegativeControls)
  ooiGroups <- data$normalApproximations |>
    filter(outcomeId > settings$nNegativeControls) |>
    group_by(outcomeId) |>
    group_split()
  
  model <- fitSystematicErrorModel(ncApproximations)
  
  estimates <- vector("list", settings$nOutcomesOfInterest)
  for (i in seq_len(settings$nOutcomesOfInterest)) {
    newData <- ooiGroups[[i]] |>
      filter(!is.na(seLogRr))
    if (bayesian) {
      estimate <- calibrateCiBayesian(model, newData)
      estimates[[i]] <- tibble(logRr = estimate$mu,
                              seLogRr = estimate$seLogRr,
                              logLb = estimate$muLb,
                              logUb = estimate$muUb,
                              logPiLb = estimate$piLb,
                              logPiUb = estimate$piUb,
                              seLogPi = (estimate$piUb - estimate$piLb) / (2 * qnorm(0.975)),
                              tau = estimate$tau,
                              tauLb = estimate$tauLb,
                              tauUb = estimate$tauUb,
                              outcomeId = newData$outcomeId[1])

        
    } else {
      estimate <- calibrateCiRandomEffects(model, newData)
      estimates[[i]] <- tibble(logRr = estimate$logRr,
                               seLogRr = estimate$seLogRr,
                               logLb = estimate$logLb,
                               logUb = estimate$logUb,
                               logPiLb = as.numeric(NA),
                               logPiUb = as.numeric(NA),
                               seLogPi = as.numeric(NA),
                               tau = sqrt(estimate$tau2),
                               tauLb = sqrt(estimate$tau2CiLower),
                               tauUb = sqrt(estimate$tau2CiUpper),
                               outcomeId = newData$outcomeId[1])
    }
  }
  estimates <- bind_rows(estimates)
  return(estimates)
}

