library(dplyr)

computePredictionInterval <- function(estimate, alpha = 0.05) {
  qmixnorm <- function(p, means, sds) {
    pmix <- function(x) {
      mean(pnorm(x, mean = means, sd = sds))
    }
    sapply(p, function(pVal) {
      if (pVal <= 0) return(-Inf)
      if (pVal >= 1) return(Inf)
      objective_function <- function(x) {
        pmix(x) - pVal
      }
      search_interval <- c(min(means) - 10 * max(sds), max(means) + 10 * max(sds))
      uniroot(objective_function, interval = search_interval)$root
    })
  }

  traces <- attr(estimate, "traces")
  predictionInterval <- HDInterval::hdi(qmixnorm, credMass = 1 - alpha, means = traces[, 1], sds = traces[, 2])
  return(predictionInterval)

  # To verify: use very large sample:
  # predictionInterval
  # predictions <- do.call(c, lapply(seq_len(nrow(traces)), function(i) rnorm(100000, traces[i, 1], traces[i, 2])))
  # predictionInterval <- HDInterval::hdi(predictions, credMass = 1 - alpha)
  # predictionInterval
  # -0.06760276  2.06989313
}

# tau = 0.1; seLogRrs = c(0.1, 0.2, 0.05, 0.15)
computeMetaAnalysisMdrr <- function(tau, seLogRrs, maxCores = 10, power = 0.8, alpha = 0.05) {

  sampleOne <- function(i, tau, seLogRrs, alpha) {
    trueLogRrs <- rnorm(length(seLogRrs), 0, tau)
    logRrs <- rnorm(length(seLogRrs), trueLogRrs, seLogRrs)
    data <- tibble(
      logRr = logRrs,
      seLogRr = seLogRrs
    )
    maEstimate <- suppressMessages(EvidenceSynthesis::computeBayesianMetaAnalysis(data, showProgressBar = FALSE, alpha = alpha))
    pi <- computePredictionInterval(maEstimate, alpha = alpha)
    maEstimate <- maEstimate |>
      mutate(piLb = pi[1],
             piUb = pi[2])
    return(maEstimate)
  }

  cluster <- ParallelLogger::makeCluster(maxCores)
  ParallelLogger::clusterRequire(cluster, "dplyr")
  parallel::clusterExport(cluster, "computePredictionInterval")
  sampleSize <- 1000
  maEstimates <- ParallelLogger::clusterApply(cluster, seq_len(sampleSize), sampleOne, tau = tau, seLogRrs = seLogRrs, alpha = alpha)
  ParallelLogger::stopCluster(cluster)
  maEstimates <- bind_rows(maEstimates)
  mdrr <- tibble(
    mdrrPi = exp(-quantile(maEstimates$piLb, 1 - power)),
    mdrrCi = exp(-quantile(maEstimates$mu95Lb, 1 - power))
  )
  return(mdrr)
}
# system.time(mdrr <- computeMetaAnalysisMdrr(tau, seLogRrs))
# user  system elapsed
# 1.486   1.557 163.295
# mdrr
# mdrrPi mdrrCi
#  1.94   1.40
#  1.96   1.40
