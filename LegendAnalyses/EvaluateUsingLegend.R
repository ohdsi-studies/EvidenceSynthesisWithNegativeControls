# Assumes you've run LegendAnalyses/DownloadLegendNcEstimates.R

# Need to have these packages installed:
# install.packages("mvtnorm")
# install.packages("EmpiricalCalibration")
# install.packages("EvidenceSynthesis")
# install.packages("ParallelLogger")

cluster <- ParallelLogger::makeCluster(10)
ParallelLogger::clusterRequire(cluster, "dplyr")
ParallelLogger::clusterRequire(cluster, "survival")
ParallelLogger::clusterRequire(cluster, "tidyr")

source("Simulations//SimulationFunctions.R")
snow::clusterExport(cluster, c("applyNaiveApproach", "applyCurrentApproach"))

source("ModelCode/SystematicErrorModel.R")
snow::clusterExport(cluster, c("fitSystematicErrorModel", "calibrateCiRandomEffects", "calibrateCiBayesian", "computePredictionInterval"))

source("LegendAnalyses/EvaluationFunctions.R")


results <- estimateLeaveOneOut(cluster, methodFunction = applyNaiveApproach, bayesian = TRUE)
evaluateLegendResults(results)
# Overall PS matched:
#   coverageCi = 0.989, precisionCi = 23.7, coveragePi = 1, precisionPi = 7.57, coverageTau = 0, precisionTau = 53
# 
# Overall unadjusted:
#   coverageCi = 0.896, precisionCi = 21.5, coveragePi = 0.967, precisionPi = 7.23, coverageTau = 0, precisionTau = 50.7


results <- estimateLeaveOneOut(cluster, methodFunction = applyCurrentApproach, bayesian = TRUE)
evaluateLegendResults(results)
# Overall PS matched:
#   coverageCi = 0.989, precisionCi = 24.1, coveragePi = 1, precisionPi = 7.59, coverageTau = NA, precisionTau = NA
# 
# Overall unadjusted:
#   coverageCi = 0.914, precisionCi = 13.8, coveragePi = 0.967, precisionPi = 8.36, coverageTau = NA, precisionTau = NA


results <- estimateLeaveOneOut(cluster, methodFunction = applyGeneralizedModel, bayesian = TRUE)
evaluateLegendResults(results)
# Overall PS matched:
#   coverageCi = 0.992, precisionCi = 21.8, coveragePi = 1, precisionPi = 7.65, coverageTau = 0, precisionTau = 53.9
# 
# Overall unadjusted:
#   coverageCi = 0.949, precisionCi = 10.6, coveragePi = 0.982, precisionPi = 5.86, coverageTau = 0, precisionTau = 62.4


ParallelLogger::stopCluster(cluster)
