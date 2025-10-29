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



results <- estimateLeaveOneOut(cluster, methodFunction = applyCurrentApproach, bayesian = TRUE)
evaluateLegendResults(results)


results <- estimateLeaveOneOut(cluster, methodFunction = applyGeneralizedModel, bayesian = TRUE)
evaluateLegendResults(results)



ParallelLogger::stopCluster(cluster)
