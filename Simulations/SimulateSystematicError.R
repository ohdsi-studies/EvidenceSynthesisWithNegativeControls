# Simulate systematic error in a network of databases. 

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
snow::clusterExport(cluster, c("simulateOne", "simulateData", "applyNaiveApproach", "applyCurrentApproach"))

source("ModelCode/SystematicErrorModel.R")
snow::clusterExport(cluster, c("fitSystematicErrorModel", "calibrateCiRandomEffects", "calibrateCiBayesian", "computePredictionInterval"))

# Simulation with lots of data, low consistency ----------------------------------------------------
settings <- createSimulationSettings(
  esSettings = EvidenceSynthesis::createSimulationSettings(
    nSites = 10,
    n = 10000
  ),
  nNegativeControls = 100,
  nOutcomesOfInterest = 10,
  cpcConsistency = 0.1
)

results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings, methodFunction = applyNaiveApproach, bayesian = FALSE)
evaluateResults(results, settings)
# coverageCi = 0.848, precisionCi = 17.5, coveragePi = 0.914, precisionPi = 10.6, coverageTau = 0.999, precisionTau = 24

results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings, methodFunction = applyCurrentApproach, approximation = "normal", bayesian = FALSE)
evaluateResults(results, settings)
# coverageCi = 0.874, precisionCi = 16.7, coveragePi = 0.89, precisionPi = 13.9, coverageTau = NA, precisionTau = NA

results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings, methodFunction = applyCurrentApproach, approximation = "normal", bayesian = TRUE)
evaluateResults(results, settings)
# coverageCi = 0.902, precisionCi = 14.1, coveragePi = 0.997, precisionPi = 4.63, coverageTau = NA, precisionTau = NA

results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings, methodFunction = applyGeneralizedModel, bayesian = FALSE)
evaluateResults(results, settings)
# coverageCi = 0.888, precisionCi = 15.5, coveragePi = NA, precisionPi = NA, coverageTau = NA, precisionTau = NA


results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings, methodFunction = applyGeneralizedModel, bayesian = TRUE)
evaluateResults(results, settings)
# coverageCi = 0.908, precisionCi = 13.1, coveragePi = 0.995, precisionPi = 4.84, coverageTau = 0, precisionTau = 35.3

ParallelLogger::stopCluster(cluster)
