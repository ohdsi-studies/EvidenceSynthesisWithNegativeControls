# Simulate systematic error in a network of databases. 

# Need to have these packages installed:
# install.packages("mvtnorm")
# install.packages("EmpiricalCalibration")
# install.packages("EvidenceSynthesis")
# install.packages("ParallelLogger")

cluster <- ParallelLogger::makeCluster(10)
ParallelLogger::clusterRequire(cluster, "dplyr")
ParallelLogger::clusterRequire(cluster, "survival")

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
# coverageCi = 0.855, precisionCi = 17.6, coveragePi = 0.915, precisionPi = 10.7, coverageTau = 1, precisionTau = 23.4

results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings, methodFunction = applyCurrentApproach, approximation = "normal", bayesian = FALSE)
evaluateResults(results, settings)
# coverageCi = 0.874, precisionCi = 16.7, coveragePi = 0.895, precisionPi = 13.6, coverageTau = NA, precisionTau = NA

results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings, methodFunction = applyCurrentApproach, approximation = "non-normal", bayesian = TRUE)
evaluateResults(results, settings)








ParallelLogger::stopCluster(cluster, settings)
