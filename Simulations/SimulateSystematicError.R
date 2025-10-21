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
    n = 10000,
    minBackgroundHazard = 2e-05,
    maxBackgroundHazard = 2e-04,
    hazardRatio = 2,
    randomEffectSd = 0.25
  ),
  nNegativeControls = 100,
  nOutcomesOfInterest = 10,
  cpcConsistency = 0.1
)
plotSystematicErrorDistributions(settings)

results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings, methodFunction = applyNaiveApproach, bayesian = FALSE)
evaluateResults(results, settings)

results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings, methodFunction = applyCurrentApproach, approximation = "normal", bayesian = FALSE)
evaluateResults(results, settings)

results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings, methodFunction = applyCurrentApproach, approximation = "normal", bayesian = TRUE)
evaluateResults(results, settings)

results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings, methodFunction = applyGeneralizedModel, bayesian = FALSE)
evaluateResults(results, settings)

results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings, methodFunction = applyGeneralizedModel, bayesian = TRUE)
evaluateResults(results, settings)

ParallelLogger::stopCluster(cluster)
