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
# coverageCi = 0.812, precisionCi = 93.7, coveragePi = 0.845, precisionPi = 16.7, coverageTau = 0.947, precisionTau = 65.1

results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings, methodFunction = applyCurrentApproach, approximation = "normal", bayesian = FALSE)
evaluateResults(results, settings)
# coverageCi = 0.945, precisionCi = 44.2, coveragePi = 0.915, precisionPi = 11, coverageTau = NA, precisionTau = NA

results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings, methodFunction = applyCurrentApproach, approximation = "normal", bayesian = TRUE)
evaluateResults(results, settings)
# coverageCi = 0.945, precisionCi = 44.3, coveragePi = 0.944, precisionPi = 9.64, coverageTau = NA, precisionTau = NA

results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings, methodFunction = applyGeneralizedModel, bayesian = FALSE)
evaluateResults(results, settings)
# coverageCi = 0.934, precisionCi = 48.7, coveragePi = NA, precisionPi = NA, coverageTau = 0.986, precisionTau = 0

results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings, methodFunction = applyGeneralizedModel, bayesian = TRUE)
evaluateResults(results, settings)
# coverageCi = 0.95, precisionCi = 41.5, coveragePi = 0.946, precisionPi = 9.7, coverageTau = 0.968, precisionTau = 83.9


# Simulation with lots of data, high consistency ----------------------------------------------------
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
  cpcConsistency = 25
)
plotSystematicErrorDistributions(settings)

results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings, methodFunction = applyNaiveApproach, bayesian = FALSE)
evaluateResults(results, settings)
# coverageCi = 0.709, precisionCi = 100, coveragePi = 0.798, precisionPi = 20.1, coverageTau = 0.956, precisionTau = 71.9

results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings, methodFunction = applyCurrentApproach, approximation = "normal", bayesian = FALSE)
evaluateResults(results, settings)
# coverageCi = 0.953, precisionCi = 30.8, coveragePi = 0.917, precisionPi = 10.6, coverageTau = NA, precisionTau = NA

results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings, methodFunction = applyCurrentApproach, approximation = "normal", bayesian = TRUE)
evaluateResults(results, settings)
# coverageCi = 0.952, precisionCi = 30.7, coveragePi = 0.922, precisionPi = 10.7, coverageTau = NA, precisionTau = NA

results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings, methodFunction = applyGeneralizedModel, bayesian = FALSE)
evaluateResults(results, settings)
# coverageCi = 0.937, precisionCi = 32.8, coveragePi = NA, precisionPi = NA, coverageTau = 0.99, precisionTau = 1.08e-08

results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings, methodFunction = applyGeneralizedModel, bayesian = TRUE)
evaluateResults(results, settings)
# coverageCi = 0.949, precisionCi = 29.9, coveragePi = 0.944, precisionPi = 9.19, coverageTau = 0.967, precisionTau = 89.8

# Simulation with less data, high consistency, more iterations --------------------------------------
settings <- createSimulationSettings(
  esSettings = EvidenceSynthesis::createSimulationSettings(
    nSites = 5,
    n = 5000,
    minBackgroundHazard = 2e-05,
    maxBackgroundHazard = 2e-04,
    hazardRatio = 2,
    randomEffectSd = 0.25
  ),
  nNegativeControls = 50,
  nOutcomesOfInterest = 10,
  cpcConsistency = 25
)
plotSystematicErrorDistributions(settings)

results <- ParallelLogger::clusterApply(cluster, 1:200, simulateOne, settings = settings, methodFunction = applyNaiveApproach, bayesian = FALSE)
evaluateResults(results, settings)
# coverageCi = 0.837, precisionCi = 34.9, coveragePi = 0.877, precisionPi = 8.76, coverageTau = 0.966, precisionTau = 19.7

results <- ParallelLogger::clusterApply(cluster, 1:200, simulateOne, settings = settings, methodFunction = applyCurrentApproach, approximation = "normal", bayesian = TRUE)
evaluateResults(results, settings)
# coverageCi = 0.916, precisionCi = 22.6, coveragePi = 0.97, precisionPi = 5.73, coverageTau = NA, precisionTau = NA

results <- ParallelLogger::clusterApply(cluster, 1:200, simulateOne, settings = settings, methodFunction = applyGeneralizedModel, bayesian = TRUE)
evaluateResults(results, settings)
# coverageCi = 0.954, precisionCi = 17, coveragePi = 0.975, precisionPi = 5.38, coverageTau = 0.997, precisionTau = 39.1

# Simulation with more databases of varying sizes, high consistency, more iterations --------------------------------------
settings <- createSimulationSettings(
  esSettings = EvidenceSynthesis::createSimulationSettings(
    nSites = 10,
    n = runif(10, 2000, 6000),
    minBackgroundHazard = 2e-05,
    maxBackgroundHazard = 2e-04,
    hazardRatio = 2,
    randomEffectSd = 0.25
  ),
  nNegativeControls = 50,
  nOutcomesOfInterest = 10,
  cpcConsistency = 25
)
plotSystematicErrorDistributions(settings)

results <- ParallelLogger::clusterApply(cluster, 1:200, simulateOne, settings = settings, methodFunction = applyNaiveApproach, bayesian = FALSE)
evaluateResults(results, settings)
# coverageCi = 0.786, precisionCi = 64.5, coveragePi = 0.805, precisionPi = 17.4, coverageTau = 0.956, precisionTau = 44.6

results <- ParallelLogger::clusterApply(cluster, 1:200, simulateOne, settings = settings, methodFunction = applyCurrentApproach, approximation = "normal", bayesian = TRUE)
evaluateResults(results, settings)
# coverageCi = 0.915, precisionCi = 32.8, coveragePi = 0.941, precisionPi = 8.89, coverageTau = NA, precisionTau = NA

results <- ParallelLogger::clusterApply(cluster, 1:200, simulateOne, settings = settings, methodFunction = applyCurrentApproach, approximation = "non-normal", bayesian = TRUE)
evaluateResults(results, settings)
# coverageCi = 0.922, precisionCi = 29.3, coveragePi = 0.945, precisionPi = 8.26, coverageTau = NA, precisionTau = NA

results <- ParallelLogger::clusterApply(cluster, 1:200, simulateOne, settings = settings, methodFunction = applyGeneralizedModel, bayesian = TRUE)
evaluateResults(results, settings)
# coverageCi = 0.926, precisionCi = 28.7, coveragePi = 0.949, precisionPi = 8.43, coverageTau = 0.997, precisionTau = 64.5

ParallelLogger::stopCluster(cluster)
