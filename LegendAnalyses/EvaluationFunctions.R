
# row = list(targetId = 102100000, comparatorId = 202100000, analysisId = 7); methodFunction = applyNaiveApproach; args = list(bayesian = TRUE)
# row = tcas[[1]]
estimateForTca <- function(row, methodFunction, ...) {
  cacheFolder <- file.path("LegendAnalyses", "cache")
  if (!dir.exists(cacheFolder)) {
    dir.create(cacheFolder)
  }
  dummy <- list(...)
  dummy$methodFunction <- methodFunction
  estimatesCacheFilename <- file.path(cacheFolder, sprintf("Estimates_t%d_c%d_a%d_m%s.rds", row$targetId, row$comparatorId, row$analysisId, rlang::hash(dummy)))
  if (file.exists(estimatesCacheFilename)) {
    maEstimates <- readRDS(estimatesCacheFilename)
  } else {
    estimates <- readRDS("LegendAnalyses/estimates.rds") |>
      inner_join(row, copy = TRUE, by = join_by(targetId, comparatorId, analysisId))
    profiles <- readRDS("LegendAnalyses/profiles.rds") |>
      inner_join(row, copy = TRUE, by = join_by(targetId, comparatorId, analysisId))
    
    # profile = profiles[1, ]
    expandProfile <- function(profile) {
      profile <- tibble(
        databaseId = profile$databaseId,
        outcomeId = profile$outcomeId, 
        analysisId = profile$analysisId,
        point = as.numeric(strsplit(profile$point, ";")[[1]]),
        value = as.numeric(strsplit(profile$value, ";")[[1]])
      )
      return(profile)
    }
    profiles <- bind_rows(lapply(split(profiles, seq_len(nrow(profiles))), expandProfile))
    
    # DatabaseId needs to be sequential integer:
    databaseIds <- estimates |>
      distinct(databaseId) |>
      pull()
    estimates <- estimates |>
      mutate(databaseId = match(databaseId, databaseIds))
    profiles <- profiles |>
      mutate(databaseId = match(databaseId, databaseIds))
    
    # Leave-one-out
    outcomeIds <- estimates |>
      distinct(outcomeId) |>
      pull(outcomeId)
    # Create dummy simulation settings so we can re-use simulation functions:
    settings <- list(
      esSettings = list(nSites = estimates |> 
                          distinct(databaseId) |> 
                          count() |>
                          pull()),
      nNegativeControls = length(outcomeIds) - 1,
      nOutcomesOfInterest = 1
    )
    
    maEstimates <- list()
    # i = 7
    for (i in seq_along(outcomeIds)) {
      outcomeId <- outcomeIds[i]
      
      # Need to renumber outcome IDs for functions, so NCs have IDs 1...n:
      ncIds <- outcomeIds[!outcomeIds == outcomeId]
      tempEstimates <- estimates |>
        mutate(outcomeId = if_else(outcomeId == !!outcomeId, length(outcomeIds), match(outcomeId, ncIds)))
      tempProfiles <- profiles |>
        mutate(outcomeId = if_else(outcomeId == !!outcomeId, length(outcomeIds), match(outcomeId, ncIds)))
      data <- list(
        normalApproximations = tempEstimates,
        nonNormalApproximations = tempProfiles
      )
      args <- list(...)
      args$data <- data
      args$settings <- settings
      maEstimates[[i]] <- do.call(methodFunction, args) 
    }
    maEstimates <- bind_rows(maEstimates) |>
      bind_cols(row)
    saveRDS(maEstimates, estimatesCacheFilename)
  }
  return(maEstimates)
}

estimateLeaveOneOut <- function(cluster, methodFunction, ...) {
  snow::clusterExport(cluster, c("estimateForTca"))
  tcas <- readRDS("LegendAnalyses/estimates.rds") |>
    distinct(targetId, targetName, comparatorId, comparatorName, analysisId, analysisName) |>
    group_by(row_number()) |>
    group_split()
  
  estimates <- ParallelLogger::clusterApply(cluster, tcas, estimateForTca, methodFunction = methodFunction, ...)
  return(estimates)
}

evaluateLegendResults <- function(results, perTca = FALSE) {
  ncSettings <-   list(esSettings = list(hazardRatio = 1, 
                                         randomEffectSd = 0))
  class(ncSettings$esSettings) <- "simulationSettings"
  
  if (perTca) {
    # result = results[[1]]
    for (result in results) {
      writeLines(sprintf("\n%s vs %s, %s:", result$targetName[1], result$comparatorName[1], result$analysisName[1]))
      evaluateResults(result, ncSettings)
    }
  }
  resultsPerAnalysis <- results |>
    bind_rows(results) |>
    group_by(analysisName) |>
    group_split()
  
  for (result in resultsPerAnalysis) {
    writeLines(sprintf("\nOverall %s:", result$analysisName[1]))
    evaluateResults(result, ncSettings)
  }
  
}


