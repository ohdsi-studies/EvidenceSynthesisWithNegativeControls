library(DatabaseConnector)
library(dplyr)

connectionDetails <- createConnectionDetails(
  dbms = "postgresql",
  server = paste(keyring::key_get("ohdsiPostgresServer"),
                 keyring::key_get("ohdsiPostgresShinyDatabase"), sep = "/"),
  user = keyring::key_get("ohdsiPostgresUser"),
  password = keyring::key_get("ohdsiPostgresPassword")
)
schema <- "legendt2dm_class_results_v2"

prettyDbNames <- tibble(
  oldName = c("CCAE",
              "MDCD",
              "MDCR",
              "OptumEHR",
              "OptumDod",
              "VA-OMOP"),
  newName = c("Merative CCAE",
              "Merative MDCD",
              "Merative MDCR",
              "Optum EHR",
              "Optum Clinformatics",
              "Veterans Affairs")
)

# executeSql(connection, "COMMIT;")
connection <- connect(connectionDetails)


sql <- "
SELECT database_id,
    target_id,
    target.exposure_name AS target_name,
    comparator_id,
    comparator.exposure_name AS comparator_name,
    cohort_method_result.outcome_id,
    outcome_name,
    analysis_id,
    CASE WHEN analysis_id = 7 THEN 'unadjusted' ELSE 'PS matched' END AS analysis_name,
    log_rr,
    se_log_rr
FROM @schema.cohort_method_result
INNER JOIN @schema.exposure_of_interest target
  ON target_id = target.exposure_id
INNER JOIN @schema.exposure_of_interest comparator
  ON comparator_id = comparator.exposure_id
INNER JOIN @schema.negative_control_outcome
  ON cohort_method_result.outcome_id = negative_control_outcome.outcome_id
WHERE se_log_rr IS NOT NULL
    AND analysis_id IN (7, 8)
    AND target.exposure_name LIKE '%main ot2'
    AND comparator.exposure_name LIKE '%main ot2'
    AND database_id NOT LIKE 'Meta-analysis%';
"
estimates <- renderTranslateQuerySql(connection = connection,
                                     sql = sql,
                                     schema = schema,
                                     snakeCaseToCamelCase = TRUE)
estimates <- estimates |>
  mutate(databaseId = as.factor(databaseId),
         analysisName = as.factor(analysisName),
         targetName = as.factor(gsub(" main ot2", "", targetName)),
         comparatorName = as.factor(gsub(" main ot2", "", comparatorName)),
         outcomeName = as.factor(gsub("_", " ", gsub("outcome/", "", outcomeName))))
estimates <- estimates |>
  inner_join(prettyDbNames, by = join_by(databaseId == oldName)) |>
  mutate(databaseId = newName) |>
  select(-newName)

sql <- "
SELECT database_id,
    target_id,
    comparator_id,
    likelihood_profile.outcome_id,
    analysis_id,
    point,
    value
FROM @schema.likelihood_profile
INNER JOIN @schema.exposure_of_interest target
  ON target_id = target.exposure_id
INNER JOIN @schema.exposure_of_interest comparator
  ON comparator_id = comparator.exposure_id
INNER JOIN @schema.negative_control_outcome
  ON likelihood_profile.outcome_id = negative_control_outcome.outcome_id
WHERE analysis_id IN (7, 8)
    AND target.exposure_name LIKE '%main ot2'
    AND comparator.exposure_name LIKE '%main ot2';
"
profiles <- renderTranslateQuerySql(connection = connection,
                                     sql = sql,
                                     schema = schema,
                                     snakeCaseToCamelCase = TRUE)
profiles <- profiles |>
  inner_join(prettyDbNames, by = join_by(databaseId == oldName)) |>
  mutate(databaseId = newName) |>
  select(-newName)

disconnect(connection)
saveRDS(estimates, "LegendAnalyses/estimates.rds")
saveRDS(profiles, "LegendAnalyses/profiles.rds")
