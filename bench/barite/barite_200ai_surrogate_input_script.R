## load a pretrained model from tensorflow file
## Use the global variable "ai_surrogate_base_path" when using file paths
## relative to the input script
initiate_model <- function() {
  init_model <- normalizePath(paste0(ai_surrogate_base_path,
                                     "model_min_max_float64.keras"))
  return(load_model_tf(init_model))
}

scale_min_max <- function(x, min, max, backtransform) {
  if (backtransform) {
    return((x * (max - min)) + min)
  } else {
    return((x - min) / (max - min))
  }
}

preprocess <- function(df, backtransform = FALSE, outputs = FALSE) {
  minmax_file <- normalizePath(paste0(ai_surrogate_base_path,
                                      "min_max_bounds.rds"))
  global_minmax <- readRDS(minmax_file)
  for (column in colnames(df)) {
    df[column] <- lapply(df[column],
                         scale_min_max,
                         global_minmax$min[column],
                         global_minmax$max[column],
                         backtransform)
  }
  return(df)
}

mass_balance <- function(predictors, prediction) {
  dBa <- abs(prediction$Ba + prediction$Barite -
             predictors$Ba - predictors$Barite)
  dSr <- abs(prediction$Sr + prediction$Celestite -
             predictors$Sr - predictors$Celestite)
  return(dBa + dSr)
}

validate_predictions <- function(predictors, prediction) {
  epsilon <- 0.000000003
  mb <- mass_balance(predictors, prediction)
  msgm("Mass balance mean:", mean(mb))
  msgm("Mass balance variance:", var(mb))
  msgm("Rows where mass balance meets threshold", epsilon, ":",
       sum(mb < epsilon))
  return(mb < epsilon)
}
