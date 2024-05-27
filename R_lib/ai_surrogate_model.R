## This file contains default function implementations for the ai surrogate.
## To load pretrained models, use pre-/postprocessing or change hyperparameters
## it is recommended to override these functions with custom implementations via
## the input script. The path to the R-file containing the functions mus be set
## in the variable "ai_surrogate_input_script". See the barite_200.R file as an
## example and the general README for more information.

library(keras)
library(tensorflow)

initiate_model <- function() {
  hidden_layers <- c(48, 96, 24)
  activation <- "relu"
  loss <- "mean_squared_error"

  input_length <- length(ai_surrogate_species)
  output_length <- length(ai_surrogate_species)
  ## Creates a new sequential model from scratch
  model <- keras_model_sequential()

  ## Input layer defined by input data shape
  model %>% layer_dense(units = input_length,
                        activation = activation,
                        input_shape = input_length,
                        dtype = "float32")

  for (layer_size in hidden_layers) {
    model %>% layer_dense(units = layer_size,
                          activation = activation,
                          dtype = "float32")
  }

  ## Output data defined by output data shape
  model %>% layer_dense(units = output_length,
                        activation = activation,
                        dtype = "float32")

  model %>% compile(loss = loss,
                    optimizer = "adam")
  return(model)
}

gpu_info <- function() {
  msgm(tf_gpu_configured())
}

prediction_step <- function(model, predictors) {
  prediction <- predict(model, as.matrix(predictors))
  colnames(prediction) <- colnames(predictors)
  return(as.data.frame(prediction))
}

preprocess <- function(df, backtransform = FALSE, outputs = FALSE) {
  return(df)
}

set_valid_predictions <- function(temp_field, prediction, validity) {
  temp_field[validity == 1, ] <- prediction[validity == 1, ]
  return(temp_field)
}

training_step <- function(model, predictor, target, validity) {
  msgm("Training:")

  x <- as.matrix(predictor)
  y <- as.matrix(target[colnames(x)])

  model %>% fit(x, y)

  model %>% save_model_tf(paste0(out_dir, "/current_model.keras"))
}
