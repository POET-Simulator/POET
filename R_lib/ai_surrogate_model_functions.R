## This file contains default function implementations for the ai surrogate.
## To use pre-/postprocessing it is recommended to override these functions
## with custom implementations via the input script. The path to the R-file
## See the barite_50.R file as an example and the general README for more
## information.

preprocess <- function(df, backtransform = FALSE, outputs = FALSE) {
  return(df)
}

postprocess <- function(df, backtransform = TRUE, outputs = TRUE) {
  return(df)
}

set_valid_predictions <- function(temp_field, prediction, validity) {
  temp_field[validity == 1, ] <- prediction[validity == 1, ]
  return(temp_field)
}

get_invalid_values <- function(df, validity) {
  return(df[validity == 0, ])
}

set_field <- function(temp_field, columns, rows, column_name_limit,
                      byrow = FALSE) {
  temp_field <- matrix(temp_field, nrow = rows, byrow = byrow)
  temp_field <- setNames(data.frame(temp_field), columns)
  temp_field <- temp_field[column_name_limit]
  return(temp_field)
}
