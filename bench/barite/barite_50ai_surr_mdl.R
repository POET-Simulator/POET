## Time-stamp: "Last modified 2024-05-30 13:27:06 delucia"

## load a pretrained model from tensorflow file
## Use the global variable "ai_surrogate_base_path" when using file paths
## relative to the input script
initiate_model <- function() {
    require(keras3)
    require(tensorflow)
    init_model <- normalizePath(paste0(ai_surrogate_base_path,
                                       "barite_50ai_all.keras"))
    Model <- keras3::load_model(init_model)
    msgm("Loaded model:")
    print(str(Model))
    return(Model)
}

scale_min_max <- function(x, min, max, backtransform) {
    if (backtransform) {
        return((x * (max - min)) + min)
    } else {
        return((x - min) / (max - min))
    }
}

minmax <- list(min = c(H = 111.012433592824, O = 55.5062185549492, Charge = -3.1028354471876e-08, 
                       Ba = 1.87312878574393e-141, Cl = 0, `S(6)` = 4.24227510643685e-07, 
                       Sr = 0.00049382996130541, Barite = 0.000999542409828586, Celestite = 0.244801877115968),
               max = c(H = 111.012433679682, O = 55.5087003521685, Charge = 5.27666636082035e-07, 
                       Ba = 0.0908849779513762, Cl = 0.195697626449355, `S(6)` = 0.000620774752665846, 
                       Sr = 0.0558680070692722, Barite = 0.756779139057097, Celestite = 1.00075422160624
                       ))

preprocess <- function(df) {
    if (!is.data.frame(df))
        df <- as.data.frame(df, check.names = FALSE)
    
    as.data.frame(lapply(colnames(df),
                         function(x) scale_min_max(x=df[x],
                                                   min=minmax$min[x],
                                                   max=minmax$max[x],
                                                   backtransform=FALSE)),
                  check.names = FALSE)
}

postprocess <- function(df) {
    if (!is.data.frame(df))
        df <- as.data.frame(df, check.names = FALSE)
    
    as.data.frame(lapply(colnames(df),
                         function(x) scale_min_max(x=df[x],
                                                   min=minmax$min[x],
                                                   max=minmax$max[x],
                                                   backtransform=TRUE)),
                  check.names = FALSE)
}

mass_balance <- function(predictors, prediction) {
    dBa <- abs(prediction$Ba + prediction$Barite -
               predictors$Ba - predictors$Barite)
    dSr <- abs(prediction$Sr + prediction$Celestite -
               predictors$Sr - predictors$Celestite)
    return(dBa + dSr)
}

validate_predictions <- function(predictors, prediction) {
    epsilon <- 1E-7
    mb <- mass_balance(predictors, prediction)
    msgm("Mass balance mean:", mean(mb))
    msgm("Mass balance variance:", var(mb))
    ret <- mb < epsilon
    msgm("Rows where mass balance meets threshold", epsilon, ":",
         sum(ret))
    return(ret)
}

training_step <- function(model, predictor, target, validity) {
  msgm("Starting incremental training:")

  ## x <- as.matrix(predictor)
  ## y <- as.matrix(target[colnames(x)])

  history <- model %>% keras3::fit(x = data.matrix(predictor),
                                   y = data.matrix(target),
                                   epochs = 10, verbose=1)

  keras3::save_model(model,
                     filepath = paste0(out_dir, "/current_model.keras"),
                     overwrite=TRUE)
  return(model)
}
