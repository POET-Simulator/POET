scale_min_max <- function(x, min, max, backtransform) {
    if (backtransform) {
        return((x * (max - min)) + min)
    } else {
        return((x - min) / (max - min))
    }
}

scale_standardizer <- function(x, mean, scale, backtransform) {
    if(backtransform){
        return(x * scale + mean)
    }
    else{
        return((x-mean) / scale)
    }
}

standard <- list(mean = c(H = 111.0124335959659, O=55.5065739707202, 'C(-4)'=1.5788555695339323e-15, 'C(4)'=0.00011905649680154037,
                         Ca= 0.00012525858032576948, Cl=0.00010368471137502122, Mg=4.5640346338857756e-05, Calcite_kin=0.0001798444527389999, 
                         Dolomite_kin=7.6152065281986634e-06), 
                scale = c(H=1.0, O=3.54850912318837e-05, 'C(-4)'=2.675559053860093e-14, 'C(4)'=1.1829735682920146e-05, Ca=1.207381343127647e-05, Cl=0.00024586541554245565,
                         Mg=0.00011794307217698012, Calcite_kin=5.946457663332385e-05, Dolomite_kin=2.688201435907049e-05))


ai_surrogate_species_input = c("H", "O", "C(-4)", "C(4)", "Ca", "Cl", "Mg", "Calcite_kin", "Dolomite_kin")
ai_surrogate_species_output = c("H", "O", "C(-4)", "C(4)", "Ca", "Mg", "Calcite_kin", "Dolomite_kin")


threshold <- list(species = "Cl", value = 2E-10)

preprocess <- function(df) {
    if (!is.data.frame(df))
        df <- as.data.frame(df, check.names = FALSE)
    
    as.data.frame(lapply(colnames(df),
                         function(x) scale_standardizer(x=df[x],
                                                   mean=standard$mean[x],
                                                   scale=standard$scale[x],
                                                   backtransform=FALSE)),
                  check.names = FALSE)
}

postprocess <- function(df) {
    if (!is.data.frame(df))
        df <- as.data.frame(df, check.names = FALSE)
    
    as.data.frame(lapply(colnames(df),
                         function(x) scale_standardizer(x=df[x],
                                                   mean=standard$mean[x],
                                                   scale=standard$scale[x],
                                                   backtransform=TRUE)),
                  check.names = FALSE)
}

mass_balance <- function(predictors, prediction) {
    dCa <- abs(prediction$Ca + prediction$Calcite_kin + prediction$Dolomite_kin -
               predictors$Ca - predictors$Calcite_kin - predictors$Dolomite_kin)
    dC <- abs(prediction$'C(-4)' + prediction$'C(4)' + prediction$Calcite_kin + 2 * prediction$Dolomite_kin 
                - predictors$'C(-4)' - predictors$'C(4)' - predictors$Calcite_kin - 2 * predictors$Dolomite_kin)
    dMg <- abs(prediction$Mg + prediction$Dolomite_kin - 
              predictors$Mg - predictors$Dolomite_kin)
    return(dCa + dC + dMg)
}

validate_predictions <- function(predictors, prediction) {
    epsilon <- 1E-8
    mb <- mass_balance(predictors, prediction)
    msgm("Mass balance mean:", mean(mb))
    msgm("Mass balance variance:", var(mb))
    ret <- mb < epsilon
    msgm("Rows where mass balance meets threshold", epsilon, ":",
         sum(ret))
    return(ret)
}