


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

standard <- list(mean = c(H = 111.01243361730982, O= 55.50673140754027, Ba= 0.0016161137065825058,
                         Cl= 0.0534503766678322, S=0.00012864849674669584, Sr=0.0252377348949622, 
                         Barite_kin=0.05292312117000998, Celestite_kin=0.9475491659328229), 
                scale = c(H=1.0, O=0.00048139729680698453, Ba=0.008945717576237102, Cl=0.03587363709464328,
                         S=0.00012035100591827131, Sr=0.01523052668095922, Barite_kin=0.21668648247230615,
                         Celestite_kin=0.21639449682671968))

minmax <- list(min = c(H = 111.012433592824, O = 55.5062185549492, Charge = -3.1028354471876e-08, 
                       Ba = 1.87312878574393e-141, Cl = 0, `S(6)` = 4.24227510643685e-07, 
                       Sr = 0.00049382996130541, Barite = 0.000999542409828586, Celestite = 0.244801877115968),
               max = c(H = 111.012433679682, O = 55.5087003521685, Charge = 5.27666636082035e-07, 
                       Ba = 0.0908849779513762, Cl = 0.195697626449355, `S(6)` = 0.000620774752665846, 
                       Sr = 0.0558680070692722, Barite = 0.756779139057097, Celestite = 1.00075422160624
                       ))

ai_surrogate_species_input = c("H", "O", "Ba", "Cl", "S", "Sr", "Barite_kin", "Celestite_kin")
ai_surrogate_species_output = c("O", "Ba", "S", "Sr", "Barite_kin", "Celestite_kin")



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
    dBa <- abs(prediction$Ba + prediction$Barite_kin -
               predictors$Ba - predictors$Barite_kin)
    dSr <- abs(prediction$Sr + prediction$Celestite_kin -
               predictors$Sr - predictors$Celestite_kin)
    dS <- abs(prediction$S + prediction$Celestite_kin + prediction$Barite_kin - 
              predictors$S - predictors$Celestite_kin - predictors$Barite_kin)
    return(dBa + dSr + dS)
}

validate_predictions <- function(predictors, prediction) {
    epsilon <- 1E-5
    mb <- mass_balance(predictors, prediction)
    msgm("Mass balance mean:", mean(mb))
    msgm("Mass balance variance:", var(mb))
    ret <- mb < epsilon
    msgm("Rows where mass balance meets threshold", epsilon, ":",
         sum(ret))
    return(ret)
}
