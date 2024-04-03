iterations <- 500
dt <- 500

list(
    timesteps = rep(dt, iterations),
    store_result = TRUE,
    out_save = seq(50, iterations, by = 50)
)
