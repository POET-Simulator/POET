iterations <- 100
dt <- 50

list(
    timesteps = rep(dt, iterations),
    store_result = TRUE,
    out_save = seq(5, iterations, by = 5)
)
