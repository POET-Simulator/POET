iterations <- 500
dt <- 50

list(
    timesteps = rep(dt, iterations),
    store_result = TRUE,
    out_save = c(5, iterations, by = 5)
)