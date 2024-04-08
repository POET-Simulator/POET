iterations <- 500
dt <- 50

out_save <- seq(5, iterations, by = 5)

list(
    timesteps = rep(dt, iterations),
    store_result = TRUE,
    out_save = out_save
)
