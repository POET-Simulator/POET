iterations <- 2000
dt <- 200

out_save <- c(1, 10, 20, seq(40, iterations, by = 40))

list(
    timesteps = rep(dt, iterations),
    store_result = TRUE,
    out_save = out_save
)
