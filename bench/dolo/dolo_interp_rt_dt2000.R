iterations <- 2000
dt <- 200

out_save <- c(1, 5, 10, seq(20, iterations, by = 20))

list(
    timesteps = rep(dt, iterations),
    store_result = TRUE,
    out_save = out_save
)
