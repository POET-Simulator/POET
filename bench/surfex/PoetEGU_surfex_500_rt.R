iterations <- 200
dt <- 1000

out_save <- c(1, 2, seq(5, iterations, by=5))
## out_save <- seq(1, iterations)

list(
    timesteps = rep(dt, iterations),
    store_result = TRUE,
    out_save = out_save
)
