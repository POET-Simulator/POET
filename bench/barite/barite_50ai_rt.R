iterations <- 1000

dt <- 200

list(
    timesteps = rep(dt, iterations),
    store_result = TRUE,
    out_save = c(1, 5, seq(20, iterations, by=20))
)
