iterations <- 50
dt <- 100

list(
    timesteps = rep(dt, iterations),
    out_save  = c(1, seq(5, iterations, by = 5)),
    store_result = TRUE
)

