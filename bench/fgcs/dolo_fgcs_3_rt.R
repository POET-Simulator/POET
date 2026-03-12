iterations <- 10000
dt <- 200

list(
    timesteps = rep(dt, iterations),
    store_result = TRUE,
    out_save = c(1,20,50, seq(100,1000, by=100), seq(2000, iterations, by=1000))
)
