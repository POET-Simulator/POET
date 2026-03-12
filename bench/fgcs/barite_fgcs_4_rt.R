iterations <- 5000

dt <- 2000

list(
    timesteps = rep(dt, iterations),
    store_result = TRUE,
    out_save = c(seq(1, 10), seq(10, 100, by= 10), seq(200, iterations, by=100))
    ## out_save = seq_len(iterations)
)
