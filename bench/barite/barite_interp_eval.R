## Time-stamp: "Last modified 2023-07-21 15:04:49 mluebke"

database <- normalizePath("../share/poet/bench/barite/db_barite.dat")
input_script <- normalizePath("../share/poet/bench/barite/barite.pqi")

#################################################################
##                          Section 1                          ##
##                     Grid initialization                     ##
#################################################################

n <- 400
m <- 200

types <- c("scratch", "phreeqc", "rds")

init_cell <- list(
  "H" = 110.0124,
  "O" = 55.5087,
  "Charge" = -1.217E-09,
  "Ba" = 1.E-10,
  "Cl" = 2.E-10,
  "S" = 6.205E-4,
  "Sr" = 6.205E-4,
  "Barite" = 0.001,
  "Celestite" = 1
)

grid <- list(
  n_cells = c(n, m),
  s_cells = c(n / 10, m / 10),
  type = types[1],
  init_cell = as.data.frame(init_cell, check.names = FALSE),
  props = names(init_cell),
  database = database,
  input_script = input_script
)


##################################################################
##                          Section 2                           ##
##         Diffusion parameters and boundary conditions         ##
##################################################################

## initial conditions

init_diffu <- list(
  # "H" = 110.0124,
  "H" = 0.00000028904,
  # "O" = 55.5087,
  "O" = 0.000000165205,
  # "Charge" = -1.217E-09,
  "Charge" = -3.337E-08,
  "Ba" = 1.E-10,
  "Cl" = 1.E-10,
  "S(6)" = 6.205E-4,
  "Sr" = 6.205E-4
)

injection_diff <- list(
  list(
    # "H" = 111.0124,
    "H" = 0.0000002890408,
    # "O" = 55.50622,
    "O" = 0.00002014464,
    # "Charge" = -3.337E-08,
    "Charge" = -3.337000004885E-08,
    "Ba" = 0.1,
    "Cl" = 0.2,
    "S(6)" = 0,
    "Sr" = 0
  )
)

## diffusion coefficients
alpha_diffu <- c(
  "H" = 1E-06,
  "O" = 1E-06,
  "Charge" = 1E-06,
  "Ba" = 1E-06,
  "Cl" = 1E-06,
  "S(6)" = 1E-06,
  "Sr" = 1E-06
)

vecinj_inner <- list(
  l1 = c(1, floor(n / 2), floor(m / 2))
  ##   l2 = c(2,80,80),
  ##   l3 = c(2,60,80)
)

boundary <- list(
  #  "N" = rep(1, n),
  "N" = rep(0, n),
  "E" = rep(0, n),
  "S" = rep(0, n),
  "W" = rep(0, n)
)

diffu_list <- names(alpha_diffu)

vecinj <- do.call(rbind.data.frame, injection_diff)
names(vecinj) <- names(init_diffu)

diffusion <- list(
  init = as.data.frame(init_diffu, check.names = FALSE),
  vecinj = vecinj,
  vecinj_inner = vecinj_inner,
  vecinj_index = boundary,
  alpha = alpha_diffu
)

#################################################################
##                          Section 3                          ##
##                  Chemistry module (Phreeqc)                 ##
#################################################################

## # Needed when using DHT
dht_species <- c(
  "H" = 10,
  "O" = 10,
  "Charge" = 3,
  "Ba" = 5,
  "Cl" = 5,
  "S(6)" = 5,
  "Sr" = 5
)

chemistry <- list(
  database = database,
  input_script = input_script,
  dht_species = dht_species
)

#################################################################
##                          Section 4                          ##
##              Putting all those things together              ##
#################################################################


iterations <- 200
dt <- 250

setup <- list(
  grid = grid,
  diffusion = diffusion,
  chemistry = chemistry,
  iterations = iterations,
  timesteps = rep(dt, iterations),
  store_result = TRUE,
  out_save = seq(1, iterations)
)
