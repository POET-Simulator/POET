## Time-stamp: "Last modified 2024-01-12 11:35:11 delucia"

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
  "Charge" = -1.216307845207e-09,
  "Ba" = 1.E-12,
  "Cl" = 2.E-12,
  "S(6)" = 6.204727095976e-04,
  "Sr" = 6.204727095976e-04,
  "Barite" = 0.001,
  "Celestite" = 1
)

grid <- list(
  n_cells = c(n, m),
  s_cells = c(1, 1),
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
    "H" = 110.0124,
    "O" = 55.5087,
    "Charge" = -1.216307845207e-09,
    "Ba" = 1.E-12,
    "Cl" = 2.E-12,
    "S(6)" = 6.204727095976e-04,
    "Sr" = 6.204727095976e-04
)

injection_diff <- list(
    list(
        "H" = 111.0124,
        "O" = 55.50622,
        "Charge" = -3.336970273297e-08,
        "Ba" = 0.1,
        "Cl" = 0.2,
        "S(6)"  = 0,
        "Sr" = 0)
)

## diffusion coefficients
alpha_diffu <- c(
  "H"  =  1E-06,
  "O"  =  1E-06,
  "Charge" = 1E-06,
  "Ba" = 1E-06,
  "Cl" = 1E-06,
  "S(6)"  = 1E-06,
  "Sr" = 1E-06
)

boundary <- list(
  "N" = c(1,1, rep(0, n-2)),
  "E" = rep(0, n),
  "S" = rep(0, n),
  "W" = c(1,1, rep(0, n-2))
)

diffu_list <- names(alpha_diffu)

vecinj <- do.call(rbind.data.frame, injection_diff)
names(vecinj) <- names(init_diffu)

diffusion <- list(
  init = as.data.frame(init_diffu, check.names = FALSE),
  vecinj = vecinj,
  vecinj_index = boundary,
  alpha = alpha_diffu
)

#################################################################
##                          Section 3                          ##
##                  Chemistry module (Phreeqc)                 ##
#################################################################

## # Needed when using DHT
dht_species <- c(
  "H" = 7,
  "O" = 7,
  "Charge" = 4,
  "Ba" = 7,
  "Cl" = 7,
  "S(6)" = 7,
  "Sr" = 7,
  "Barite" = 4,
  "Celestite" = 4
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
