## Time-stamp: "Last modified 2023-01-10 13:51:40 delucia"

database <- normalizePath("../share/poet/examples/phreeqc_kin.dat")
input_script <- normalizePath("../share/poet/bench/dolo_inner.pqi")

#################################################################
##                          Section 1                          ##
##                     Grid initialization                     ##
#################################################################

n <- 100
m <- 100

types <- c("scratch", "phreeqc", "rds")

init_cell <- list(
  "H" = 110.683,
  "O" = 55.3413,
  "Charge" = -5.0822e-19,
  "C" = 1.2279E-4,
  "Ca" = 1.2279E-4,
  "Cl" = 0,
  "Mg" = 0,
  "O2g" = 0.499957,
  "Calcite" = 2.07e-4,
  "Dolomite" = 0
)

grid <- list(
  n_cells = c(n, m),
  s_cells = c(1, 1),
  type = types[1],
  init_cell = as.data.frame(init_cell),
  props = names(init_cell),
  database = database,
  input_script = input_script
)


##################################################################
##                          Section 2                           ##
##         Diffusion parameters and boundary conditions         ##
##################################################################

## initial conditions
init_diffu <- c(
  "H" = 110.683,
  "O" = 55.3413,
  "Charge" = -5.0822e-19,
  "C" = 1.2279E-4,
  "Ca" = 1.2279E-4,
  "Cl" = 0,
  "Mg" = 0
)

## diffusion coefficients
alpha_diffu <- c(
  "H" = 1E-6,
  "O" = 1E-6,
  "Charge" = 1E-6,
  "C" = 1E-6,
  "Ca" = 1E-6,
  "Cl" = 1E-6,
  "Mg" = 1E-6
)

## list of boundary conditions/inner nodes
vecinj_diffu <- list(
    list(
        "H" = 110.683,
        "O" = 55.3413,
        "Charge" = 1.90431e-16,
        "C" = 0,
        "Ca" = 0,
        "Cl" = 0.002,
        "Mg" = 0.001
    ),
    list(
        "H" = 110.683,
        "O" = 55.3413,
        "Charge" = 1.90431e-16,
        "C" = 0,
        "Ca" = 0.0,
        "Cl" = 0.004,
        "Mg" = 0.002
    )
)

vecinj_inner <- list(
  l1 = c(1,20,20),
  l2 = c(2,80,80),
  l3 = c(2,60,80)
)

boundary <- list(
#  "N" = c(1, rep(0, n-1)),
  "N" = rep(0, n),
  "E" = rep(0, n),
  "S" = rep(0, n),
  "W" = rep(0, n)
)

diffu_list <- names(alpha_diffu)

diffusion <- list(
  init = init_diffu,
  vecinj = do.call(rbind.data.frame, vecinj_diffu),
  vecinj_inner = vecinj_inner,
  vecinj_index = boundary,
  alpha = alpha_diffu
)

#################################################################
##                          Section 3                          ##
##                  Chemistry module (Phreeqc)                 ##
#################################################################


## # Needed when using DHT
signif_vector <- c(10, 10, 2, 5, 5, 5, 5, 0, 5, 5)
prop_type <- c("", "", "", "act", "act", "act", "act", "ignore", "", "")
prop <- names(init_cell)

chemistry <- list(
  database = database,
  input_script = input_script
)

#################################################################
##                          Section 4                          ##
##              Putting all those things together              ##
#################################################################


iterations <- 1000
dt <- 200

setup <- list(
  grid = grid,
  diffusion = diffusion,
  chemistry = chemistry,
  iterations = iterations,
  timesteps = rep(dt, iterations),
  store_result = TRUE,
  out_save = c(5, iterations)
)
