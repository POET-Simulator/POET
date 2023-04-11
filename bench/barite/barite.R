## Time-stamp: "Last modified 2023-04-11 14:29:51 delucia"

database <- normalizePath("./db_barite.dat")
input_script <- normalizePath("./barite.pqi")

#################################################################
##                          Section 1                          ##
##                     Grid initialization                     ##
#################################################################

n <- 20
m <- 20

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
  "H" = 110.0124,
  "O" = 55.5087,
  "Charge" = -1.217E-09,
  "Ba" = 1.E-10,
  "Cl" = 1.E-10,
  "S" = 6.205E-4,
  "Sr" = 6.205E-4
)

injection_diff <- list(
    list(
        "H" = 111.0124,
        "O" = 55.50622,
        "Charge" = -3.337E-08,
        "Ba" = 0.1,
        "Cl" = 0.2,
        "S"  = 0,
        "Sr" = 0)
)

## diffusion coefficients
alpha_diffu <- c(
  "H"  =  1E-06,
  "O"  =  1E-06,
  "Charge" = 1E-06,
  "Ba" = 1E-06,
  "Cl" = 1E-06,
  "S"  = 1E-06,
  "Sr" = 1E-06
)

## vecinj_inner <- list(
##   l1 = c(1,20,20),
##   l2 = c(2,80,80),
##   l3 = c(2,60,80)
## )

boundary <- list(
  "N" = rep(1, n),
##  "N" = rep(0, n),
  "E" = rep(0, n),
  "S" = rep(0, n),
  "W" = rep(0, n)
)

diffu_list <- names(alpha_diffu)

diffusion <- list(
  init = init_diffu,
  vecinj = do.call(rbind.data.frame, injection_diff),
#  vecinj_inner = vecinj_inner,
  vecinj_index = boundary,
  alpha = alpha_diffu
)

#################################################################
##                          Section 3                          ##
##                  Chemistry module (Phreeqc)                 ##
#################################################################


## # Needed when using DHT
signif_vector <- c(9, 9, 10, 5, 5, 5, 5, 5, 5)
prop_type <- c("", "", "", "act", "act", "act", "act", "", "")
prop <- names(init_cell)

chemistry <- list(
  database = database,
  input_script = input_script
)

#################################################################
##                          Section 4                          ##
##              Putting all those things together              ##
#################################################################


iterations <- 4
dt <- 100

setup <- list(
  grid = grid,
  diffusion = diffusion,
  chemistry = chemistry,
  iterations = iterations,
  timesteps = rep(dt, iterations),
  store_result = TRUE,
  out_save = seq(1, iterations)
)
