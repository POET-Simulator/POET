## Time-stamp: "Last modified 2023-08-16 14:57:25 mluebke"

database <- normalizePath("../share/poet/bench/dolo/phreeqc_kin.dat")
input_script <- normalizePath("../share/poet/bench/dolo/dolo_inner.pqi")

#################################################################
##                          Section 1                          ##
##                     Grid initialization                     ##
#################################################################

n <- 400
m <- 200

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
  s_cells = c(5, 2.5),
  type = types[1]
)


##################################################################
##                          Section 2                           ##
##         Diffusion parameters and boundary conditions         ##
##################################################################

## initial conditions
init_diffu <- list(
  "H" = 1.110124E+02,
  "O" = 5.550833E+01,
  "Charge" = -1.216307659761E-09,
  "C(4)" = 1.230067028174E-04,
  "Ca" = 1.230067028174E-04,
  "Cl" = 0,
  "Mg" = 0
)

## diffusion coefficients
alpha_diffu <- c(
  "H" = 1E-6,
  "O" = 1E-6,
  "Charge" = 1E-6,
  "C(4)" = 1E-6,
  "Ca" = 1E-6,
  "Cl" = 1E-6,
  "Mg" = 1E-6
)

## list of boundary conditions/inner nodes
vecinj_diffu <- list(
  list(
    "H" = 1.110124E+02,
    "O" = 5.550796E+01,
    "Charge" = -3.230390327801E-08,
    "C(4)" = 0,
    "Ca" = 0,
    "Cl" = 0.002,
    "Mg" = 0.001
  ),
  list(
    "H" = 110.683,
    "O" = 55.3413,
    "Charge" = 1.90431e-16,
    "C(4)" = 0,
    "Ca" = 0.0,
    "Cl" = 0.004,
    "Mg" = 0.002
  ),
  init_diffu
)

vecinj_inner <- list(
  l1 = c(1, floor(n / 2), floor(m / 2))
  # l2 = c(2,1400,800),
  # l3 = c(2,1600,800)
)

boundary <- list(
  #  "N" = c(1, rep(0, n-1)),
  "N" = rep(3, n),
  "E" = rep(3, m),
  "S" = rep(3, n),
  "W" = rep(3, m)
)

diffu_list <- names(alpha_diffu)

vecinj <- do.call(rbind.data.frame, vecinj_diffu)
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


## # optional when using DHT
dht_species <- c(
  "H" = 3,
  "O" = 3,
  "Charge" = 3,
  "C(4)" = 6,
  "Ca" = 6,
  "Cl" = 3,
  "Mg" = 5,
  "Calcite" = 4,
  "Dolomite" = 4
)

## # Optional when using Interpolation (example with less key species and custom
## # significant digits)

# pht_species <- c(
#  "C(4)" = 3,
#  "Ca" = 3,
#  "Mg" = 2,
#  "Calcite" = 2,
#  "Dolomite" = 2
# )

check_sign_cal_dol_dht <- function(old, new) {
  if ((old["Calcite"] == 0) != (new["Calcite"] == 0)) {
    return(TRUE)
  }
  if ((old["Dolomite"] == 0) != (new["Dolomite"] == 0)) {
    return(TRUE)
  }
  return(FALSE)
}

fuzz_input_dht_keys <- function(input) {
  return(input[names(dht_species)])
}

check_sign_cal_dol_interp <- function(to_interp, data_set) {
  data_set <- as.data.frame(do.call(rbind, data_set), check.names = FALSE, optional = TRUE)
  names(data_set) <- names(dht_species)
  cal <- (data_set$Calcite == 0) == (to_interp["Calcite"] == 0)
  dol <- (data_set$Dolomite == 0) == (to_interp["Dolomite"] == 0)

  cal_dol_same_sig <- cal == dol
  return(rev(which(!cal_dol_same_sig)))
}

check_neg_cal_dol <- function(result) {
  neg_sign <- (result["Calcite"] < 0) || (result["Dolomite"] < 0)
  return(neg_sign)
}

hooks <- list(
  dht_fill = check_sign_cal_dol_dht,
  dht_fuzz = fuzz_input_dht_keys,
  interp_pre_func = check_sign_cal_dol_interp,
  interp_post_func = check_neg_cal_dol
)

chemistry <- list(
  database = database,
  input_script = input_script,
  dht_species = dht_species,
  hooks = hooks
  #  pht_species = pht_species
)

#################################################################
##                          Section 4                          ##
##              Putting all those things together              ##
#################################################################


iterations <- 20000
dt <- 200

setup <- list(
  grid = grid,
  diffusion = diffusion,
  chemistry = chemistry,
  iterations = iterations,
  timesteps = rep(dt, iterations),
  store_result = TRUE,
  out_save = c(1, seq(50, iterations, by = 50))
)
