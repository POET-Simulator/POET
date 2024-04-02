## Time-stamp: "Last modified 2023-09-06 10:58:16 mluebke"

database <- normalizePath("../share/poet/bench/dolo/phreeqc_kin.dat")
input_script <- normalizePath("../share/poet/bench/dolo/dolo_inner.pqi")

#################################################################
##                          Section 1                          ##
##                     Grid initialization                     ##
#################################################################

n <- 1500
m <- 500

types <- c("scratch", "phreeqc", "rds")

init_cell <- list(
  "H" = 110.683,
  "O" = 55.3413,
  "Charge" = -5.0822e-19,
  "C(4)" = 1.2279E-4,
  "Ca" = 1.2279E-4,
  "Cl" = 0,
  "Mg" = 0,
  "O2g" = 0.499957,
  "Calcite" = 2.07e-4,
  "Dolomite" = 0
)

grid <- list(
  n_cells = c(n, m),
  s_cells = c(n, m),
  type = types[1]
)

##################################################################
##                          Section 2                           ##
##                     Advection parameters                     ##
##################################################################

## initial conditions

## HACK: We need the chemical initialization here, as chem module initialization
## depends on transport until now. This will change in the future.
init_adv <- c(
  "H" = 110.124,
  "O" = 55.0622,
  "Charge" = -1.216307659761E-09,
  "C(4)" = 1.230067028174E-04,
  "Ca" = 1.230067028174E-04,
  "Cl" = 0,
  "Mg" = 0
)

## list of boundary conditions/inner nodes
vecinj_adv <- list(
  list(
    "H" = 110.683,
    "O" = 55.3413,
    "Charge" = 1.90431e-16,
    "C(4)" = 0,
    "Ca" = 0,
    "Cl" = 0.002,
    "Mg" = 0.001
  ),
  list(
    # "H" = 110.124,
    # "O" = 55.0622,
    # "Charge" = -1.216307659761E-09,
    # "C(4)" = 1.230067028174E-04,
    # "Ca" = 1.230067028174E-04,
    # "Cl" = 0,
    # "Mg" = 0
    "H" = 110.124,
    "O" = 55.0622,
    "Charge" = -1.217e-09,
    "C(4)" = 0,
    "Ca" = 0,
    "Cl" = 0,
    "Mg" = 0
  )
)

vecinj_inner <- list(
  # l1 = c(2, 1, 1)
)

# Create a list to store grid cell information
flux_list <- list()

# Function to get the index of a grid cell given its row and column
get_index <- function(row, col) {
  index <- (row - 1) * m + col
  if (index < 1) {
    index <- -1
  } else if (index > n * m) {
    index <- -1
  }
  return(index)
}

flux_val <- 0.005
# Loop through each row and column to populate the flux_list
for (row in 1:n) {
  for (col in 1:m) {
    index <- get_index(row, col)

    # Initialize the connections for the current cell
    flux <- c()

    # Check and add connections to the east, south, west, and north cells
    # east
    flux <- c(flux, -flux_val)

    # south
    flux <- c(flux, -flux_val)

    # west
    flux <- c(flux, flux_val)

    # north
    flux <- c(flux, flux_val)

    # Store the connections in the flux_list
    flux_list[[index]] <- flux
  }
}

vecinj <- do.call(rbind.data.frame, vecinj_adv)
names(vecinj) <- names(init_adv)

advection <- list(
  init = init_adv,
  vecinj = vecinj,
  vecinj_inner = vecinj_inner,
  const_flux = flux_list
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
  neg_sign <- (result["Calcite"] <- 0) || (result["Dolomite"] < 0)
  return(any(neg_sign))
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
)

#################################################################
##                          Section 4                          ##
##              Putting all those things together              ##
#################################################################


iterations <- 500
dt <- 600

out_save <- c(1, iterations)

setup <- list(
  grid = grid,
  advection = advection,
  chemistry = chemistry,
  iterations = iterations,
  timesteps = rep(dt, iterations),
  store_result = TRUE,
  out_save = out_save
)
