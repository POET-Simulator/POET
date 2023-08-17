## Time-stamp: "Last modified 2023-08-17 16:28:57 mluebke"

database <- normalizePath("../share/poet/bench/dolo/phreeqc_kin.dat")
input_script <- normalizePath("../share/poet/bench/dolo/dolo_inner.pqi")

#################################################################
##                          Section 1                          ##
##                     Grid initialization                     ##
#################################################################

n <- 50
m <- 50

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
  s_cells = c(1, 1),
  type = types[1]
)

##################################################################
##                          Section 2                           ##
##                     Advection parameters                     ##
##################################################################

## initial conditions
init_adv <- c(
  "H" = 110.683,
  "O" = 55.3413,
  "Charge" = -5.0822e-19,
  "C(4)" = 1.2279E-4,
  "Ca" = 1.2279E-4,
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
  )
)

vecinj_inner <- list(
  l1 = c(1, 1, 1)
)

# Create a list to store grid cell information
grid_list <- vector("list", n * m)

# Function to get the index of a grid cell given its row and column
get_index <- function(row, col) {
  return((row - 1) * m + col)
}

# Loop through each row and column to populate the grid_list
for (row in 1:n) {
  for (col in 1:m) {
    index <- get_index(row, col)

    # Initialize the connections for the current cell
    connections <- c()

    # Check and add connections to the east, south, west, and north cells
    if (col < m) {
      connections <- c(connections, get_index(row, col + 1))
    }
    if (row < n) {
      connections <- c(connections, get_index(row + 1, col))
    }
    if (col > index) {
      connections <- c(connections, get_index(row, col - 1))
    }
    if (row > index) {
      connections <- c(connections, get_index(row - 1, col))
    }

    # Store the connections in the grid_list
    grid_list[[index]] <- connections
  }
}

vecinj <- do.call(rbind.data.frame, vecinj_adv)
names(vecinj) <- names(init_adv)

advection <- list(
  init = init_adv,
  vecinj = vecinj,
  vecinj_inner = vecinj_inner,
  grid = grid_list
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
  "C(4)" = 5,
  "Ca" = 5,
  "Cl" = 5,
  "Mg" = 5,
  "Calcite" = 5,
  "Dolomite" = 5
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


iterations <- 10
dt <- 200

setup <- list(
  grid = grid,
  advection = advection,
  chemistry = chemistry,
  iterations = iterations,
  timesteps = rep(dt, iterations),
  store_result = TRUE
)
