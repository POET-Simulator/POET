## Time-stamp: "Last modified 2023-08-22 12:44:53 mluebke"

#################################################################
##                          Section 1                          ##
##                     Grid initialization                     ##
#################################################################

n <- 5
m <- 5

grid <- list(
  n_cells = c(n, m),
  s_cells = c(n, m)
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
    "H" = 110.124,
    "O" = 55.0622,
    "Charge" = -1.217e-09,
    "C(4)" = 0,
    "Ca" = 0,
    "Cl" = 0,
    "Mg" = 0
  ),
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
  l1 = c(2, 1, 1)
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

flux_val <- 1

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
##                          Section 4                          ##
##              Putting all those things together              ##
#################################################################

setup <- list(
  grid = grid,
  advection = advection
)
