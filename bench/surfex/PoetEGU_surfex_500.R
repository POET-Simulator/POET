rows <- 500
cols <- 200

grid_left <- matrix(1, nrow = rows, ncol = cols/2)
grid_rght <- matrix(2, nrow = rows, ncol = cols/2)
grid_def <- cbind(grid_left, grid_rght)


# Define grid configuration for POET model
grid_setup <- list(
  pqc_in_file = "./SurfexEGU.pqi",
  pqc_db_file = "./SMILE_2021_11_01_TH.dat", # Path to the database file for Phreeqc
  grid_def = grid_def,    # Definition of the grid, containing IDs according to the Phreeqc input script
  grid_size = c(10, 4),   # Size of the grid in meters
  constant_cells = c()    # IDs of cells with constant concentration
)

bound_def <- list(
  "type" = rep("constant", cols),
  "sol_id" = rep(3, cols),
  "cell" = seq(1, cols)
)

diffusion_setup <- list(
  boundaries = list(
    "N" = bound_def
  ),
  alpha_x = matrix(runif(rows*cols))*1e-8, 
  alpha_y = matrix(runif(rows*cols))*1e-9## ,1e-10
)


chemistry_setup <- list()

# Define a setup list for simulation configuration
setup <- list(
  Grid = grid_setup, # Parameters related to the grid structure
  Diffusion = diffusion_setup, # Parameters related to the diffusion process
  Chemistry = chemistry_setup # Parameters related to the chemistry process
)
