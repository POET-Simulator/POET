rows <- 100
cols <- 100

grid_def <- matrix(1, nrow = rows, ncol = cols)

# Define grid configuration for POET model
grid_setup <- list(
  pqc_in_file = "./SurfExBase.pqi",
  pqc_db_file = "./SMILE_2021_11_01_TH.dat", # Path to the database file for Phreeqc
  grid_def = grid_def, # Definition of the grid, containing IDs according to the Phreeqc input script
  grid_size = c(1, 1), # Size of the grid in meters
  constant_cells = c() # IDs of cells with constant concentration
)

bound_def <- list(
  "type" = rep("constant", cols),
  "sol_id" = rep(2, cols),
  "cell" = seq(1, cols)
)

diffusion_setup <- list(
  boundaries = list(
    "N" = bound_def
  ),
  alpha_x = 1e-6,
  alpha_y = 1e-6
)


chemistry_setup <- list()

# Define a setup list for simulation configuration
setup <- list(
  Grid = grid_setup, # Parameters related to the grid structure
  Diffusion = diffusion_setup, # Parameters related to the diffusion process
  Chemistry = chemistry_setup # Parameters related to the chemistry process
)