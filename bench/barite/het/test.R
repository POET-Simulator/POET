grid_def <- matrix(c(2, 3), nrow = 2, ncol = 5)

# Define grid configuration for POET model
grid_setup <- list(
    pqc_in = "./barite_het.pqi",
    pqc_db = "./db_barite.dat", # Path to the database file for Phreeqc
    grid_def = grid_def, # Definition of the grid, containing IDs according to the Phreeqc input script
    grid_size = c(ncol(grid_def), nrow(grid_def)), # Size of the grid in meters
    constant_cells = c() # IDs of cells with constant concentration
)

# Define a setup list for simulation configuration
setup <- list(
    grid = grid_setup # Parameters related to the grid structure
)
