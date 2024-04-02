grid_def <- matrix(2, nrow = 200, ncol = 200)

# Define grid configuration for POET model
grid_setup <- list(
    pqc_in_file = "./dol.pqi",
    pqc_db_file = "./phreeqc_kin.dat", # Path to the database file for Phreeqc
    grid_def = grid_def, # Definition of the grid, containing IDs according to the Phreeqc input script
    grid_size = c(ncol(grid_def), nrow(grid_def)), # Size of the grid in meters
    constant_cells = c() # IDs of cells with constant concentration
)

bound_size <- 2

diffusion_setup <- list(
    boundaries = list(
        "W" = list(
            "type" = rep("constant", bound_size),
            "sol_id" = rep(3, bound_size),
            "cell" = seq(1, bound_size)
        ),
        "N" = list(
            "type" = rep("constant", bound_size),
            "sol_id" = rep(3, bound_size),
            "cell" = seq(1, bound_size)
        )
    ),
    alpha_x = 1e-6,
    alpha_y = 1e-6
)

# Define a setup list for simulation configuration
setup <- list(
    Grid = grid_setup, # Parameters related to the grid structure
    Diffusion = diffusion_setup # Parameters related to the diffusion process
)
