grid_def <- matrix(c(2, 3), nrow = 2, ncol = 5)

# Define grid configuration for POET model
grid_setup <- list(
    pqc_in_file = "./barite_het.pqi",
    pqc_db_file = "./db_barite.dat", # Path to the database file for Phreeqc
    grid_def = grid_def, # Definition of the grid, containing IDs according to the Phreeqc input script
    grid_size = c(ncol(grid_def), nrow(grid_def)), # Size of the grid in meters
    constant_cells = c() # IDs of cells with constant concentration
)

diffusion_setup <- list(
    boundaries = list(
        "W" = list(
            "type" = rep("constant", nrow(grid_def)),
            "sol_id" = rep(4, nrow(grid_def)),
            "cell" = seq_len(nrow(grid_def))
        )
    ),
    alpha_x = 1e-6,
    alpha_y = matrix(runif(10, 1e-8, 1e-7),
        nrow = nrow(grid_def),
        ncol = ncol(grid_def)
    )
)

# Define a setup list for simulation configuration
setup <- list(
    Grid = grid_setup, # Parameters related to the grid structure
    Diffusion = diffusion_setup, # Parameters related to the diffusion process
    Chemistry = list()
)
