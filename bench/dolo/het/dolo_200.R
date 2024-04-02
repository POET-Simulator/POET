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
    return(input[names(dht_species)])
}

chemistry_setup <- list(
    dht_species = c(
        "H" = 3,
        "O" = 3,
        "Charge" = 3,
        "C(4)" = 6,
        "Ca" = 6,
        "Cl" = 3,
        "Mg" = 5,
        "Calcite" = 4,
        "Dolomite" = 4
    ),
    hooks = list(
        dht_fill = check_sign_cal_dol_dht,
        dht_fuzz = fuzz_input_dht_keys
    )
)

# Define a setup list for simulation configuration
setup <- list(
    Grid = grid_setup, # Parameters related to the grid structure
    Diffusion = diffusion_setup, # Parameters related to the diffusion process
    Chemistry = chemistry_setup # Parameters related to the chemistry process
)
