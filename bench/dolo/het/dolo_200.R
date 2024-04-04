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

check_sign_cal_dol_interp <- function(to_interp, data_set) {
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

# Optional when using Interpolation (example with less key species and custom
# significant digits)

pht_species <- c(
    "C(4)" = 3,
    "Ca" = 3,
    "Mg" = 2,
    "Calcite" = 2,
    "Dolomite" = 2
)

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
    pht_species = pht_species,
    hooks = list(
        dht_fill = check_sign_cal_dol_dht,
        dht_fuzz = fuzz_input_dht_keys,
        interp_pre = check_sign_cal_dol_interp,
        interp_post = check_neg_cal_dol
    )
)

# Define a setup list for simulation configuration
setup <- list(
    Grid = grid_setup, # Parameters related to the grid structure
    Diffusion = diffusion_setup, # Parameters related to the diffusion process
    Chemistry = chemistry_setup # Parameters related to the chemistry process
)
