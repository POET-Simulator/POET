## Time-stamp: "Last modified 2024-05-30 13:34:14 delucia"
cols <- 50
rows <- 50

s_cols <- 0.25
s_rows <- 0.25

grid_def <- matrix(2, nrow = rows, ncol = cols)

# Define grid configuration for POET model
grid_setup <- list(
  pqc_in_file = "./barite.pqi",
  pqc_db_file = "./db_barite.dat", ## Path to the database file for Phreeqc
  grid_def = grid_def, ## Definition of the grid, containing IDs according to the Phreeqc input script
  grid_size = c(s_rows, s_cols), ## Size of the grid in meters
  constant_cells = c() ## IDs of cells with constant concentration
)

bound_length <- 2

bound_def <- list(
  "type" = rep("constant", bound_length),
  "sol_id" = rep(3, bound_length),
  "cell" = seq(1, bound_length)
)

homogenous_alpha <- 1e-8

diffusion_setup <- list(
  boundaries = list(
    "W" = bound_def,
    "N" = bound_def
  ),
  alpha_x = homogenous_alpha,
  alpha_y = homogenous_alpha
)

dht_species <- c(
  "H"         = 3,
  "O"         = 3,
  "Charge"    = 3,
  "Ba"        = 6,
  "Cl"        = 6,
  "S"         = 6,
  "Sr"        = 6,
  "Barite"    = 5,
  "Celestite" = 5
)

chemistry_setup <- list(
  dht_species = dht_species,
  ai_surrogate_input_script = "./barite_50ai_surr_mdl.R"
)

# Define a setup list for simulation configuration
setup <- list(
  Grid = grid_setup, # Parameters related to the grid structure
  Diffusion = diffusion_setup, # Parameters related to the diffusion process
  Chemistry = chemistry_setup
)
