## Time-stamp: "Last modified 2024-12-20 15:30:05 delucia"

cols <- 400
rows <- 400

dim_cols <- 80
dim_rows <- 80

rad   <- cols/40  ## circle radius (in nodes)
boxs  <- cols/16  ## box side (in nodes)

cc <- seq(boxs/2, cols - boxs/2, by = boxs)
length(cc)

colcen <- rep(cc, length(cc))
rowcen <- rep(cc, each = length(cc))
centers <- cbind(rowcen, colcen)

mr <- matrix(rep(seq_len(cols), each = rows), byrow = TRUE, nrow = rows)
mc <- matrix(rep(seq_len(cols), rows), byrow = TRUE, nrow = rows)

tmpl <- lapply(seq_len(nrow(centers)), function(x) which((mr-centers[x, 1])^2 + (mc-centers[x, 2])^2 < rad^2, arr.ind = TRUE))

inds <-  do.call(rbind, tmpl)

grid <- matrix(1, nrow = rows, ncol = cols)
grid[inds] <- 2

alpha <- matrix( 1e-5, ncol = cols, nrow = rows)
alpha[inds] <- 1e-8

## PoetUtils::PlotField(grid, cols = cols, rows = rows, contour = FALSE, scale = FALSE, las=1)

## Define grid configuration for POET model
grid_setup <- list(
  pqc_in_file    = "./barite_fgcs_3.pqi",
  pqc_db_file    = "../barite/db_barite.dat",     ## database file
  grid_def       = grid,                  ## grid definition, IDs according to the Phreeqc input 
  grid_size      = c(dim_cols, dim_rows), ## grid size in meters
  constant_cells = c()                    ## IDs of cells with constant concentration
)

bound_length <- cols/10

bound_N <- list(
  "type"   = rep("constant", bound_length),
  "sol_id" = rep(3, bound_length),
  "cell"   = seq(1, bound_length)
)

bound_W <- list(
  "type"   = rep("constant", bound_length),
  "sol_id" = rep(3, bound_length),
  "cell"   = seq(1, bound_length)
)
bound_E <- list(
  "type"   = rep("constant", bound_length),
  "sol_id" = rep(4, bound_length),
  "cell"   = seq(rows-bound_length+1, rows)
)

bound_S <- list(
  "type"   = rep("constant", bound_length),
  "sol_id" = rep(4, bound_length),
  "cell"   = seq(cols-bound_length+1, cols)
)

diffusion_setup <- list(
  boundaries = list(
    "W" = bound_W,
    "N" = bound_N,
    "E" = bound_E,
    "S" = bound_S
  ),
  alpha_x = alpha,
  alpha_y = alpha
)

dht_species <- c(
  "H"         = 7,
  "O"         = 7,
  "Ba"        = 7,
  "Cl"        = 7,
  "S"         = 7,
  "Sr"        = 7,
  "Barite"    = 4,
  "Celestite" = 4
)

pht_species <- c(
  "Ba"        = 3,
  "Cl"        = 3,
  "S"         = 3,
  "Sr"        = 3,
  "Barite"    = 1,
  "Celestite" = 1
)

chemistry_setup <- list(
  dht_species = dht_species,
  pht_species = pht_species
)

## Define a setup list for simulation configuration
setup <- list(
  Grid = grid_setup, ## Parameters related to the grid structure
  Diffusion = diffusion_setup, ## Parameters related to the diffusion process
  Chemistry = chemistry_setup
)
