database <- normalizePath("../data/phreeqc_kin.dat")
input_script <- normalizePath("../data/dol.pqi")

#################################################################
##                          Section 1                          ##
##                     Grid initialization                     ##
#################################################################

n <- 5
m <- 5

types <- c("scratch", "phreeqc", "rds")

# initsim <- c("SELECTED_OUTPUT", "-high_precision true",
#              "-reset false",
#              "USER_PUNCH",
#              "-head C Ca Cl Mg pH pe O2g Calcite Dolomite",
#              "10 PUNCH TOT(\"C\"), TOT(\"Ca\"), TOT(\"Cl\"), TOT(\"Mg\"), -LA(\"H+\"), -LA(\"e-\"), EQUI(\"O2g\"), EQUI(\"Calcite\"), EQUI(\"Dolomite\")",
#              "SOLUTION 1",
#              "units mol/kgw",
#              "temp 25.0", "water 1",
#              "pH 9.91 charge",
#              "pe 4.0",
#              "C   1.2279E-04",
#              "Ca  1.2279E-04",
#              "Cl 1E-12",
#              "Mg 1E-12",
#              "PURE 1",
#              "O2g     -0.6788 10.0",
#              "Calcite  0.0 2.07E-3",
#              "Dolomite 0.0 0.0",
#              "END")

# needed if init type is set to "scratch"
# prop <- c("C", "Ca", "Cl", "Mg", "pH", "pe", "O2g", "Calcite", "Dolomite")

init_cell <- list(
  "H" = 110.683,
  "O" = 55.3413,
  "Charge" = -5.0822e-19,
  "C" = 1.2279E-4,
  "Ca" = 1.2279E-4,
  "Cl" = 0,
  "Mg" = 0,
  "O2g" = 0.499957,
  "Calcite" = 2.07e-4,
  "Dolomite" = 0
)

grid <- list(
  n_cells = c(n, m),
  s_cells = c(n, m),
  type = types[1],
  init_cell = as.data.frame(init_cell),
  props = names(init_cell),
  database = database,
  input_script = input_script
)


##################################################################
##                          Section 2                           ##
##         Diffusion parameters and boundary conditions         ##
##################################################################

init_diffu <- c(
  "H" = 110.683,
  "O" = 55.3413,
  "Charge" = -5.0822e-19,
  "C" = 1.2279E-4,
  "Ca" = 1.2279E-4,
  "Cl" = 0,
  "Mg" = 0
)

alpha_diffu <- c(
  "H" = 1E-4,
  "O" = 1E-4,
  "Charge" = 1E-4,
  "C" = 1E-4,
  "Ca" = 1E-4,
  "Cl" = 1E-4,
  "Mg" = 1E-4
)

vecinj_diffu <- list(
  list(
    "H" = 110.683,
    "O" = 55.3413,
    "Charge" = 1.90431e-16,
    "C" = 0,
    "Ca" = 0,
    "Cl" = 0.002,
    "Mg" = 0.001
  )
)

#inner_index <- c(5, 15, 25)
#inner_vecinj_index <- rep(1, 3)
#
#vecinj_inner <- cbind(inner_index, inner_vecinj_index)
vecinj_inner <- list(
  l1 = c(1,2,2)
)


boundary <- list(
  "N" = rep(0, n),
  "E" = rep(0, n),
  "S" = rep(0, n),
  "W" = rep(0, n)
)

diffu_list <- names(alpha_diffu)

diffusion <- list(
  init = init_diffu,
  vecinj = do.call(rbind.data.frame, vecinj_diffu),
  vecinj_inner = vecinj_inner,
  vecinj_index = boundary,
  alpha = alpha_diffu
)

#################################################################
##                          Section 3                          ##
##                  Chemitry module (Phreeqc)                  ##
#################################################################

# db <- RPhreeFile(system.file("extdata", "phreeqc_kin.dat",
#  package = "RedModRphree"
# ), is.db = TRUE)
#
# phreeqc::phrLoadDatabaseString(db)

# NOTE: This won't be needed in the future either. Could also be done in a. pqi
# file
base <- c(
  "SOLUTION 1",
  "units mol/kgw",
  "temp 25.0",
  "water 1",
  "pH 9.91 charge",
  "pe 4.0",
  "C   1.2279E-04",
  "Ca  1.2279E-04",
  "Mg 0.001",
  "Cl 0.002",
  "PURE 1",
  "O2g -0.1675 10",
  "KINETICS 1",
  "-steps 100",
  "-step_divide 100",
  "-bad_step_max 2000",
  "Calcite", "-m      0.000207",
  "-parms  0.0032",
  "Dolomite",
  "-m      0.0",
  "-parms  0.00032",
  "END"
)

selout <- c(
  "SELECTED_OUTPUT", "-high_precision true", "-reset false",
  "-time", "-soln", "-temperature true", "-water true",
  "-pH", "-pe", "-totals C Ca Cl Mg",
  "-kinetic_reactants Calcite Dolomite", "-equilibrium O2g"
)


# Needed when using DHT
signif_vector <- c(7, 7, 7, 7, 7, 7, 7, 5, 5, 7)
prop_type <- c("act", "act", "act", "act", "logact", "logact", "ignore", "act", "act", "act")
prop <- names(init_cell)

chemistry <- list(
  database = database,
  input_script = input_script
)

#################################################################
##                          Section 4                          ##
##              Putting all those things together              ##
#################################################################


iterations <- 10

setup <- list(
  # bound = myboundmat,
  base = base,
  first = selout,
  # initsim = initsim,
  # Cf = 1,
  grid = grid,
  diffusion = diffusion,
  chemistry = chemistry,
  prop = prop,
  immobile = c(7, 8, 9),
  kin = c(8, 9),
  ann = list(O2g = -0.1675),
  # phase = "FLUX1",
  # density = "DEN1",
  reduce = FALSE,
  # snapshots = demodir, ## directory where we will read MUFITS SUM files
  # gridfile = paste0(demodir, "/d2ascii.run.GRID.SUM")
  # init = init,
  # vecinj = vecinj,
  # cinj = c(0,1),
  # boundary = boundary,
  # injections = FALSE,
  iterations = iterations,
  timesteps = rep(10, iterations)
)
