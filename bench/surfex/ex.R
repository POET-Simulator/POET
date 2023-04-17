## Time-stamp: "Last modified 2023-04-17 12:29:27 mluebke"

database <- normalizePath("./SMILE_2021_11_01_TH.dat")
input_script <- normalizePath("./ExBase.pqi")

cat(paste(":: R This is a test 1\n"))

#################################################################
##                          Section 1                          ##
##                     Grid initialization                     ##
#################################################################

n <- 100
m <- 100

types <- c("scratch", "phreeqc", "rds")

init_cell <- list(H = 1.476571028625e-01,
                  O = 7.392297218936e-02,
                  Charge = -1.765225732724e-18,
                  `C(-4)` = 2.477908970828e-21,
                  `C(4)` = 2.647623016916e-06,
                  Ca = 2.889623169138e-05,
                  Cl = 4.292806181039e-04,
                  `Fe(2)` =1.908142472666e-07,
                  `Fe(3)` =3.173306589931e-12,
                  `H(0)` =2.675642675119e-15,
                  K = 2.530134809667e-06,
                  Mg =2.313806319294e-05,
                  Na =3.674633059628e-04,
                  `S(-2)` = 8.589766637180e-15,
                  `S(2)` = 1.205284362720e-19,
                  `S(4)` = 9.108958772790e-18,
                  `S(6)` = 2.198092329098e-05,
                  Sr = 6.012080128154e-07,
                  `U(4)` = 1.039668623852e-14,
                  `U(5)` = 1.208394829796e-15,
                  `U(6)` = 2.976409147150e-12)

grid <- list(
  n_cells = c(n, m),
  s_cells = c(1, 1),
  type = "scratch",
  init_cell = as.data.frame(init_cell, check.names = FALSE),
  props = names(init_cell),
  database = database,
  input_script = input_script
)


##################################################################
##                          Section 2                           ##
##         Diffusion parameters and boundary conditions         ##
##################################################################

vecinj_diffu <- list(
    list(H = 0.147659686316291,
         O = 0.0739242798146046,
         Charge = 7.46361643222701e-20,
         `C(-4)` = 2.92438561098248e-21,
         `C(4)` = 2.65160558871092e-06,
         Ca = 2.89001071336443e-05,
         Cl = 0.000429291158114428,
         `Fe(2)` = 1.90823391198114e-07,
         `Fe(3)` = 3.10832423034763e-12,
         `H(0)` = 2.7888235127385e-15,
         K = 2.5301787e-06,
         Mg = 2.31391999937907e-05,
         Na = 0.00036746969,
         `S(-2)` = 1.01376078438546e-14,
         `S(2)` = 1.42247026981542e-19,
         `S(4)` = 9.49422092568557e-18,
         `S(6)` = 2.19812504654191e-05,
         Sr = 6.01218519999999e-07,
         `U(4)` = 4.82255946569383e-12,
         `U(5)` = 5.49050615347901e-13,
         `U(6)` = 1.32462838991902e-09)
)

vecinj <- do.call(rbind.data.frame, vecinj_diffu)
names(vecinj) <- grid$props

## diffusion coefficients
alpha_diffu <- c(H = 1E-6, O = 1E-6, Charge = 1E-6, `C(-4)` = 1E-6,
                 `C(4)` = 1E-6, Ca = 1E-6, Cl = 1E-6, `Fe(2)` = 1E-6,
                 `Fe(3)` = 1E-6, `H(0)` = 1E-6, K = 1E-6, Mg = 1E-6,
                 Na = 1E-6, `S(-2)` = 1E-6, `S(2)` = 1E-6,
                 `S(4)` = 1E-6, `S(6)` = 1E-6, Sr = 1E-6,
                 `U(4)` = 1E-6, `U(5)` = 1E-6, `U(6)` = 1E-6)

## list of boundary conditions/inner nodes

## vecinj_inner <- list(
##    list(1,1,1)
## )

boundary <- list(
  "N" = rep(1, n),
  "E" = rep(0, n),
  "S" = rep(0, n),
  "W" = rep(0, n)
)

diffu_list <- names(alpha_diffu)

vecinj <- do.call(rbind.data.frame, vecinj_diffu)
names(vecinj) <- names(init_cell)

diffusion <- list(
  init = as.data.frame(init_cell, check.names = FALSE),
  vecinj = vecinj,
#  vecinj_inner = vecinj_inner,
  vecinj_index = boundary,
  alpha = alpha_diffu
)

#################################################################
##                          Section 3                          ##
##                  Chemistry module (Phreeqc)                 ##
#################################################################


chemistry <- list(
  database = database,
  input_script = input_script
)

#################################################################
##                          Section 4                          ##
##              Putting all those things together              ##
#################################################################


iterations <- 10
dt <- 200

setup <- list(
  grid = grid,
  diffusion = diffusion,
  chemistry = chemistry,
  iterations = iterations,
  timesteps = rep(dt, iterations),
  store_result = TRUE
)
