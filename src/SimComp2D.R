## chemical database
db <- RPhreeFile("mdl_quint_kin.dat", is.db=TRUE)
phreeqc::phrLoadDatabaseString(db)

## only the directory
demodir <- system.file("extdata", "demo_rtwithmufits", package="Rmufits")

prop <- c("Al", "C","Ca","Cl","Fe", "K", "Mg","Na", "Si", "pH", ## "pe",
          "Albite", "Calcite", "Chlorite", "Illite", "Kaolinite")

signif_vector <- c(7,7,7,7,7,7,7,7,7,6, 5,5,5,5,5)
prop_type <- rep("normal", length(signif_vector))

base <- c("SOLUTION 1",
          "units  mol/kgw",
          "pH     6.77",
          "temp  35",
          "-water 1",
          "Al  8.06386e-09",
          "C   0.0006108294",
          "Ca  0.09709463",
          "Cl  4.340042",
          "Fe  1.234357e-05",
          "K   0.01117434",
          "Mg  0.0406959",
          "Na  4.189209",
          "Si  0.0001935754",
          "INCREMENTAL_REACTIONS true",
          "KINETICS 1 ",
          "-steps 86400",
          "-bad_step_max 10000",
          "-cvode  true",
          "Albite",
          "-m      8.432165", ## 1540.0",
          "-parms  01.54  100",
          "Calcite",
          "-m   0.0",
          "-parms  10  100",
          "Chlorite",
          "-m      1.106585", ## 202.100",
          "-parms  64.84  100",
          "Illite",
          "-m      0.9549153", ## 174.400",
          "-parms  43.38  100",
          "Kaolinite",
          "-m      0.0",
          "-parms  29.17  100",
          "END")

selout <- c("KNOBS",
            "-convergence_tolerance 1E-6",
            "SELECTED_OUTPUT",
            "-reset            false",
            "USER_PUNCH",
            "-head  Al C Ca Cl Fe K Mg Na Si pH Albite Calcite Chlorite Illite Kaolinite", ## pe
            "10 PUNCH TOT(\"Al\"), TOT(\"C\"), TOT(\"Ca\"), TOT(\"Cl\"), TOT(\"Fe\"), TOT(\"K\"), TOT(\"Mg\"), TOT(\"Na\"), TOT(\"Si\"), -LA(\"H+\"), KIN(\"Albite\"), KIN(\"Calcite\"), KIN(\"Chlorite\"), KIN(\"Illite\"), KIN(\"Kaolinite\")" )

## Define initial conditions as equilibrium with primary minerals
ipr <- c(Al  = 8.689e-10,
         C   = 0.0006108,
         Ca  = 0.09709,
         Cl  = 4.34,
         Fe  = 1.802e-06,
         K   = 0.01131,
         Mg  = 0.04074,
         Na  = 4.189,
         Si  = 7.653e-05,
         pH  = 6.889,
         Albite    =  5.0,
         Calcite   =  0.0,
         Chlorite  = 10.0,
         Illite    =  2.0,
         Kaolinite =  0.0
         )

initstate <- matrix(rep(ipr, 2500), byrow=TRUE, ncol=length(ipr))
colnames(initstate) <- names(ipr)

vecinj <- c(Al= 8.694e-10,
            C = 8.182e-01,
            Ca= 9.710e-02,
            Cl= 4.340e+00,
            Fe= 1.778e-06,
            K = 1.131e-02,
            Mg= 4.074e-02,
            Na= 4.189e+00,
            Si= 7.652e-05,
            pH= 2.556228)


## setup boundary conditions for transport - we have already read the
## GRID with the following code:
## grid <- Rmufits::ReadGrid(paste0(demodir,"/d2ascii.run.GRID.SUM"))
## cbound <- which(grid$cell$ACTNUM == 2)
## dput(cbound)
cbound <- c(1L, 50L, 100L, 150L, 200L, 250L, 300L, 350L, 400L, 450L, 500L,
            550L, 600L, 650L, 700L, 750L, 800L, 850L, 900L, 950L, 1000L,
            1050L, 1100L, 1150L, 1200L, 1250L, 1300L, 1350L, 1400L, 1450L,
            1500L, 1550L, 1600L, 1650L, 1700L, 1750L, 1800L, 1850L, 1900L,
            1950L, 2000L, 2050L, 2100L, 2150L, 2200L, 2250L, 2300L, 2350L,
            2400L, 2450L, 2451L, 2452L, 2453L, 2454L, 2455L, 2456L, 2457L,
            2458L, 2459L, 2460L, 2461L, 2462L, 2463L, 2464L, 2465L, 2466L,
            2467L, 2468L, 2469L, 2470L, 2471L, 2472L, 2473L, 2474L, 2475L,
            2476L, 2477L, 2478L, 2479L, 2480L, 2481L, 2482L, 2483L, 2484L,
            2485L, 2486L, 2487L, 2488L, 2489L, 2490L, 2491L, 2492L, 2493L,
            2494L, 2495L, 2496L, 2497L, 2498L, 2499L, 2500L)

boundary_matrix <- matrix(rep(ipr[1:10], length(cbound)), byrow=TRUE, nrow=length(cbound))
colnames(boundary_matrix) <- names(ipr[1:10])
boundary_matrix[1, ] <- vecinj

boundary_matrix <- cbind(cbound,boundary_matrix)

setup <- list(n = 2500,
              bound = boundary_matrix,
              base  = base,
              first = selout,
              initsim = initstate,
              Cf    = 1,
              prop  = prop,
              immobile = seq(11,15),
              kin      = seq(11,15),
              phase   = "FLUX1",
              density = "DEN1",
              reduce = FALSE,
              snapshots = demodir, ## directory where we will read MUFITS SUM files
              gridfile  = paste0(demodir,"/d2ascii.run.GRID.SUM")
              )
