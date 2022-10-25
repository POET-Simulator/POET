## chemical database
db <- RPhreeFile(system.file("extdata", "phreeqc_kin.dat",
                             package="RedModRphree"), is.db=TRUE)

phreeqc::phrLoadDatabaseString(db)

## only the directory
demodir <- system.file("extdata", "demo_rtwithmufits", package="Rmufits")

prop <- c("C","Ca","Cl","Mg","pH","pe","O2g", "Calcite","Dolomite")

signif_vector <- c(7,7,7,7,7,7,7,5,5)
prop_type <- c("act","act","act","act","logact","logact","ignore","act","act")


base <- c("SOLUTION 1",
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
          "END")

selout <- c("SELECTED_OUTPUT", "-high_precision true", "-reset false",
            "-time", "-soln", "-temperature true", "-water true",
            "-pH", "-pe", "-totals C Ca Cl Mg",
            "-kinetic_reactants Calcite Dolomite", "-equilibrium O2g")

initsim <- c("SELECTED_OUTPUT", "-high_precision true",
             "-reset false",
             "USER_PUNCH",
             "-head C Ca Cl Mg pH pe O2g Calcite Dolomite",
             "10 PUNCH TOT(\"C\"), TOT(\"Ca\"), TOT(\"Cl\"), TOT(\"Mg\"), -LA(\"H+\"), -LA(\"e-\"), EQUI(\"O2g\"), EQUI(\"Calcite\"), EQUI(\"Dolomite\")",
             "SOLUTION 1",
             "units mol/kgw",
             "temp 25.0", "water 1",
             "pH 9.91 charge",
             "pe 4.0",
             "C   1.2279E-04",
             "Ca  1.2279E-04",
             "Cl 1E-12",
             "Mg 1E-12",
             "PURE 1",
             "O2g     -0.6788 10.0",
             "Calcite  0.0 2.07E-3",
             "Dolomite 0.0 0.0",
             "END")

vecinj <- c("C"= 0,
            "Ca"  = 0,
            "Cl"  = 0.002,
            "Mg"  = 0.001,
            "pe"   = 4,
            "pH"   = 7)

init <- c("C(4)"= 1.2279E-4,
          "Ca"  =1.2279E-4,
          "Cl"  =0,        
          "Mg"  =0,
          "pe"  =4,
          "pH"  =7,
          "Calcite"= 2.07e-4,
          "Dolomite"= 0)


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

boundinit <- matrix(rep(init[-c(7,8)], length(cbound)), byrow=TRUE, nrow=length(cbound))
myboundmat <- cbind(cbound,boundinit)
myboundmat[cbound==1, c(2:7)] <- vecinj
colnames(myboundmat) <- c("cbound", names(vecinj))

# TODO: dt and iterations

setup <- list(n=2500,
              bound=myboundmat,   
              base=base,
              first=selout,
              initsim=initsim,
              Cf=1,
              prop=prop,
              immobile=c(7,8,9),
              kin= c(8,9),
              ann=list(O2g=-0.1675), 
              phase="FLUX1",
              density="DEN1",
              reduce=FALSE,
              snapshots=demodir, ## directory where we will read MUFITS SUM files
              gridfile=paste0(demodir,"/d2ascii.run.GRID.SUM")
              )
