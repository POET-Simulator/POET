# library(RedModRphree)
# library(Rmufits)
# library(RcppVTK)

db <- RPhreeFile(system.file("extdata", "phreeqc_kin.dat",
                             package="RedModRphree"), is.db=TRUE)

phreeqc::phrLoadDatabaseString(db)

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

selout <- c("SELECTED_OUTPUT",
            "-high_precision true",
            "-reset false",
            "-time",
            "-soln",
            "-temperature true",
            "-water true",
            "-pH",
            "-pe",
            "-totals C Ca Cl Mg",
            "-kinetic_reactants Calcite Dolomite",
            "-equilibrium O2g")

initsim <- c("SELECTED_OUTPUT",
             "-high_precision true",
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

bound_elm <- c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 13L, 14L, 
               15L, 16L, 17L, 18L, 19L, 20L, 21L, 22L, 23L, 24L, 25L, 26L, 27L, 
               28L, 29L, 30L, 31L, 32L, 33L, 34L, 35L, 36L, 37L, 38L, 39L, 40L, 
               41L, 42L, 43L, 44L, 45L, 46L, 47L, 48L, 49L, 50L, 51L, 52L, 53L, 
               54L, 55L, 56L, 57L, 58L, 59L, 60L, 61L, 62L, 63L, 64L, 65L, 66L, 
               67L, 68L, 69L, 70L, 71L, 72L, 73L, 74L, 75L, 76L, 77L, 78L, 79L, 
               80L, 81L, 82L, 83L, 84L, 85L, 86L, 87L, 88L, 89L, 90L, 91L, 92L, 
               93L, 94L, 95L, 96L, 97L, 98L, 99L, 100L, 101L, 102L, 103L, 104L, 
               105L, 106L, 107L, 108L, 214L, 215L)

inj_elm <- c(7426L, 18233L, 29040L, 39847L, 
             50654L, 61461L, 72268L, 83075L, 93882L, 104689L, 115496L, 126303L, 
             137110L, 147917L, 158724L, 169531L, 180338L, 191145L, 201952L, 
             212759L, 223566L, 234373L, 245180L, 255987L, 266794L, 277601L, 
             288408L, 299215L, 310022L, 320829L, 331636L, 342443L, 353250L, 
             364057L, 374864L, 385671L, 396478L, 407285L, 418092L, 428899L, 
             439706L, 450513L, 461320L, 472127L, 482934L, 493741L, 504548L, 
             515355L)

cbound <- inj_elm

boundinit <- matrix(rep(vecinj, length(cbound)), ncol=length(vecinj), byrow=TRUE) 
myboundmat <- cbind(cbound,boundinit)

## distinguish between injection and real boundaries
colnames(myboundmat) <- c("cbound", names(vecinj))

setup <- list(n=648420,
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
              density="DENS",
              reduce=FALSE,
              snapshots="snaps/AllSnaps_cmp_v3.rds",
              gridfile ="snaps/GridKtz_cmp_v3.rds"
              )
