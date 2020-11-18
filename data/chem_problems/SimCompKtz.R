db <- RPhreeFile("mdl_quint_kin.dat", is.db=TRUE)
phreeqc::phrLoadDatabaseString(db)

## only the directory
demodir <- "./snaps/"
         
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

initstate <- matrix(rep(ipr, 648420), byrow=TRUE, ncol=length(ipr))
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

            
            

prop <- c("Al", "C","Ca","Cl","Fe", "K", "Mg","Na", "Si", "pH", ## "pe",
          "Albite", "Calcite", "Chlorite", "Illite", "Kaolinite")

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
              base=base, 
              bound=myboundmat, 
              first=selout, 
              initsim=initstate,  
              Cf=1, 
              prop=prop, 
              immobile=seq(11,15),
              kin= seq(11,15),
              phase="FLUX1",
              density="DENS",
              reduce=TRUE,
              snapshots="snaps/AllSnaps_cmp_v3.rds",
              gridfile ="snaps/GridKtz_cmp_v3.rds")

         
