## Time-stamp: "Last modified 2023-08-15 11:58:23 delucia"

### Copyright (C) 2018-2023 Marco De Lucia, Max Luebke (GFZ Potsdam)
###
### POET is free software; you can redistribute it and/or modify it under the
### terms of the GNU General Public License as published by the Free Software
### Foundation; either version 2 of the License, or (at your option) any later
### version.
###
### POET is distributed in the hope that it will be useful, but WITHOUT ANY
### WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
### A PARTICULAR PURPOSE. See the GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License along with
### this program; if not, write to the Free Software Foundation, Inc., 51
### Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

master_init <- function(setup, out_dir, init_field) {
    ## Setup the directory where we will store the results
    if (!dir.exists(out_dir)) {
        dir.create(out_dir)
        msgm("created directory ", out_dir)
    } else {
        msgm("dir ", out_dir, " already exists, I will overwrite!")
    }
    if (is.null(setup$store_result)) {
        msgm("store_result doesn't exist!")
    } else {
        msgm("store_result is ", setup$store_result)
    }

    setup$iter <- 1
    setup$timesteps <- setup$timesteps
    setup$maxiter <- length(setup$timesteps)
    setup$iterations <- setup$maxiter
    setup$simulation_time <- 0

    if (is.null(setup[["store_result"]])) {
        setup$store_result <- TRUE
    }

    if (setup$store_result) {
        init_field_out <- paste0(out_dir, "/iter_0.rds")
        init_field <- data.frame(init_field, check.names = FALSE)
        saveRDS(init_field, file = init_field_out)
        msgm("Stored initial field in ", init_field_out)
        if (is.null(setup[["out_save"]])) {
            setup$out_save <- seq(1, setup$iterations)
        }
    }

    setup$out_dir <- out_dir

    return(setup)
}

## This function, called only by master, stores on disk the last
## calculated time step if store_result is TRUE and increments the
## iteration counter
master_iteration_end <- function(setup, state_T, state_C) {
    iter <- setup$iter
    # print(iter)
    ## max digits for iterations
    dgts <- as.integer(ceiling(log10(setup$maxiter)))
    ## string format to use in sprintf
    fmt <- paste0("%0", dgts, "d")

    ## Write on disk state_T and state_C after every iteration
    ## comprised in setup$out_save
    if (setup$store_result) {
        if (iter %in% setup$out_save) {
            nameout <- paste0(setup$out_dir, "/iter_", sprintf(fmt = fmt, iter), ".rds")
            state_T <- data.frame(state_T, check.names = FALSE)
            state_C <- data.frame(state_C, check.names = FALSE)
            
            ai_surrogate_info <- list(
                prediction_time = if(exists("ai_prediction_time")) as.integer(ai_prediction_time) else NULL,
                training_time = if(exists("ai_training_time")) as.integer(ai_training_time) else NULL,
                valid_predictions = if(exists("validity_vector")) validity_vector else NULL)
            saveRDS(list(
                T = state_T,
                C = state_C,
                simtime = as.integer(setup$simulation_time),
                totaltime = as.integer(totaltime),
                ai_surrogate_info = ai_surrogate_info
            ), file = nameout)
            msgm("results stored in <", nameout, ">")
        }
    }
    ## Add last time step to simulation time
    setup$simulation_time <- setup$simulation_time + setup$timesteps[iter]

    msgm("done iteration", iter, "/", length(setup$timesteps))
    setup$iter <- setup$iter + 1
    return(setup)
}


## Attach the name of the calling function to the message displayed on
## R's stdout
msgm <- function(...) {
    prefix <- paste0("R: ")
    cat(paste(prefix, ..., "\n"))
    invisible()
}


## Function called by master R process to store on disk all relevant
## parameters for the simulation
StoreSetup <- function(setup, filesim, out_dir) {
    to_store <- vector(mode = "list", length = 4)
    ## names(to_store) <- c("Sim", "Flow", "Transport", "Chemistry", "DHT")
    names(to_store) <- c("Sim", "Transport", "DHT", "Cmdline")

    ## read the setup R file, which is sourced in kin.cpp
    tmpbuff <- file(filesim, "r")
    setupfile <- readLines(tmpbuff)
    close.connection(tmpbuff)

    to_store$Sim <- setupfile

    ## to_store$Flow <- list(
    ##     snapshots  = setup$snapshots,
    ##     gridfile   = setup$gridfile,
    ##     phase      = setup$phase,
    ##     density    = setup$density,
    ##     dt_differ  = setup$dt_differ,
    ##     prolong    = setup$prolong,
    ##     maxiter    = setup$maxiter,
    ##     saved_iter = setup$iter_output,
    ##     out_save   = setup$out_save )

    to_store$Transport <- setup$diffusion

    ## to_store$Chemistry <- list(
    ##    nprocs   = n_procs,
    ##    wp_size  = work_package_size,
    ##    base     = setup$base,
    ##    first    = setup$first,
    ##    init     = setup$initsim,
    ##    db       = db,
    ##    kin      = setup$kin,
    ##    ann      = setup$ann)

    if (dht_enabled) {
        to_store$DHT <- list(
            enabled   = dht_enabled,
            log       = dht_log
            ## signif    = dht_final_signif,
            ## proptype  = dht_final_proptype
        )
    } else {
        to_store$DHT <- FALSE
    }

    if (dht_enabled) {
        to_store$DHT <- list(
            enabled   = dht_enabled,
            log       = dht_log
            # signif    = dht_final_signif,
            # proptype  = dht_final_proptype
        )
    } else {
        to_store$DHT <- FALSE
    }

    saveRDS(to_store, file = paste0(fileout, "/setup.rds"))
    msgm("initialization stored in ", paste0(fileout, "/setup.rds"))
}

GetWorkPackageSizesVector <- function(n_packages, package_size, len) {
    ids <- rep(1:n_packages, times = package_size, each = 1)[1:len]
    return(as.integer(table(ids)))
}
