#ifndef POET_H
#define POET_H

#include <string>
static const char *poet_version = "tmp/ml/v0.2-107-g9990c12";

// using the Raw string literal to avoid escaping the quotes
static inline std::string kin_r_library = R"(## Time-stamp: "Last modified 2023-08-15 11:58:23 delucia"

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


## Simple function to check file extension. It is needed to check if
## the GridFile is SUM (MUFITS format) or rds/RData
FileExt <- function(x) {
    pos <- regexpr("\\.([[:alnum:]]+)$", x)
    ifelse(pos > -1L, substring(x, pos + 1L), "")
}

master_init <- function(setup) {
    msgm("Process with rank 0 reading GRID properties")

    ## Setup the directory where we will store the results
    verb <- FALSE
    if (local_rank == 0) {
        verb <- TRUE ## verbosity loading MUFITS results
        if (!dir.exists(fileout)) {
            dir.create(fileout)
            msgm("created directory ", fileout)
        } else {
            msgm("dir ", fileout, " already exists, I will overwrite!")
        }
        if (!exists("store_result")) {
            msgm("store_result doesn't exist!")
        } else {
            msgm("store_result is ", store_result)
        }
    } else {
        
    }

    setup$iter <- 1
    setup$maxiter <- setup$iterations
    setup$timesteps <- setup$timesteps
    setup$simulation_time <- 0

    if (is.null(setup[["store_result"]])) {
        setup$store_result <- TRUE
    }
    
    if (setup$store_result) {
        if (is.null(setup[["out_save"]])) {
            setup$out_save <- seq(1, setup$iterations)
        }
    }
    
    return(setup)
}

## This function, called only by master, stores on disk the last
## calculated time step if store_result is TRUE and increments the
## iteration counter
master_iteration_end <- function(setup) {
    iter <- setup$iter
    ## max digits for iterations
    dgts <-  as.integer(ceiling(log10(setup$iterations + 1)))
    ## string format to use in sprintf 
    fmt <- paste0("%0", dgts, "d")
    
    ## Write on disk state_T and state_C after every iteration
    ## comprised in setup$out_save
    if (setup$store_result) {
        if (iter %in% setup$out_save) {
            nameout <- paste0(fileout, "/iter_", sprintf(fmt=fmt, iter), ".rds")
            info <- list(
                tr_req_dt = as.integer(setup$req_dt)
                ## tr_allow_dt = setup$allowed_dt,
                ## tr_inniter = as.integer(setup$inniter)
            )
            saveRDS(list(
                T = setup$state_T, C = setup$state_C,
                simtime = as.integer(setup$simtime),
                tr_info = info
            ), file = nameout)
            msgm("results stored in <", nameout, ">")
        }
    }
    msgm("done iteration", iter, "/", setup$maxiter)
    setup$iter <- setup$iter + 1
    return(setup)
}

## function for the workers to compute chemistry through PHREEQC
slave_chemistry <- function(setup, data) {
    base <- setup$base
    first <- setup$first
    prop <- setup$prop
    immobile <- setup$immobile
    kin <- setup$kin
    ann <- setup$ann
    
    iter <- setup$iter
    timesteps <- setup$timesteps
    dt <- timesteps[iter]
    
    state_T <- data ## not the global field, but the work-package
    
    ## treat special H+/pH, e-/pe cases
    state_T <- RedModRphree::Act2pH(state_T)
    
    ## reduction of the problem
    if (setup$reduce) {
        reduced <- ReduceStateOmit(state_T, omit = setup$ann)
    } else {
        reduced <- state_T
    }
    
    ## form the PHREEQC input script for the current work package
    inplist <- SplitMultiKin(
        data = reduced, procs = 1, base = base, first = first,
        ann = ann, prop = prop, minerals = immobile, kin = kin, dt = dt
    )
    
    ## if (local_rank==1 & iter==1)
    ##         RPhreeWriteInp("FirstInp", inplist)
    
    tmpC <- RunPQC(inplist, procs = 1, second = TRUE)
    
    ## recompose after the reduction
    if (setup$reduce) {
        state_C <- RecomposeState(tmpC, reduced)
    } else {
        state_C <- tmpC
    }
    
    ## the next line is needed since we don't need all columns of
    ## PHREEQC output
    return(state_C[, prop])
}

## This function, called by master
master_chemistry <- function(setup, data) {
    state_T <- setup$state_T
    
    msgm(" chemistry iteration", setup$iter)
    
    ## treat special H+/pH, e-/pe cases
    state_T <- RedModRphree::Act2pH(state_T)
    
    ## reduction of the problem
    if (setup$reduce) {
        reduced <- ReduceStateOmit(state_T, omit = setup$ann)
    } else {
        reduced <- state_T
    }
    
    ## inject data from workers
    res_C <- data
    
    rownames(res_C) <- NULL
    
    ## print(res_C)
    
    if (nrow(res_C) > nrow(reduced)) {
        res_C <- res_C[seq(2, nrow(res_C), by = 2), ]
    }
    
    ## recompose after the reduction
    if (setup$reduce) {
        state_C <- RecomposeState(res_C, reduced)
    } else {
        state_C <- res_C
    }
    
    setup$state_C <- state_C
    setup$reduced <- reduced
    
    return(setup)
}


## Adapted version for "reduction"
ReduceStateOmit <- function(data, omit = NULL, sign = 6) {
    require(mgcv)
    
    rem <- colnames(data)
    if (is.list(omit)) {
        indomi <- match(names(omit), colnames(data))
        datao <- data[, -indomi]
    } else {
        datao <- data
    }
    
    datao <- signif(datao, sign)
    red <- mgcv::uniquecombs(datao)
    inds <- attr(red, "index")
    now <- ncol(red)
    
    
    ## reattach the omitted column(s)
    ## FIXME: control if more than one ann is present
    if (is.list(omit)) {
        red <- cbind(red, rep(data[1, indomi], nrow(red)))

        colnames(red)[now + 1] <- names(omit)

        ret <- red[, colnames(data)]
    } else {
        ret <- red
    }
    rownames(ret) <- NULL
    attr(ret, "index") <- inds
    return(ret)
}



## Attach the name of the calling function to the message displayed on
## R's stdout
msgm <- function(...) {
    if (local_rank == 0) {
        fname <- as.list(sys.call(-1))[[1]]
        prefix <- paste0("R: ", fname, " ::")
        cat(paste(prefix, ..., "\n"))
    }
    invisible()
}


## Function called by master R process to store on disk all relevant
## parameters for the simulation
StoreSetup <- function(setup) {

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
      #signif    = dht_final_signif,
      #proptype  = dht_final_proptype
    )
  } else {
    to_store$DHT <- FALSE
  }

  saveRDS(to_store, file = paste0(fileout, "/setup.rds"))
  msgm("initialization stored in ", paste0(fileout, "/setup.rds"))
}

GetWorkPackageSizesVector <- function(n_packages, package_size, len) {
    ids <- rep(1:n_packages, times=package_size, each = 1)[1:len]
    return(as.integer(table(ids)))
}
)";

static inline std::string init_r_library = R"(input <- readRDS("/home/max/Storage/poet/build/apps/out.rds")

grid_def <- matrix(c(2, 3), nrow = 2, ncol = 5)

library(data.table)

pqc_to_grid <- function(pqc_in, grid) {
    # Convert the input DataFrame to a data.table
    dt <- data.table(pqc_in)

    # Flatten the matrix into a vector
    id_vector <- as.vector(t(grid))

    # Initialize an empty data.table to store the results
    result_dt <- data.table()

    # Iterate over each ID in the vector
    for (id_mat in id_vector) {
        # Find the matching row in the data.table
        matching_dt <- dt[id == id_mat]

        # Append the matching data.table row to the result data.table
        result_dt <- rbind(result_dt, matching_dt)
    }

    # Remove all columns which only contain NaN
    # result_dt <- result_dt[, colSums(is.na(result_dt)) != nrow(result_dt)]

    res_df <- as.data.frame(result_dt)

    return(res_df[, colSums(is.na(res_df)) != nrow(res_df)])
}

pqc_init <- pqc_to_grid(input, grid_def)
test <- pqc_init

modify_module_sizes <- function(mod_sizes, pqc_mat, init_grid) {
    # Find all unique IDs in init_grid
    unique_ids <- unique(as.vector(init_grid$id))

    # remove rows from pqc_mat that are not in init_grid
    pqc_mat <- as.data.frame(pqc_mat)
    pqc_mat <- pqc_mat[pqc_mat$id %in% unique_ids, ]

    # Find the column indices where all rows are NA
    na_cols <- which(colSums(is.na(pqc_mat)) == nrow(pqc_mat))

    # Build cumsum over mod_sizes
    cum_mod_sizes <- cumsum(mod_sizes)

    # Find the indices where the value of na_cols is equal to the value of cum_mod_sizes
    idx <- which(cum_mod_sizes %in% na_cols)

    # Set the value of mod_sizes to 0 at the indices found in the previous step
    mod_sizes[idx] <- 0

    return(mod_sizes)
}

# mod_sizes <- c(7, 0, 4, 2, 0)

# unique_ids <- unique(as.vector(pqc_init$id))

# tmp <- as.data.frame(input)
# pqc_test <- tmp[tmp$id %in% unique_ids, ]
# na_cols <- which(colSums(is.na(pqc_test)) == nrow(pqc_test))
# cum_mod_sizes <- cumsum(mod_sizes)

# # Get the indices of the columns of cum_mod_sizes where the value of the column is equal to the value of na_cols
# idx <- which(cum_mod_sizes %in% na_cols)
# mod_sizes[idx] <- 0

# idx <- which(na_cols == cum_mod_sizes)



# idx <- which(na_cols[1] >= cum_mod_sizes)

# mod_sizes <- modify_module_sizes(mod_sizes, pqc_init, pqc_init)

# remove column with all NA
)";

#endif // POET_H
