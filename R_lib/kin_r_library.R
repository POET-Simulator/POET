### Copyright (C) 2018-2022 Marco De Lucia, Max Luebke (GFZ Potsdam)
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
  return(setup)
}

## This function, called only by master, stores on disk the last
## calculated time step if store_result is TRUE and increments the
## iteration counter
master_iteration_end <- function(setup) {
  iter <- setup$iter
  ## MDL Write on disk state_T and state_C after every iteration
  ## comprised in setup$out_save
  # if (store_result) {
  #    if (iter %in% setup$out_save) {
  nameout <- paste0(fileout, "/iter_", sprintf("%03d", iter), ".rds")
  info <- list(
    tr_req_dt = as.integer(setup$req_dt)
#    tr_allow_dt = setup$allowed_dt,
#    tr_inniter = as.integer(setup$inniter)
  )
  saveRDS(list(
    T = setup$state_T, C = setup$state_C,
    simtime = as.integer(setup$simtime),
    tr_info = info
  ), file = nameout)
  msgm("results stored in <", nameout, ">")
  #    }
  # }
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

  ### inject data from workers
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
  to_store <- vector(mode = "list", length = 3)
  # names(to_store) <- c("Sim", "Flow", "Transport", "Chemistry", "DHT")
  names(to_store) <- c("Sim", "Transport", "DHT")

  ## read the setup R file, which is sourced in kin.cpp
  tmpbuff <- file(filesim, "r")
  setupfile <- readLines(tmpbuff)
  close.connection(tmpbuff)

  to_store$Sim <- setupfile

  # to_store$Flow <- list(
  #     snapshots  = setup$snapshots,
  #     gridfile   = setup$gridfile,
  #     phase      = setup$phase,
  #     density    = setup$density,
  #     dt_differ  = setup$dt_differ,
  #     prolong    = setup$prolong,
  #     maxiter    = setup$maxiter,
  #     saved_iter = setup$iter_output,
  #     out_save   = setup$out_save )

  to_store$Transport <- setup$diffusion
  # to_store$Chemistry <- list(
  #    nprocs   = n_procs,
  #    wp_size  = work_package_size,
  #    base     = setup$base,
  #    first    = setup$first,
  #    init     = setup$initsim,
  #    db       = db,
  #    kin      = setup$kin,
  #    ann      = setup$ann)

  if (dht_enabled) {
    to_store$DHT <- list(
      enabled   = dht_enabled,
      log       = dht_log,
      signif    = dht_final_signif,
      proptype  = dht_final_proptype
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
