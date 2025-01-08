/*
** Copyright (C) 2018-2021 Alexander Lindemann, Max Luebke (University of
** Potsdam)
**
** Copyright (C) 2018-2023 Marco De Lucia, Max Luebke (GFZ Potsdam)
**
** Copyright (C) 2023-2024 Max Luebke (University of Potsdam)
**
** POET is free software; you can redistribute it and/or modify it under the
** terms of the GNU General Public License as published by the Free Software
** Foundation; either version 2 of the License, or (at your option) any later
** version.
**
** POET is distributed in the hope that it will be useful, but WITHOUT ANY
** WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
** A PARTICULAR PURPOSE. See the GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License along with
** this program; if not, write to the Free Software Foundation, Inc., 51
** Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#pragma once

#include <cstdint>
#include <set>
#include <string>
#include <vector>

#include <Rcpp.h>

#define SRC_DIR "/mnt/scratch/signer/poet"
#define CUDA_SRC_DIR ""
#define USE_NAA "USE_NAA;ON"

// target and source IP address for NAA support
// #ifdef USE_NAA
#define SOURCE_IP "10.3.10.41"
#define TARGET_IP "10.3.10.42"
// #endif

static const char *poet_version = "naaice/v0.3-105-g13ad41d";

// using the Raw string literal to avoid escaping the quotes
static const inline std::string kin_r_library = R"(### Copyright (C) 2018-2023 Marco De Lucia, Max Luebke (GFZ Potsdam)
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

    dgts <- as.integer(ceiling(log10(setup$maxiter)))
    ## string format to use in sprintf
    fmt <- paste0("%0", dgts, "d")

    if (is.null(setup[["store_result"]])) {
        setup$store_result <- TRUE
    }

    if (setup$store_result) {
        init_field_out <- paste0(out_dir, "/iter_", sprintf(fmt = fmt, 0), ".", setup$out_ext)
        init_field <- data.frame(init_field, check.names = FALSE)
        SaveRObj(x = init_field, path = init_field_out)
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
    dgts <- as.integer(ceiling(log10(setup$maxiter + 1)))
    ## string format to use in sprintf
    fmt <- paste0("%0", dgts, "d")

    ## Write on disk state_T and state_C after every iteration
    ## comprised in setup$out_save
    if (setup$store_result) {
        if (iter %in% setup$out_save) {
            nameout <- paste0(setup$out_dir, "/iter_", sprintf(fmt = fmt, iter), ".", setup$out_ext)
            state_T <- data.frame(state_T, check.names = FALSE)
            state_C <- data.frame(state_C, check.names = FALSE)

            ai_surrogate_info <- list(
                prediction_time = if (exists("ai_prediction_time")) ai_prediction_time else NULL,
                predictions_validity = if (exists("validity_vector")) validity_vector else NULL,
                cluster_labels = if (exists("cluster_labels")) cluster_labels else NULL,
                predictions = if (exists("predictions")) predictions else NULL
            )

            SaveRObj(x = list(
                T = state_T,
                C = state_C,
                simtime = as.integer(setup$simulation_time),
                totaltime = as.integer(totaltime),
                ai_surrogate_info = ai_surrogate_info
            ), path = nameout)
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


## Handler to read R objs from binary files using either builtin
## readRDS() or qs::qread() based on file extension
ReadRObj <- function(path) {
    ## code borrowed from tools::file_ext()
    pos <- regexpr("\\.([[:alnum:]]+)$", path)
    extension <- ifelse(pos > -1L, substring(path, pos + 1L), "")

    switch(extension,
        rds = readRDS(path),
        qs  = qs::qread(path)
    )
}

## Handler to store R objs to binary files using either builtin
## saveRDS() or qs::qsave() based on file extension
SaveRObj <- function(x, path) {
    msgm("Storing to", path)
    ## code borrowed from tools::file_ext()
    pos <- regexpr("\\.([[:alnum:]]+)$", path)
    extension <- ifelse(pos > -1L, substring(path, pos + 1L), "")

    switch(extension,
        rds = saveRDS(object = x, file = path),
        qs  = qs::qsave(x = x, file = path)
    )
}
)";
static const inline std::string init_r_library = R"(### Copyright (C) 2018-2024 Marco De Lucia, Max Luebke (GFZ Potsdam, University of Potsdam)
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

##' @param pqc_mat matrix, containing IDs and PHREEQC outputs 
##' @param grid matrix, zonation referring to pqc_mat$ID 
##' @return a data.frame
pqc_to_grid <- function(pqc_mat, grid) {
    # Convert the input DataFrame to a matrix
    pqc_mat <- as.matrix(pqc_mat)

    # Flatten the matrix into a vector
    id_vector <- as.integer(t(grid))

    # Find the matching rows in the matrix
    row_indices <- match(id_vector, pqc_mat[, "ID"])

    # Extract the matching rows from pqc_mat to size of grid matrix
    result_mat <- pqc_mat[row_indices, ]

    # Convert the result matrix to a data frame
    res_df <- as.data.frame(result_mat)

    # Remove all columns which only contain NaN
    res_df <- res_df[, colSums(is.na(res_df)) != nrow(res_df)]

    # Remove row names
    rownames(res_df) <- NULL

    return(res_df)
}


##' @param pqc_mat matrix, 
##' @param transport_spec column name of species in pqc_mat
##' @param id
##' @title 
##' @return 
resolve_pqc_bound <- function(pqc_mat, transport_spec, id) {
    df <- as.data.frame(pqc_mat, check.names = FALSE)
    value <- df[df$ID == id, transport_spec]

    if (is.nan(value)) {
        value <- 0
    }

    return(value)
}

##' @title 
##' @param init_grid 
##' @param new_names 
##' @return 
add_missing_transport_species <- function(init_grid, new_names) {
    # add 'ID' to new_names front, as it is not a transport species but required
    new_names <- c("ID", new_names)
    sol_length <- length(new_names)

    new_grid <- data.frame(matrix(0, nrow = nrow(init_grid), ncol = sol_length))
    names(new_grid) <- new_names

    matching_cols <- intersect(names(init_grid), new_names)

    # Copy matching columns from init_grid to new_grid
    new_grid[, matching_cols] <- init_grid[, matching_cols]


    # Add missing columns to new_grid
    append_df <- init_grid[, !(names(init_grid) %in% new_names)]
    new_grid <- cbind(new_grid, append_df)

    return(new_grid)
}
)";
static const inline std::string ai_surrogate_r_library =
    R"(## This file contains default function implementations for the ai surrogate.
## To use pre-/postprocessing it is recommended to override these functions
## with custom implementations via the input script. The path to the R-file
## See the barite_50.R file as an example and the general README for more
## information.

preprocess <- function(df) {
  return(df)
}

postprocess <- function(df) {
  return(df)
}

set_valid_predictions <- function(temp_field, prediction, validity) {
  temp_field[validity == 1, ] <- prediction[validity == 1, ]
  return(temp_field)
}

get_invalid_values <- function(df, validity) {
  return(df[validity == 0, ])
}

set_field <- function(temp_field, columns, rows, column_name_limit,
                      byrow = FALSE) {
  temp_field <- matrix(temp_field, nrow = rows, byrow = byrow)
  temp_field <- setNames(data.frame(temp_field), columns)
  temp_field <- temp_field[column_name_limit]
  return(temp_field)
}
)";
static const inline std::string r_runtime_parameters = "mysetup";

struct RuntimeParameters {
  std::string out_dir;
  std::vector<double> timesteps;

  Rcpp::List init_params;

  bool as_rds = false;
  std::string out_ext; // MDL added to accomodate for qs::qsave/qread

  bool print_progress = false;

  static constexpr std::uint32_t WORK_PACKAGE_SIZE_DEFAULT = 32;
  std::uint32_t work_package_size;

  bool use_dht = false;
  static constexpr std::uint32_t DHT_SIZE_DEFAULT = 1.5E3;
  std::uint32_t dht_size;
  static constexpr std::uint8_t DHT_SNAPS_DEFAULT = 0;
  std::uint8_t dht_snaps;

  bool use_interp = false;
  static constexpr std::uint32_t INTERP_SIZE_DEFAULT = 100;
  std::uint32_t interp_size;
  static constexpr std::uint32_t INTERP_MIN_ENTRIES_DEFAULT = 5;
  std::uint32_t interp_min_entries;
  static constexpr std::uint32_t INTERP_BUCKET_ENTRIES_DEFAULT = 20;
  std::uint32_t interp_bucket_entries;

  /*AI surriogate configuration*/
  bool use_ai_surrogate = false; // Can be set with command line flag ---ai-surrogate
  bool disable_training = false; // Can be set in the R input script
  bool use_clustering = false; // Can be set in the R input script
  bool use_Keras_predictions = false; // Can be set in the R input script
  bool train_only_invalid = false; // Can be set in the R input script
  int batch_size = 2560; // default value determined in test on the UP Turing cluster
  int training_epochs = 20;; // Can be set in the R input script
  int training_data_size; // Can be set in the R input script
  bool use_naa = true;
  std::string save_model_path = ""; // Can be set in the R input script
  std::string cuda_src_dir = CUDA_SRC_DIR; // From CMake
};
