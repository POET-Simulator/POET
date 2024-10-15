### Copyright (C) 2018-2024 Marco De Lucia, Max Luebke (GFZ Potsdam, University of Potsdam)
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
# pqc_to_grid <- function(pqc_mat, grid) {
#     # Convert the input DataFrame to a matrix
#     pqc_mat <- as.matrix(pqc_mat)

#     # Flatten the matrix into a vector
#     id_vector <- as.integer(t(grid))

#     # Find the matching rows in the matrix
#     row_indices <- match(id_vector, pqc_mat[, "ID"])

#     # Extract the matching rows from pqc_mat to size of grid matrix
#     result_mat <- pqc_mat[row_indices, ]

#     # Convert the result matrix to a data frame
#     res_df <- as.data.frame(result_mat)

#     # Remove all columns which only contain NaN
#     res_df <- res_df[, colSums(is.na(res_df)) != nrow(res_df)]

#     # Remove row names
#     rownames(res_df) <- NULL

#     return(res_df)
# }

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
    # res_df <- res_df[, colSums(is.na(res_df)) != nrow(res_df)]

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
