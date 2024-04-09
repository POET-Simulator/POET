pqc_to_grid <- function(pqc_in, grid) {
    # Convert the input DataFrame to a matrix
    pqc_in <- as.matrix(pqc_in)

    # Flatten the matrix into a vector
    id_vector <- as.numeric(t(grid))

    # Find the matching rows in the matrix
    row_indices <- match(id_vector, pqc_in[, "ID"])

    # Extract the matching rows from pqc_in to size of grid matrix
    result_mat <- pqc_in[row_indices, ]

    # Convert the result matrix to a data frame
    res_df <- as.data.frame(result_mat)

    # Remove all columns which only contain NaN
    res_df <- res_df[, colSums(is.na(res_df)) != nrow(res_df)]

    # Remove row names
    rownames(res_df) <- NULL

    return(res_df)
}

resolve_pqc_bound <- function(pqc_mat, transport_spec, id) {
    df <- as.data.frame(pqc_mat, check.names = FALSE)
    value <- df[df$ID == id, transport_spec]

    if (is.nan(value)) {
        value <- 0
    }

    return(value)
}

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