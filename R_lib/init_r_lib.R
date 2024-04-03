has_foreach <- require(foreach)
has_doParallel <- require(doParallel)

seq_pqc_to_grid <- function(pqc_in, grid) {
    # Convert the input DataFrame to a matrix
    dt <- as.matrix(pqc_in)

    # Flatten the matrix into a vector
    id_vector <- as.vector(t(grid))

    # Initialize an empty matrix to store the results
    result_mat <- matrix(nrow = 0, ncol = ncol(dt))

    # Iterate over each ID in the vector
    for (id_mat in id_vector) {
        # Find the matching row in the matrix
        matching_row <- dt[dt[, "ID"] == id_mat, ]

        # Append the matching row to the result matrix
        result_mat <- rbind(result_mat, matching_row)
    }

    # Convert the result matrix to a data frame
    res_df <- as.data.frame(result_mat)

    # Remove all columns which only contain NaN
    res_df <- res_df[, colSums(is.na(res_df)) != nrow(res_df)]

    # Remove row names
    rownames(res_df) <- NULL

    return(res_df)
}

par_pqc_to_grid <- function(pqc_in, grid) {
    # Convert the input DataFrame to a matrix
    dt <- as.matrix(pqc_in)

    # Flatten the matrix into a vector
    id_vector <- as.vector(t(grid))

    # Initialize an empty matrix to store the results
    # result_mat <- matrix(nrow = 0, ncol = ncol(dt))

    # Set up parallel processing
    num_cores <- detectCores()
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)

    # Iterate over each ID in the vector in parallel
    result_mat <- foreach(id_mat = id_vector, .combine = rbind) %dopar% {
        # Find the matching row in the matrix
        matching_row <- dt[dt[, "ID"] == id_mat, ]

        # Return the matching row
        matching_row
    }

    # Stop the parallel processing
    stopCluster(cl)

    # Convert the result matrix to a data frame
    res_df <- as.data.frame(result_mat)

    # Remove all columns which only contain NaN
    res_df <- res_df[, colSums(is.na(res_df)) != nrow(res_df)]

    # Remove row names
    rownames(res_df) <- NULL

    return(res_df)
}

pqc_to_grid <- function(pqc_in, grid) {
    if (has_doParallel && has_foreach) {
        print("Using parallel grid creation")
        return(par_pqc_to_grid(pqc_in, grid))
    } else {
        print("Using sequential grid creation")
        return(seq_pqc_to_grid(pqc_in, grid))
    }
}

resolve_pqc_bound <- function(pqc_mat, transport_spec, id) {
    df <- as.data.frame(pqc_mat, check.names = FALSE)
    value <- df[df$ID == id, transport_spec]

    if (is.nan(value)) {
        value <- 0
    }

    return(value)
}

add_column_after_position <- function(df, new_col, pos, new_col_name) {
    # Split the data frame into two parts
    df_left <- df[, 1:(pos)]
    df_right <- df[, (pos + 1):ncol(df)]

    # Add the new column to the left part
    df_left[[new_col_name]] <- new_col

    # Combine the left part, new column, and right part
    df_new <- cbind(df_left, df_right)

    return(df_new)
}

add_missing_transport_species <- function(init_grid, new_names, old_size) {
    # skip the ID column
    column_index <- old_size + 1

    for (name in new_names) {
        init_grid <- add_column_after_position(init_grid, rep(0, nrow(init_grid)), column_index, name)
        column_index <- column_index + 1
    }

    return(init_grid)
}