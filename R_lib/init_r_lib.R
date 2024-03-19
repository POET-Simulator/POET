pqc_to_grid <- function(pqc_in, grid) {
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
