input <- readRDS("/home/max/Storage/poet/build/apps/out.rds")

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
        matching_dt <- dt[dt$id == id_mat]

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
    unique_ids <- unique(init_grid$id)

    # remove rows from pqc_mat that are not in init_grid
    pqc_mat <- as.data.frame(pqc_mat)
    pqc_mat <- pqc_mat[pqc_mat$id %in% unique_ids, ]

    # Find the column indices where all rows are NaN
    na_cols <- which(sapply(pqc_mat, function(x) all(is.na(x))))
    # na_cols <- which(colSums(is.nan(pqc_mat)) == nrow(pqc_mat))

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
