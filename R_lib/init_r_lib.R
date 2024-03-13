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

# remove column with all NA
