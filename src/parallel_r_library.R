distribute_work_packages <- function(len, package_size) 
{
    ## Check if work_package is a divisor of grid length and act
    ## accordingly
    if ((len %% package_size) == 0) {
        n_packages <- len/package_size
    } else {
        n_packages <- floor(len/package_size) + 1
    }

    ids <- rep(1:n_packages, times=package_size, each = 1)[1:len]
    return(ids)
}

compute_wp_sizes <- function(ids)
{
    as.integer(table(ids))
}

shuffle_field <- function(data, send_ord) 
{
    shf <- data[send_ord,]
    return(shf)
}

unshuffle_field <- function(data, send_ord) 
{
    data[send_ord,] <- data
    rownames(data) <- NULL
    return(data)
}

stat_wp_sizes <- function(sizes)
{
    res <- as.data.frame(table(sizes))
    if (nrow(res)>1) {
        msgm("Chosen work_package_size is not a divisor of grid length. ")
        colnames(res) <- c("Size", "N")
        print(res)
    } else {
        msgm("All work packages of length ", sizes[1])
    }
}
