### Copyright (C) 2018-2021 Marco De Lucia (GFZ Potsdam)
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
