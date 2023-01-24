## Simple library of functions to assess and visualize the results of the coupled simulations

## Time-stamp: "Last modified 2023-01-18 16:02:58 delucia"

require(RedModRphree)
require(Rmufits)  ## essentially for PlotCartCellData 
require(Rcpp)

curdir <- dirname(sys.frame(1)$ofile) ##path.expand(".")

print(paste("RFun_Eval.R is in ", curdir))
sourceCpp(file = paste0(curdir, "/interpret_keys.cpp"))

# Wrapper around previous sourced Rcpp function
ConvertDHTKey <- function(value) {
  rcpp_key_convert(value)
}

## function which reads all simulation results in a given directory
ReadRTSims <- function(dir) {
    files_full <- list.files(dir, pattern="iter.*rds", full.names=TRUE)
    files_name <- list.files(dir, pattern="iter.*rds", full.names=FALSE)
    res <- lapply(files_full, readRDS)
    names(res) <- gsub(".rds","",files_name, fixed=TRUE)
    return(res)
}

## function which reads all successive DHT stored in a given directory
ReadAllDHT <- function(dir, new_scheme = TRUE) {
    files_full <- list.files(dir, pattern="iter.*dht", full.names=TRUE)
    files_name <- list.files(dir, pattern="iter.*dht", full.names=FALSE)
    res <- lapply(files_full, ReadDHT, new_scheme = new_scheme)
    names(res) <- gsub(".rds","",files_name, fixed=TRUE)
    return(res)
}

## function which reads one .dht file and gives a matrix
ReadDHT <- function(file, new_scheme = TRUE) {
    conn <- file(file, "rb")  ## open for reading in binary mode
    if (!isSeekable(conn))
        stop("Connection not seekable")

    ## we first reposition ourselves to the end of the file...
    tmp <- seek(conn, where=0, origin = "end")
    ## ... and then back to the origin so to store the length in bytes
    flen <- seek(conn, where=0, origin = "start")

    ## we read the first 2 integers (4 bytes each) containing dimensions in bytes
    dims <- readBin(conn, what="integer", n=2)

    ## compute dimensions of the data
    tots <- sum(dims)
    ncol <- tots/8
    nrow <- (flen - 8)/tots ## 8 here is 2*sizeof("int")
    buff <- readBin(conn, what="double", n=ncol*nrow)
    ## close connection
    close(conn)
    res <- matrix(buff, nrow=nrow, ncol=ncol, byrow=TRUE)

    if (new_scheme) {
      nkeys <- dims[1] / 8
      keys <- res[, 1:nkeys]

      conv <- apply(keys, 2, ConvertDHTKey)
      res[, 1:nkeys] <- conv
    }

    return(res)
}

## Scatter plots of each variable in the iteration 
PlotScatter <- function(sam1, sam2, which=NULL, labs=c("NO DHT", "DHT"), pch=".", cols=3, ...) {
    if ((!is.data.frame(sam1)) & ("T" %in% names(sam1)))
        sam1 <- sam1$C
    if ((!is.data.frame(sam2)) & ("T" %in% names(sam2)))
        sam2 <- sam2$C
    if (is.numeric(which))
        inds <- which
    else if (is.character(which))
        inds <- match(which, colnames(sam1))
    else if (is.null(which))
        inds <- seq_along(colnames(sam1))

    rows <- nrow(matrix(seq_along(inds), ncol=cols))
    par(mfrow=c(rows, cols))
    a <- lapply(inds, function(x) {
        plot(sam1[,x], sam2[,x], main=colnames(sam1)[x],  xlab=labs[1], ylab=labs[2], pch=pch, col="red", ...)
        abline(0,1, col="grey", cex=1.5)
    })
    invisible()
}

##### Some metrics for relative comparison

## Using range as norm
RranRMSE <- function(pred, obs)
    sqrt(mean((pred - obs)^2))/abs(max(pred) - min(pred))

## Using max val as norm
RmaxRMSE <- function(pred, obs)
    sqrt(mean((pred - obs)^2)/abs(max(pred)))

## Using sd as norm
RsdRMSE <- function(pred, obs)
    sqrt(mean((pred - obs)^2))/sd(pred)

## Using mean as norm
RmeanRMSE <- function(pred, obs)
    sqrt(mean((pred - obs)^2))/mean(pred)

## Using mean as norm
RAEmax <- function(pred, obs)
    mean(abs(pred - obs))/max(pred)

## Max absolute error
MAE <- function(pred, obs)
    max(abs(pred - obs))

## workhorse function for ComputeErrors and its use with mapply
AppliedFun <- function(a, b, .fun) mapply(.fun, as.list(a$C), as.list(b$C))

## Compute the diffs between two simulation, iter by iter,
## with a given metric (passed in form of function name to this function)
ComputeErrors <- function(sim1, sim2, FUN=RMSE) {
    if (length(sim1)!= length(sim2)) {
        cat("The simulations do not have the same length, subsetting to the shortest\n")
        a <- min(length(sim1), length(sim2))
        sim1 <- sim1[1:a]
        sim2 <- sim2[1:a]
    }
    if (!is.function(match.fun(FUN))) {
        cat("Invalid function\n")

    }

    t(mapply(AppliedFun, sim1, sim2, MoreArgs=list(.fun=FUN)))
}

## Function to display the error progress between 2 simulations
ErrorProgress <- function(mat, ignore, colors, metric, ...) {
    if (missing(colors))
        colors <- sample(rainbow(ncol(mat)))

    if (missing(metric))
        metric <- "Metric"

    ## if the optional argument "ignore" (a character vector) is
    ## passed, we remove the matching column names
    if (!missing(ignore)) {
        to_remove <- match(ignore, colnames(mat))
        mat <- mat[, -to_remove]
        colors <- colors[-to_remove]
    }
    
    yc <- mat[nrow(mat),]
    par(mar=c(5,4,2,6))
    matplot(mat, type="l", lty=1, lwd=2, col=colors, xlab="iteration", ylab=metric, ...)
    mtext(colnames(mat), side = 4, line = 2, outer = FALSE, at = yc, adj = 0.5, col = colors, las=2)
}

## Function which exports all simulations to ParaView's .vtu Requires
## package RcppVTK
ExportToParaview <- function(vtu, nameout, results) {
    require(RcppVTK)
    n <- length(results)
    vars <- colnames(results[[1]])
    ## strip eventually present ".vtu" from nameout
    nameout <- sub(".vtu", "", nameout, fixed=TRUE)
    namesteps <- paste0(nameout,  ".", sprintf("%04d",seq(1,n)), ".vtu")
    for (step in seq_along(results)) {
        file.copy(from=vtu, to=namesteps[step], overwrite = TRUE)
        cat(paste("Saving step ", step, " in file ", namesteps[step], "\n"))
        ret <- ExportMatrixToVTU (fin=vtu, fout=namesteps[step], names=colnames(results[[step]]), mat=results[[step]])
    }
    invisible(ret)
}


## Version of Rmufits::PlotCartCellData with the ability to fix the
## "breaks" for color coding of 2D simulations
Plot2DCellData <- function (data, grid, nx, ny, contour = TRUE,
                            nlevels = 12, breaks, palette = "heat.colors",
                            rev.palette = TRUE, scale = TRUE, plot.axes=TRUE, ...) {
    if (!missing(grid)) {
        xc <- unique(sort(grid$cell$XCOORD))
        yc <- unique(sort(grid$cell$YCOORD))
        nx <- length(xc)
        ny <- length(yc)
        if (!length(data) == nx * ny) 
            stop("Wrong nx, ny or grid")
    } else {
        xc <- seq(1, nx)
        yc <- seq(1, ny)
    }
    z <- matrix(round(data, 6), ncol = nx, nrow = ny, byrow = TRUE)
    pp <- t(z[rev(seq(1, nrow(z))), ])

    if (missing(breaks)) {
        breaks <- pretty(data, n = nlevels)
    }
    
    breakslen <- length(breaks)
    colors <- do.call(palette, list(n = breakslen - 1))
    if (rev.palette) 
        colors <- rev(colors)
    if (scale) {
        par(mfrow = c(1, 2))
        nf <- layout(matrix(c(1, 2), 1, 2, byrow = TRUE), widths = c(4, 
            1))
    }
    par(las = 1, mar = c(5, 5, 3, 1))
    image(xc, yc, pp, xlab = "X [m]", ylab = "Y[m]", las = 1, asp = 1,
          breaks = breaks, col = colors, axes = FALSE, ann=plot.axes,
          ...)

    if (plot.axes) {
        axis(1)
        axis(2)
    }
    if (contour) 
        contour(unique(sort(xc)), unique(sort(yc)), pp, breaks = breaks, 
            add = TRUE)
    if (scale) {
        par(las = 1, mar = c(5, 1, 5, 5))
        PlotImageScale(data, breaks = breaks, add.axis = FALSE, 
            axis.pos = 4, col = colors)
        axis(4, at = breaks)
    }
    invisible(pp)
}
