## Time-stamp: "Last modified 2023-08-18 21:46:35 delucia"

library(ReacTran)
library(deSolve)
library(rootSolve)
options(width=114)


Centre <- function(N) {
    N2 <- ceiling(N/2)
    N2*(N+1)-N
}

# The model equations
HydraulicCharge <- function(t, y, parms) {
    CONC  <- matrix(nrow = parms["Np"], ncol = parms["Np"], data = y)
    ## cat("CONC[N2,N2]=",CONC[N2,N2],"\n")
    CONC[ parms["N2p"], parms["N2p"] ] <- parms["injp"]
    dCONC <- tran.2D(CONC, dx = parms["dxp"], dy = parms["dxp"],
                     D.x=parms["Kp"],
                     D.y=parms["Kp"],
                     C.x.up   = rep(parms["boundp"], parms["Np"]),
                     C.x.down = rep(parms["boundp"], parms["Np"]),
                     C.y.up   = rep(parms["boundp"], parms["Np"]),
                     C.y.down = rep(parms["boundp"], parms["Np"]))$dC
    dCONC[ Centre(parms["Np"]) ] <- 0
    return(list(dCONC))
}


DarcyVec <- function(v, bound, delta) {
    aug <- c(bound, v, bound)
    grad <- diff(aug)/delta 
}


Darcy <- function(h, dx, boundary, K) {
    nx <- ncol(h) 
    ny <- nrow(h)
    ## nx+1 x-components of \vec{U}
    Ux <- -K*apply(h, 1, DarcyVec, bound=boundary, delta=dx)
    ## ny+1 y-components of \vec{U}
    Uy <- -K*apply(h, 2, DarcyVec, bound=boundary, delta=dx)
    list(Ux=t(Ux), Uy=Uy)
}

## To reassemble Darcy velocities from cell interfaces to cell centres
RollSums <- function(vec) {
    vec[-length(vec)] + vec[-1]
}


##' @title Compute Darcy velocities assuming steady flow
##' @param N number of cells per side in 2D square grid
##' @param L domain size (m)
##' @param inj_h fixed hydraulic charge at domain center
##' @param bound_h fixed hydraulic charge at all boundaries
##' @param K permeability
##' @return 
##' @author 
GenerateVelocities <- function(N, L, inj_h=20, bound_h=1, K=1E-4) {
    require(ReacTran)
    require(deSolve)
    
    ## Construction of the 2D grid
    x.grid <- setup.grid.1D(x.up = 0, L = L, N = N)
    y.grid <- setup.grid.1D(x.up = 0, L = L, N = N)
    grid2D <- setup.grid.2D(x.grid, y.grid)
    
    dx <- L/N

    ## rough "center" of square domain
    N2    <- ceiling(N/2)

    ## initial condition: "boundary" everywhere...
    y <- matrix(nrow = N, ncol = N, data = bound_h)
    ## except in the central point!
    y[N2, N2] <- inj_h

    pars <- c(Kp=K, boundp = bound_h, Np=N, N2p=N2, injp=inj_h, dxp=dx)
    cat("\n:: calling ode.2D...")
    outc <- deSolve::ode.2D(y = y, func = HydraulicCharge,
                            times = c(0, 1E7), parms = pars,
                            dim = c(N, N), lrw = 71485484)
    
    cat(" [ DONE ]\n")
    charge <- matrix(outc[2,-1], nrow = N, ncol = N)


    cat(":: Computing charge gradient...")
    U <- Darcy(charge, dx=dx, boundary=bound_h, K=K) 
    cat(" [ DONE ]\n")

    ## Reassemble the total velocities per cell centre
    fx <- -t(apply(U$Ux, 1, RollSums))
    fy <-   -apply(U$Uy, 2, RollSums)

    x <-  y <- seq(dx/2, L - dx/2, length=N)
    xx <- rep(x, N)
    yy <- rep(y, each=N)

    norms <- sqrt(fx**2+fy**2)
    tab <- data.frame(x=xx, y=yy, norm=as.numeric(norms), ux=as.numeric(fx), uy=as.numeric(fy))

    list(Charge=charge, outc=outc, U=U, xy=tab)
}


PlotFlow <- function(tab, skip=1, scale=TRUE, arr.len=0.05,
                     arr.lwd = 1, ...) {
    if (scale) {
        tab$ux <- scale(tab$ux, center=FALSE)
        tab$uy <- scale(tab$uy, center=FALSE)
    } 

    ngrid <- sqrt(nrow(tab))
    dx    <- 2*tab$x[1]
    fieldlen <- dx*ngrid

    plot(0, 0, "n", xlim = c(0,fieldlen), ylim = c(0,fieldlen),
         asp = 1, xlab="EASTING", ylab="NORTHING", las=1, ...)
    
    arrows(x0 = tab$x, y0 = tab$y,
           x1 = tab$x - tab$uy, y1 = tab$y - tab$ux,
           length=arr.len, lwd=arr.lwd)
    
    invisible()
}


## square domain of 100x100 m discretized in 101x101 cells
aa <- GenerateVelocities(N=201, L=100, inj_h=20, bound_h=1, K=1E-2)

## plot hydraulic charge (via deSolve's method!)
image(aa$outc, ask = FALSE, main ="Hydraulic charge",
      legend = TRUE, add.contour = FALSE, subset = time == 1E7,
      xlab="", ylab="", axes=FALSE, asp=1)

x11()

PlotFlow(aa$xy, skip=1, scale=TRUE, arr.lwd = 0.5)


image(matrix(log10(aa$xy$norm), ncol=201), asp=1)
