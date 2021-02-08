## Simple function to check file extension. It is needed to check if
## the GridFile is SUM (MUFITS format) or rds/RData
FileExt <- function (x)
{
    pos <- regexpr("\\.([[:alnum:]]+)$", x)
    ifelse(pos > -1L, substring(x, pos + 1L), "")
}

### This function is only called by the master process. It sets up
### the data structures for the coupled simulations and performs the
### first timestep or "iteration zero", possibly involving the phreeqc
### calculation of the initial state (if it is not already provided as
### matrix) and the first Advection iteration
master_init <- function(setup)
{
    msgm("Process with rank 0 reading GRID and MUFITS_sims")

    ## MDL TODO: actually grid and all the MUFITS snapshots should be
    ## read and stored only by the master process!!

    ## Setup the directory where we will store the results
    verb <- FALSE
    if (local_rank==0) {
        verb <- TRUE ## verbosity loading MUFITS results
        if (!dir.exists(fileout)) {
            dir.create(fileout)
            msgm("created directory ", fileout)
        } else {
            msgm("dir ", fileout," already exists, I will overwrite!")
        }
        if (!exists("store_result"))
            msgm("store_result doesn't exist!")
        else
            msgm("store_result is ", store_result)

    } else {
        
    }

    ## to enhance flexibility, gridfile can be now given in SUM or,
    ## already processed, in rds format. We check for extension first. 
    gridext <- FileExt(setup$gridfile)
    if (gridext=="SUM") {
        msgm("Generating grid from MUFITS SUM file...")
        mufits_grid <- ReadGrid(setup$gridfile, verbose=verb)
    } else {
        msgm("Reading grid from rds file...")
        mufits_grid <- readRDS(setup$gridfile)
    }

    ## load all the snapshots at once. The setup argument "snapshots"
    ## now can be either a *directory* where the .SUM files are stored
    ## or a rds file.
    if (dir.exists(setup$snapshots)) {
        msgm("<snapshots> points to a directory; reading from SUM files in there...")
        mufits_sims <- LoadMufitsSumRes(dir=setup$snapshots, verbose=verb)
    } else {
        msgm("<snapshots> points to a file. Reading as rds...")
        mufits_sims <- readRDS(setup$snapshots)
    }
    
    ## Since this function is evaluated by the R process called from
    ## the C++ program, we need to make available all these variables
    ## to the R parent frame!
    assign("mufits_grid", mufits_grid, pos=parent.frame())
    assign("mufits_sims", mufits_sims, pos=parent.frame())
    ## cat(local_rank, "assignement complete\n")


    ## we calculate the *coupling* iterations
    nstep <- length(mufits_sims)
    msgm("using", nstep,"flow snapshots")
    ## cat(local_rank, "nstep:", nstep, "\n")
    ## cat(local_rank, "names:", paste(names(mufits_sims), collate=","), "\n")

    ## we have these snapshots; we output the results of the coupling
    ## after each timesteps
    timesteps <- diff(sapply(mufits_sims, function(x) {return(x$info)}))
    ## cat(local_rank, "timesteps:", paste0(timesteps, collate=" "), "\n")

    maxiter <- length(timesteps)

    ## steady case not treated in this version!!!!
    steady <- FALSE

    setup$iter <- 1
    setup$maxiter <- maxiter
    setup$timesteps <- timesteps
    setup$simulation_time <- 0
    
    setup$steady <- steady
    
    if (nrow(setup$bound)==1) {
        boundmatAct <- t(RedModRphree::pH2Act(setup$bound))
        msg("formed correct matrix from setup$bound:")
        print(boundmatAct)
    } else {
        boundmatAct <- RedModRphree::pH2Act(setup$bound)
    }

    setup$boundmatAct <- boundmatAct
    return(setup)
}

## This function is called by all processes.
## Since the worker processes need state_C in worker.cpp -> worker_function()
## even the worker need to run this code.
## TODO: worker.cpp -> worker_function() --> state_C is only needed to obtain headers 
## and size of dataframe ... if this will be adjusted, this function will be 
## unnecessary for workers
init_chemistry <- function(setup) 
{
    ## setup the chemistry
    if (!is.matrix(setup$initsim)) {
        msgm("initial state defined through PHREEQC simulation, assuming homogeneous medium!")
        tmpfirstrun <- RunPQC_Warnings(setup$initsim, second=TRUE)
        
        ## if (local_rank==0){
        ##     print(tmpfirstrun)
        ## }
        remcolnames <- colnames(tmpfirstrun)
        if (nrow(tmpfirstrun) > 1) {
            ## msgm("Iter 0 selecting second row")
            firstrun <- matrix(tmpfirstrun[2,], nrow=1, byrow=TRUE)
            colnames(firstrun) <- remcolnames
        } else {
            firstrun  <- tmpfirstrun
        }
        state_C <- matrix(rep(firstrun,setup$n), byrow=TRUE, ncol=length(setup$prop))
        colnames(state_C) <- colnames(firstrun)

        ## if (local_rank==0)
        ##     saveRDS(state_C, "initstate.rds")
    } else {
        msgm("given initial state")
        state_C <- setup$initsim
    }
    ## save state_T and state_C
    setup$state_C <- state_C

    return(setup)
}

## relic code, may be useful in future
finit <- function(setup)
{
    ## saved <- setup$saved
    ## state_C <- setup$state_C
    ## timesteps <- setup$timesteps
    ## out_save <- setup$out_save
    ## to_save <- setup$to_save

    ## if (out_save) {
    ##     msgm("saving <saved>")
    ##     attr(saved,"savedtimesteps") <- timesteps[save]
    ##     attr(saved,"timesteps") <- timesteps
    ##     return(saved)
    ## } else {    
    ##     attr(saved,"timesteps") <- timesteps
    ##     return(state_C)
    ## }
}

## This function, called only by the master, computes the *inner*
## timesteps for transport given the requested simulation timestep and
## the current snapshot (through the Courant Fredrich Levy condition)
## NB: we always iterate: 1)T 2)C
master_iteration_setup <- function(setup)
{
    ## in this version, "grid" and "mufits_sims" are variables already
    ## living in the environment of the R process with local_rank 0

    ## if (local_rank == 0) {
    ##     cat("Process ", local_rank, "iteration_start\n")
    ## } else {
    ##     cat("Process ", local_rank, "SHOULD NOT call iteration_start!\n")
    ## }


    ## in C++, "iter" starts from 1
    iter <- setup$iter
    
    next_dt <- setup$timesteps[iter]

    ## setting the current snapshot
    snap <- mufits_sims[[iter]]

    msgm("Current simulation time:", round(setup$simulation_time),"[s] or ", round(setup$simulation_time/3600/24,2), "[day]")
    msgm("Time step for current iteration:", round(next_dt),"[s] or ", round(next_dt/3600/24,2), "[day]")
    setup$simulation_time <- setup$simulation_time+next_dt
    msgm("Target simulation time:", round(setup$simulation_time),"[s] or ", round((setup$simulation_time)/3600/24,2), "[day]")

    
    ## since the phase and density may change in MUFITS
    ## simulations/snapshots, we MUST tell their name at setup. Here
    ## we match the right column for both variables
    nphase   <- match(setup$phase,   colnames(snap$conn))
    ndensity <- match(setup$density, colnames(snap$cell))

    ## the "new cfl"
    flux <- snap$conn[, nphase]
    cellporvol <- mufits_grid$cell$PORO * mufits_grid$cell$CELLVOL
    
    ## in MUFITS, positive flux means from CELLID1 to CELLID2
    pos_flux <- ifelse(flux>0, snap$conn$CONNID1, snap$conn$CONNID2)
    poro <- mufits_grid$cell$PORO[pos_flux]
    vol <- mufits_grid$cell$CELLVOL[snap$conn$CONNID1]

    ## extract the right density for the right fluxes
    density <- snap$cell[pos_flux, ndensity]
    
    ## transform flux from Mufits ton/day to m3/s. This is a VOLUMETRIC FLUX
    fluxv <- flux/3600/24*1000/density
    
    ## now we want the darcy velocity
    excl0 <- which(abs(fluxv)<.Machine$double.eps)
    ## msgm("excl0"); print(excl0)
    ## msgm("length(excl0)"); print(length(excl0))
    ## The CFL condition is expressed in terms of total flux and total
    ## pore volume
    cfl <- as.numeric(min(abs(vol*poro/fluxv)[-excl0]))

    if (!is.finite(cfl))
        stop(msgm("CFL is ", cfl,"; quitting"))
             
    allowed_dt <- setup$Cf*cfl
    requested_dt <- next_dt ## target_time - setup$simulation_time
    msgm("CFL allows dt of <", round(cfl, 2)," [s]; multiplicator is ", setup$Cf,
         "; requested_dt is: ", round(requested_dt, 2))
        
    if (requested_dt > allowed_dt) {
        inniter <- floor(requested_dt/allowed_dt)
        inner_timesteps <- c(rep(allowed_dt, inniter), requested_dt - allowed_dt * inniter)
        msgm("Performing ", inniter, " inner iterations")
    } else {
        inner_timesteps <- requested_dt
        inniter <- 1
    }

    setup$inner_timesteps <- inner_timesteps
    setup$requested_dt <- requested_dt
    setup$allowed_dt <- allowed_dt
    setup$inniter <- inniter
    setup$iter <- iter
    setup$fluxv <- fluxv
    setup$no_transport_cells <- excl0
    setup$cellporvol <- cellporvol
    
    return(setup)
}


## This function, called only by master, stores on disk the last
## calculated time step if store_result is TRUE and increments the
## iteration counter
master_iteration_end <- function(setup) {
    iter    <- setup$iter
    ## MDL Write on disk state_T and state_C after every iteration
    if (store_result) {
        nameout <- paste0(fileout, '/iter_', sprintf('%03d', iter), '.rds')
        saveRDS(list(T=setup$state_T, C=setup$state_C), file=nameout)
        msgm("results stored in <", nameout, ">")
    }
    msgm("done iteration", iter, "/", setup$maxiter)
    setup$iter <- setup$iter + 1
    return(setup)
}

## master: compute advection
master_advection <- function(setup) {
    msgm("requested dt=", setup$requested_dt)
    concmat <- RedModRphree::pH2Act(setup$state_C)
    inner_timesteps <- setup$inner_timesteps
    Cf <- setup$Cf 
    iter <- setup$iter
    excl <- setup$no_transport_cells
    ## msgm("excl"); print(excl)
    immobile <- setup$immobile
    boundmat <- setup$boundmatAct
    snap <- mufits_sims[[setup$iter]]

    ## conc is a matrix with colnames 
    if (is.matrix(concmat)) {
        if (length(immobile)>0)
            val <- concmat[,-immobile]
        else
            val <- concmat
        
        ## sort the columns of the matrix matching the names of boundary 
        sptr <- colnames(val)
        inflow <- boundmat[, colnames(boundmat)!="cbound", drop=FALSE]
        spin <- colnames(inflow) 
        ind <- match(spin,sptr)
        if (TRUE %in% is.na(ind)) {
            msgm("Missing boundary conditions for some species")
            val <- cbind(val,matrix(0,ncol=sum(is.na(ind)),nrow=nrow(val) ))
            colnames(val) <- c(sptr,spin[is.na(ind)])
        }
        
        sptr <- colnames(val)
        ind <- match(sptr,spin)
        ## msgm("ind"); print(ind)
        
        cnew <- val

        msgm("Computing transport of ", ncol(val), " species")
        ## if (local_rank==0) {
        ##     saveRDS(list(conc=newconc, times=times, activeCells=grid$cell$CELLID[-exclude_cell], fluxv=fluxv),
        ##             file=paste0("TranspUpwind_", local_rank, ".RData"))
        ## }
        ## cnew[ boundmat[, 1] ] <- boundmat[, 2]

        ## cells where we won't update the concentrations: union of
        ## boundary and no-flow cells
        if (length(excl) == 0)
            exclude_cell <- sort(c(boundmat[, 1]))
        else
            exclude_cell <- sort(c(boundmat[, 1], excl))

        ## msgm("mufits_grid$cell$CELLID[-excludecell]:")
        ## print(mufits_grid$cell$CELLID[-exclude_cell])

        for (i in seq(1, ncol(val))) {
            ## form a 2-column matrix with cell id and boundary
            ## concentration for those elements
            newboundmat <- boundmat[,c(1, ind[i]+1)]
            ## vector with the old concentrations
            concv <- val[,i]
            ## apply boundary conditions to the concentration vector
            ## (they should stay constant but better safe than sorry)
            concv[newboundmat[, 1]] <- newboundmat[, 2]

            ## call the function
            cnew[,i] <- CppTransportUpwindIrr(concv = concv,
                                              times = inner_timesteps, 
                                              activeCells = mufits_grid$cell$CELLID[-exclude_cell],
                                              fluxv = setup$fluxv, 
                                              listconn = mufits_grid$listconn,
                                              porVol = setup$cellporvol)
        }
        
        colnames(cnew) <- colnames(val)
        
        if ( length(immobile) > 0)  { 
            res <- cbind(cnew, concmat[,immobile])
        } else {
            res <- cnew
        }

        ## check for negative values. This SHOULD NOT OCCUR and may be
        ## the result of numerical dispersion (or bug in my transport
        ## code!)
        if (any(res <0 )) { 
            a <- length(res[res < 0 ])
            res[res < 0 ] <- 0
            msgm("-> ", a, "concentrations were negative")
        }
        ## msgm("state_T after iteration", setup$iter, ":")
        ## print(head(res))
    } else {
        msgm("state_C at iteration", setup$iter, " is not a Matrix, doing nothing!")
    }

    setup$state_T <- res
    
    return(setup)
}


## function for the workers to compute chemistry through PHREEQC
slave_chemistry <- function(setup, data)
{
    base     <- setup$base
    first    <- setup$first
    prop     <- setup$prop
    immobile <- setup$immobile
    kin      <- setup$kin
    ann      <- setup$ann

    if (local_rank == 0) {
        iter     <- setup$iter
        timesteps <- setup$timesteps
        dt       <- timesteps[iter]
    }

    state_T <- data ## not the global field, but the work-package

    ## treat special H+/pH, e-/pe cases
    state_T <- RedModRphree::Act2pH(state_T)
    
    ## reduction of the problem
    if(setup$reduce)
        reduced <- ReduceStateOmit(state_T, omit=setup$ann)
    else
        reduced <- state_T

    ## if (local_rank==1) {
    ##     msg("worker", local_rank,"; iter=", iter, "dt=", dt)
    ##     msg("reduce is", setup$reduce)
    ##     msg("data:")
    ##     print(reduced)
    ##     msg("base:")
    ##     print(base)
    ##     msg("first:")
    ##     print(first)
    ##     msg("ann:")
    ##     print(ann)
    ##     msg("prop:")
    ##     print(prop)
    ##     msg("immobile:")
    ##     print(immobile)
    ##     msg("kin:")
    ##     print(kin)
    ## }

    ## form the PHREEQC input script for the current work package
    inplist <- splitMultiKin(data=reduced, procs=1, base=base, first=first,
                             ann=ann, prop=prop, minerals=immobile, kin=kin, dt=dt)

    ## if (local_rank==1 & iter==1) 
    ##         RPhreeWriteInp("FirstInp", inplist)

    tmpC <- RunPQC_Warnings(inplist, procs=1, second=TRUE)

    ## recompose after the reduction
    if (setup$reduce)
        state_C <- RecomposeState(tmpC, reduced)
    else {
        state_C <- tmpC
    }

    ## the next line is needed since we don't need all columns of
    ## PHREEQC output
    return(state_C[, prop])
} 



## This function, called by master
master_chemistry <- function(setup, data)
{
    state_T <- setup$state_T 

    msgm(" chemistry iteration", setup$iter)

    ## treat special H+/pH, e-/pe cases
    state_T <- RedModRphree::Act2pH(state_T)
    
    ## reduction of the problem
    if(setup$reduce)
        reduced <- ReduceStateOmit(state_T, omit=setup$ann)
    else
        reduced <- state_T
    
    ### inject data from workers
    res_C <- data

    rownames(res_C) <- NULL
    
    ## print(res_C)
    
    if (nrow(res_C) > nrow(reduced)) {
        res_C <- res_C[seq(2,nrow(res_C), by=2),]
    }
    
    ## recompose after the reduction
    if (setup$reduce)
        state_C <- RecomposeState(res_C, reduced)
    else {
        state_C <- res_C
    }
    
    setup$state_C <- state_C
    setup$reduced <- reduced

    return(setup)
} 
    

## Adapted version for "reduction"
ReduceStateOmit <- function (data, omit=NULL, sign=6) 
{
    require(mgcv)

    rem <- colnames(data)
    if (is.list(omit)) {
        indomi <- match(names(omit), colnames(data))
        datao <- data[, -indomi]
    } else datao <- data

    datao <- signif(datao, sign)
    red <- mgcv::uniquecombs(datao)
    inds <- attr(red, "index")
    now <- ncol(red)

    
    ## reattach the omitted column(s)
    ## FIXME: control if more than one ann is present
    if (is.list(omit)) {
        red <- cbind(red, rep( data[ 1, indomi], nrow(red)))
    
        colnames(red)[now+1] <- names(omit)

        ret <- red[, colnames(data)]
    } else {
        ret <- red
    }
    rownames(ret) <- NULL
    attr(ret, "index") <- inds
    return(ret)
}



## Attach the name of the calling function to the message displayed on
## R's stdout
msgm <- function (...)
{
    
    if (local_rank==0) {
        fname <- as.list(sys.call(-1))[[1]]
        prefix <- paste0("R: ", fname, " ::")
        cat(paste(prefix, ..., "\n"))
    }
    invisible()
}


RunPQC_Warnings <- function (input, procs = 1, second = TRUE) 
{

    .runPQC <- function(input, onlysecond) {
        require(phreeqc)
        phreeqc::phrSetErrorStringsOn(TRUE)
        error_count <- phreeqc::phrRunString(input)
        warnings <- phreeqc::phrGetWarningStrings()
        checked_warn <- CheckWarnings(warnings)

        if (is.null(error_count)) {
            out <- phreeqc::phrGetSelectedOutput()[[1]]
            nso <- colnames(out)
            nso <- sub("..mol.kgw.", "", nso, fixed = TRUE)
            nso <- sub(".mol.kgw.", "", nso, fixed = TRUE)
            nso <- sub("^k_", "", nso)
            nso <- sub(".g.", "(g)", nso, fixed = TRUE)
            if (second)
                out <- out[seq(2, nrow(out), by = 2), ]
            colnames(out) <- nso
            return(data.matrix(out))
        } else {
            msgm("Error in simmulation!")
            cat(phrGetErrorStrings(), sep="\n")
        }
    }
    
    if (!is.list(input)) {
        if (is.character(input)) {
            res <- .runPQC(input, onlysecond = second)
        } else {
            stopmsg("something wrong with the input, dying!")
        }
    }
    else {
        res <- foreach(i = seq_along(input), .combine = rbind) %dopar%
            .runPQC(input[[i]], onlysecond = second)
    }
    return(res)
}


CheckWarnings <- function(warnings)
{
    ## get rid of "Negative moles in solution" warnings
    recover <- grep("Negative moles", warnings, fixed=TRUE)
    warnings <- warnings[-recover]

    lw <- length(warnings)
    if (lw>0) {
        ## nfile <- paste0("warnings_", local_rank, ".txt")
        nfile <- "Warnings.txt"
        cat(warnings, file=nfile, sep="\n", append=TRUE)
    }
    return(lw)
}
