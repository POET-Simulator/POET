## Example script for assessing simulation results with "RFun_Eval.R"
## Time-stamp: "Last modified 2020-02-05 03:36:05 delucia"

## Read all the needed functions
source("RFun_Eval.R")



## read the simulation results of a SimDol2D variant with setup$prolong
sd <- ReadRTSims("test_prol4")

simtimes <- sapply(sd, "[","simtime")


## we still need to read the grid to perform 2D visualization
demodir <- system.file("extdata", "demo_rtwithmufits", package="Rmufits")
grid <- ReadGrid(paste0(demodir,"/d2ascii.run.GRID.SUM"))

## fix the color scale
br <- seq(0, 0.0012, length=13)


## we want animation
library(animation)

## workhorse function to be used with package "animation"
PlotAn <- function(tot, prop, grid, breaks)
{
    for (step in seq(1, length(tot))) {
        snap <- tot[[step]]$C
        time <- tot[[step]]$simtime/3600/24
        ind <- match(prop, colnames(snap))
        Plot2DCellData(snap[,ind], grid=grid, contour=FALSE, breaks=breaks, nlevels=length(breaks), scale=TRUE, main=paste0(prop," after ", time, "days"))
    }
}

### simple animations ###

## As gif
saveGIF({
    PlotAn(tot=sd, prop="Dolomite", grid=grid, breaks=br)
      }, img.name = 'Dolomite_plot', movie.name="Dol2Dbreaks.gif",
          interval = 0.5, nmax = length(sd))

## as HTML
saveHTML({
    PlotAn(tot=sd, prop="Dolomite", grid=grid, breaks=br)
      }, img.name = 'Dolomite_plot', htmlfile="dolomite.html",
          interval = 0.5, nmax = length(sd))

## For inclusion in latex
saveLatex({
    PlotAn(tot=sd, prop="Dolomite", grid=grid, breaks=br)
      }, img.name = 'Dolomite_plot',
          latex.filename = 'dolomite_prolong.tex',
          interval = 0.5, nmax = length(sd))



sd <- ReadRTSims("test_prol4")

source("RFun_Eval.R")

library(RedModRphree)
library(Rmufits)
library(RcppVTK)


#### Evaluation of discrepancies of KtzDol_6p200_1_nodht_180 / KtzDol_6p200_1_dhtlog_180
dht <- ReadRTSims("KtzDol_6p200_1_dhtlog_180")
ref <- ReadRTSims("KtzDol_6p200_1_nodht_180")

rmse <- ComputeErrors(dht, ref, FUN=RMSE)

rae  <- ComputeErrors(dht, ref, FUN=RAEmax)

mae  <- ComputeErrors(dht, ref, FUN=MAE)

rel_max_rmse <- ComputeErrors(dht, ref, FUN=RmaxRMSE)

## Visualize the 2 computed relative errors
## start by defining a nice color palette using RColorBrewer
mycol <- RColorBrewer::brewer.pal(8, "Dark2")

## uncomment the next line to save this image to a pdf file
cairo_pdf("Ref_VS_DHT.pdf", height=4, width=12)
par(mfrow=c(1,2), family="serif")

ErrorProgress(rmse,  colors=mycol, ignore=c("O2g", "pe", "Cl", "C"), log="y",
              ylim=c(1E-12, 1E-3), las=1, main="Mean RMSE, reference vs. dht", metric="RMSE")
abline(h=10^-seq(11,2), col="grey", lwd=1, lty="dashed")

ErrorProgress(mae,  colors=mycol, ignore=c("O2g", "pe", "Cl", "C"), log="y",
              ylim=c(1E-11, 1E-2), las=1, main="Max Absolute Error, reference vs. dht", metric="MAE")
abline(h=10^-seq(10,3), col="grey", lwd=1, lty="dashed")
dev.off()

cairo_pdf("Scatter_Ref_VS_DHT.pdf", height=5, width=9)
par(family="serif")
PlotScatter(ref[[200]], dht[[200]],
            which=c("Calcite", "Dolomite", "Ca", "Mg", "C", "pH"),
            cols=3, labs=c("Reference", "DHT"), pch=".", cex=2)
dev.off()

time_ref <- readRDS("KtzDol_6p200_1_nodht_180/timings.rds")
time_dht <- readRDS("KtzDol_6p200_1_dhtlog_180/timings.rds")


## export to paraview
resk <- lapply(res, function(x) return(data.matrix(x$C)))
ExportToParaview("template_vtu_ketzinlarge.vtu", "KtzDol_6p200_1/paraview", results=resk)




