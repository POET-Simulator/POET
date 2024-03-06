## Time-stamp: "Last modified 2022-12-16 20:26:03 delucia"

source("../../../util/data_evaluation/RFun_Eval.R")

sd <- ReadRTSims("naaice_2d")

sd <- ReadRTSims("Sim2D")


sd <- ReadRTSims("inner")

tim <- readRDS("inner/timings.rds")


simtimes <- sapply(sd, "[","simtime")

## workhorse function to be used with package "animation"
PlotAn <- function(tot, prop, grid, breaks) {
    for (step in seq(1, length(tot))) {
        snap <- tot[[step]]$C
        time <- tot[[step]]$simtime/3600/24
        ind <- match(prop, colnames(snap))
        Plot2DCellData(snap[,ind], grid=grid, contour=FALSE, breaks=breaks, nlevels=length(breaks), scale=TRUE, main=paste0(prop," after ", time, "days"))
    }
}


options(width=110)
library(viridis)

Plot2DCellData(sd$iter_050$C$Cl, nx=1/100, ny=1/100, contour = TRUE,
               nlevels = 12, palette = "heat.colors",
               rev.palette = TRUE, scale = TRUE, main="Cl")

Plot2DCellData(sd$iter_050$C$Dolomite, nx=100, ny=100, contour = FALSE,
               nlevels = 12, palette = "heat.colors",
               rev.palette = TRUE, scale = TRUE, )

cairo_pdf("naaice_inner_Dolo.pdf", width=8, height = 6, family="serif")
Plot2DCellData(sd$iter_100$C$Dolomite, nx=100, ny=100, contour = FALSE,
               nlevels = 12, palette = "viridis",
               rev.palette = TRUE, scale = TRUE, plot.axes = FALSE,
               main="2D Diffusion - Dolomite after 2E+4 s (100 iterations)")
dev.off()

cairo_pdf("naaice_inner_Mg.pdf", width=8, height = 6, family="serif")
Plot2DCellData(sd$iter_100$C$Mg, nx=100, ny=100, contour = FALSE,
               nlevels = 12, palette = "terrain.colors",
               rev.palette = TRUE, scale = TRUE, plot.axes=FALSE,
               main="2D Diffusion - Mg after 2E+4 s (100 iterations)")
dev.off()
