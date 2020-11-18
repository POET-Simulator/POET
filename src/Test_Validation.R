## Example script for assessing simulation results with "RFun_Eval.R"
## Time-stamp: "Last modified 2020-01-27 17:23:25 delucia"

## Read all the needed functions
source("RFun_Eval.R")


## Read all .rds in a given directory
trds <- ReadRTSims("TestDHT_Log")
## Read all .dht in a given directory
tdht <- ReadAllDHT("TestDHT_Log")

## control the content of a .dht
rc <- trds[[1]]$C
rt <- trds[[1]]$T
dht <- tdht[[1]][, 10:18]
## extract the unique rows of the first iteration (3 rows expected, as per dht)
rtu <- mgcv::uniquecombs(rt)
rcu <- mgcv::uniquecombs(rc)
all.equal(attr(rtu,"index"), attr(rcu,"index"))
colnames(dht) <- colnames(rcu)


## compare _Log _NoLog and reference
## First read all the simulations
rlog <- ReadRTSims("TestDHT_Log")   ## all default values, np 4, with --dht-log
nlog <- ReadRTSims("TestDHT_NoLog") ## all default values, np 4, without --dht-log
ref  <- ReadRTSims("TestNoDHT")     ## all default values, np 4

## Compute absolute RMSE for each variable
rmse_lnl <- ComputeErrors(rlog, nlog, FUN=RMSE)
## Compute relative RMSE for each variable using sd as norm
rel_sd_rmse_lnl <- ComputeErrors(rlog, nlog, FUN=RsdRMSE)
## Compute relative RMSE for each variable using MAX as norm
rel_max_rmse_lnl <- ComputeErrors(rlog, nlog, FUN=RmaxRMSE)

## Visualize the 2 computed relative errors
## start by defining a nice color palette using RColorBrewer
mycol <- RColorBrewer::brewer.pal(8, "Dark2")

## uncomment the next line to save this image to a pdf file
## cairo_pdf("RelErr_Dol_dhtlog_vs_nolog.pdf", height=9, width=12)
par(mfrow=c(1,2))
ErrorProgress(rel_sd_rmse_lnl,  colors=mycol, ignore=c("O2g"), log="y", ylim=c(1E-08, 1E-3), las=1, main="Rel.err dht-log vs dht, norm=SD", metric="RMSE_SD")
## add a grid
abline(h=10^-seq(8,3), col="grey", lwd=0.5, lty="dashed")
ErrorProgress(rel_max_rmse_lnl, colors=mycol, ignore=c("O2g"), log="y", ylim=c(1E-08, 1E-3), las=1, main="Rel.err dht-log vs dht, norm=MAX", metric="RMSE_MAX")
abline(h=10^-seq(8,3), col="grey", lwd=0.5, lty="dashed")

## uncomment the next line when saving to pdf
## dev.off()


## Visualize Scatter Plot of each variable
PlotScatter(rlog[[20]], nlog[[20]], labs=c("DHT LOG", "DHT"))

PlotScatter(rlog[[4]], nlog[[4]], labs=c("DHT LOG", "DHT"))

## Same as before but between dht-log and ref
rmse_logref <- ComputeErrors(rlog, ref, FUN=RMSE)
rel_sd_rmse_logref <- ComputeErrors(rlog, ref, FUN=RsdRMSE)
rel_max_rmse_logref <- ComputeErrors(rlog, ref, FUN=RmaxRMSE)

cairo_pdf("RelErr_Dol_ref_vs_dhtlog.pdf", height=9, width=12)
par(mfrow=c(1,2))
ErrorProgress(rel_sd_rmse_logref,  colors=mycol, ignore=c("O2g"), log="y", ylim=c(1E-08, 1E-3), las=1, main="Rel.err dht-log vs Ref, norm=SD", metric="RMSE_SD")
abline(h=10^-seq(8,3), col="grey", lwd=0.5, lty="dashed")
ErrorProgress(rel_max_rmse_logref,  colors=mycol, ignore=c("O2g"), log="y", ylim=c(1E-08, 1E-3), las=1, main="Rel.err dht-log vs Ref, norm=MAX", metric="RMSE_MAX")
abline(h=10^-seq(8,3), col="grey", lwd=0.5, lty="dashed")

dev.off()

x11(); PlotScatter(rlog[[20]]$C, ref[[20]]$C, labs=c("DHT log", "ref"))

rmse_nlogref <- ComputeErrors(nlog, ref, FUN=RMSE)
rel_sd_rmse_nlogref <- ComputeErrors(nlog, ref, FUN=RsdRMSE)
ErrorProgress(rel_sd_rmse_nlogref, log="y")
PlotScatter(nlog[[20]]$C, ref[[20]]$C, labs=c("DHT no log", "ref"))

## Same as before but between dht-nolog and ref
rmse_nlogref <- ComputeErrors(nlog, ref, FUN=RMSE)
rel_sd_rmse_nlogref <- ComputeErrors(nlog, ref, FUN=RsdRMSE)
rel_max_rmse_nlogref <- ComputeErrors(nlog, ref, FUN=RmaxRMSE)

cairo_pdf("RelErr_Dol_ref_vs_dht.pdf", height=9, width=12)
par(mfrow=c(1,2))
ErrorProgress(rel_sd_rmse_nlogref,  colors=mycol, ignore=c("O2g"), log="y", ylim=c(1E-08, 1E-3), las=1, main="Rel.err dht-nolog vs Ref, norm=SD", metric="RMSE_SD")
abline(h=10^-seq(8,3), col="grey", lwd=0.5, lty="dashed")
ErrorProgress(rel_max_rmse_nlogref,  colors=mycol, ignore=c("O2g"), log="y", ylim=c(1E-08, 1E-3), las=1, main="Rel.err dht-nolog vs Ref, norm=MAX", metric="RMSE_MAX")
abline(h=10^-seq(8,3), col="grey", lwd=0.5, lty="dashed")

dev.off()
