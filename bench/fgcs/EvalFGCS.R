## Time-stamp: "Last modified 2024-12-11 23:21:25 delucia"

library(PoetUtils)
library(viridis)


res <- ReadPOETSims("./res_fgcs2_96/")

pp <- PlotField(res$iter_200$C$Barite, rows = 200, cols = 200, contour = FALSE,
                nlevels=12, palette=terrain.colors)

cairo_pdf("fgcs_Celestite_init.pdf", family="serif")
par(mar=c(0,0,0,0))
pp <- PlotField((res$iter_000$Celestite), rows = 200, cols = 200,
                contour = FALSE, breaks=c(-0.5,0.5,1.5),
                palette = grey.colors, plot.axes = FALSE, scale = FALSE,
                main="Initial Celestite crystals")
dev.off()


cairo_pdf("fgcs_Ba_init.pdf", family="serif")
par(mar=c(0,0,0,0))
pp <- PlotField(log10(res$iter_001$C$Cl), rows = 200, cols = 200,
                contour = FALSE, 
                palette = terrain.colors, plot.axes = FALSE, scale = FALSE,
                main="log10(Ba)")
dev.off()



pp <- PlotField(log10(res$iter_002$C$Ba), rows = 200, cols = 200,
                contour = FALSE, palette = viridis, rev.palette = FALSE,
                main = "log10(Ba) after 5 iterations")

pp <- PlotField(log10(res$iter_200$C$`S(6)`), rows = 200, cols = 200, contour = FALSE)


str(res$iter_00)

res$iter_178$C$Barite

pp <- res$iter_043$C$Barite

breaks <- pretty(pp, n = 5)

br <- c(0, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1)

pp <- PlotField(res$iter_200$C$Barite, rows = 200, cols = 200, contour = FALSE,
                breaks = br, palette=terrain.colors)



cairo_pdf("fgcs_Barite_200.pdf", family="serif")
pp <- PlotField(log10(res$iter_200$C$Barite), rows = 200, cols = 200,
                contour = FALSE, palette = terrain.colors, plot.axes = FALSE,
                rev.palette = FALSE, main = "log10(Barite) after 200 iter")
dev.off()

ref <- ReadPOETSims("./res_fgcs_2_ref")

rei <- ReadPOETSims("./res_fgcs_2_interp1/")


timref <- ReadRObj("./res_fgcs_2_ref/timings.qs")
timint <- ReadRObj("./res_fgcs_2_interp1/timings.qs")

timref

timint

wch <- c("H","O", "Ba", "Sr","Cl", "S(6)")

rf <- data.matrix(ref$iter_001$C[, wch])
r1 <- data.matrix(rei$iter_001$C[, wch])

r1[is.nan(r1)] <- NA
rf[is.nan(rf)] <- NA

cairo_pdf("fgcs_interp_1.pdf", family="serif", width = 10, height = 7)
PlotScatter(rf, r1, which = wch, labs = c("ref", "interp"), cols = 3, log="", las = 1, pch=4)
dev.off()



head(r1)

head(rf)

rf$O
r1$O
