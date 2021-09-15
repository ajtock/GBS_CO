#!/applications/R/R-3.3.2/bin/Rscript

# Plot histograms of COs per individual

# Usage:
# ./hist_COs_per_ind.R 

library(GenomicRanges)

inDir <- "/projects/ajt200/GBS_CO/All_081018/coordinates/pooled_lanes/"
plotDir <- paste0(inDir, "hist/")
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

filePath <- system(paste0("ls ", inDir, "*_pooled_GRanges.RData"),
                   intern = T)
popName <- lapply(seq_along(filePath), function(x) {
  substr(filePath[x],
         start = 61,
         stop = length(strsplit(filePath[x], "")[[1]])-21)
})

pdf(paste0(plotDir, "hist_COs_per_ind_xlim55.pdf"),
    height = 1.5*round((length(filePath)/2)), width = 5.5)
par(mfrow = c(round(length(filePath)/2), 2))
par(mar=c(3.2, 3.2, 2.1, 2.1))
for(x in 1:length(filePath)) {
  load(filePath[x])
  COs_per_ind <- NULL
  for(j in 1:length(unique(popCOsGR_pooled$lib))) {
     ind_COsGR <- popCOsGR_pooled[popCOsGR_pooled$lib == unique(popCOsGR_pooled$lib)[j]]
     COs_per_ind <- c(COs_per_ind, length(ind_COsGR))
  }
  hist(COs_per_ind,
       xlim = c(0, 55),
       breaks = 10000,
       main = popName[[x]],
       cex.main = 0.9,
       xlab = "",
       ylab = "")
  mtext(side = 1, line = 2, cex = 0.7, text = "Crossovers")
  mtext(side = 2, line = 2, cex = 0.7, text = "Individuals")
  abline(v = mean(COs_per_ind), col = "red", lty = 2)
}
dev.off()

pdf(paste0(plotDir, "hist_COs_per_ind_xlim100.pdf"),
    height = 1.5*round((length(filePath)/2)), width = 5.5)
par(mfrow = c(round(length(filePath)/2), 2))
par(mar=c(3.2, 3.2, 2.1, 2.1))
for(x in 1:length(filePath)) {
  load(filePath[x])
  COs_per_ind <- NULL
  for(j in 1:length(unique(popCOsGR_pooled$lib))) {
     ind_COsGR <- popCOsGR_pooled[popCOsGR_pooled$lib == unique(popCOsGR_pooled$lib)[j]]
     COs_per_ind <- c(COs_per_ind, length(ind_COsGR))
  }
  hist(COs_per_ind,
       xlim = c(0, 100),
       breaks = 10000,
       main = popName[[x]],
       cex.main = 0.9,
       xlab = "",
       ylab = "")
  mtext(side = 1, line = 2, cex = 0.7, text = "Crossovers")
  mtext(side = 2, line = 2, cex = 0.7, text = "Individuals")
  abline(v = mean(COs_per_ind), col = "red", lty = 2)
}
dev.off()

