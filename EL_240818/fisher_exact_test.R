#!/applications/R/R-3.4.0/bin/Rscript

# ./fisher_exact_test.R EL_wt/wt_CO_density_genome_1Mb.txt EL_taf4b/taf4b_CO_density_genome_1Mb.txt 1Mb

#inDir <- "/projects/ajt200/GBS_CO/EL_240818/"
#wtFile <- "EL_wt/wt_CO_density_genome_1Mb.txt"
#mutantFile <- "EL_taf4b/taf4b_CO_density_genome_1Mb.txt"
#winName <- "1Mb"

args <- commandArgs(trailingOnly = TRUE)
inDir <- "./"
wtFile <- args[1]
mutantFile <- args[2]
winName <- as.character(args[3])

wt <- read.table(paste0(inDir, wtFile))
mutant <- read.table(paste0(inDir, mutantFile))

contingencyTable <- function(wtInWin, wtOutWin, mutantInWin, mutantOutWin) {
  conTab <- matrix(c(wtInWin, wtOutWin, mutantInWin, mutantOutWin), ncol = 2)
  rownames(conTab) <- c("InWin", "OutWin")
  colnames(conTab) <- c("wt", "mutant")
  conTab
}

fisherPvals <- lapply(1:dim(wt)[1], function(x) {
  fisher.test(contingencyTable(wtInWin = wt$COsChrWin[x],
                               wtOutWin = sum(wt$COsChrWin)-wt$COsChrWin[x],
                               mutantInWin = mutant$COsChrWin[x],
                               mutantOutWin = sum(mutant$COsChrWin)-mutant$COsChrWin[x]),
              alternative = "two.sided")$p.value
})

# Genomic definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# Pericentromeric regions are as defined in Supplemental Table S26
# of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

# Make cumulative genomes
sumchr <- cumsum(c(0, chrLens))
print(sumchr)
#[1] 0  30427671  50125960  73585790  92170846 119146348
sumchr_tot <- sumchr[length(sumchr)]
print(sumchr_tot)

centromeres <- sapply(seq_along(centromeres), function(x) {
  centromeres[x] + sumchr[x]
})
print(centromeres)
pericenStart <- sapply(seq_along(pericenStart), function(x) {
  pericenStart[x] + sumchr[x]
})
print(pericenStart)
pericenEnd <- sapply(seq_along(pericenEnd), function(x) {
  pericenEnd[x] + sumchr[x]
})
print(pericenEnd)

# Function to plot genome-scale profile of wt vs mutant (one Y-axis) 
wtVmutantGenomePlot <- function(xplot,
                                wt,
                                mutant,
                                Ylabel,
                                legendLoc,
                                legendLabs,
                                wtCol,
                                mutantCol) {
  plot(xplot, wt, type = "l", lwd = 1.5, col = wtCol,
       ylim = c(min(c(wt, mutant)),
                max(c(wt, mutant))),
       xlab = "",
       ylab = "", yaxt = "n",
       main = bquote(italic("r"[s]) ~ " = " ~ .(round(cor(wt, mutant, method = "spearman", use = "pairwise.complete.obs"), digits = 2))))
  lines(xplot, mutant, type = "l", lwd = 1.5, col = mutantCol)
  axis(side = 2, at = pretty(c(wt, mutant)), lwd.tick = 1.5)
  mtext(side = 2, line = 2.25, cex = 1, text = Ylabel, col = "black")
  abline(v = sumchr, lty = 1, lwd = 0.75)
  abline(v = centromeres, lty = 5, lwd = 0.75, col = "black")
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(wtCol, mutantCol),
         text.col = c(wtCol, mutantCol),
         text.font = c(1, 3),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
}

# Function to plot -log10-transformed P-values derived from Fisher's exact tests
minusLog10PvalPlot <- function(xplot, Pvals) {
  plot(xplot, -log10(Pvals), type = "l", lwd = 1.5, col = "black",
       ylim = c(min(-log10(Pvals)),
                max(-log10(Pvals))),
       xlab = "",
       ylab = "",
       main = bquote("Fisher's exact tests (heterozygote versus "*italic("taf4b-1")*" crossover counts per "*.(winName)*" window)"))
  mtext(side = 2, line = 2.25, cex = 1, text = bquote("-Log"[10]*"("*italic("P")*"-value) "*.(winName)^-1), col = "black")
  abline(v = sumchr, lty = 1, lwd = 0.75)
  abline(v = centromeres, lty = 5, lwd = 0.75, col = "black")
  abline(h = -log10(0.05), lty = 5, lwd = 1, col = "red")
  box(lwd = 1.5)
}

pdf(file = paste0(inDir, "Fisher_exact_test_Pvals_genomeProfile_", winName, "_v240818.pdf"),
    height = 6, width = 9)
par(mfrow = c(2, 1))
par(mar = c(2.1, 4.1, 2.1, 4.1))
par(mgp = c(3, 1, 0))

minusLog10PvalPlot(xplot = wt$cumWindows,
                   Pvals = unlist(fisherPvals))
wtVmutantGenomePlot(xplot = wt$cumWindows,
                    wt = wt$COsChrWin,
                    mutant = mutant$COsChrWin,
                    Ylabel = bquote("Crossover density "*.(winName)^-1),
                    legendLoc = "topleft",
                    legendLabs = c("Heterozygote", "taf4b-1"),
                    wtCol = "blueviolet",
                    mutantCol = "dodgerblue3")
dev.off()

