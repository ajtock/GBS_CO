#!/applications/R/R-3.5.0/bin/Rscript

# Perform a Fisher's exact test on crossover counts in each genomic window
# for two genotypes
# Argument 3 specifies the number of proportionally scaled windows
# that each chromosome has been divided into

# Usage from within directory containing windowed CO counts files:
# ./3_fisher_exact_tests_scaled_win_TelCen_Emma.R EL_wt_ColBur_F2_2019 EL_taf4b_ColBur_F2_2019 10 0.1

library(GenomicRanges)

args <- commandArgs(trailingOnly = TRUE)
pop1Name <- args[1]
pop2Name <- args[2]
prop <- as.numeric(args[3])
propName <- paste0(as.character(prop), "ths")
FDR <- as.numeric(args[4])
FDRname <- paste0("FDR", as.character(FDR))

# Load crossover interval GRanges objects to determine
# number of individuals per population
load(paste0("./", pop1Name, "_pooled_GRanges.RData"))
if( exists("popCOsGR") == TRUE ) {
  popCOsGR_pooled <- popCOsGR
}
print(paste0(pop1Name, " individuals with crossovers:"))
print(length(unique(popCOsGR_pooled$lib)))
pop1Ind <- length(unique(popCOsGR_pooled$lib))
popCOsGR_pooled <- NULL
popCOsGR <- NULL

load(paste0("./", pop2Name, "_pooled_GRanges.RData"))
if( length(popCOsGR) != 0 ) {
  popCOsGR_pooled <- popCOsGR
}
print(paste0(pop2Name, " individuals with crossovers:"))
print(length(unique(popCOsGR_pooled$lib)))
pop2Ind <- length(unique(popCOsGR_pooled$lib))
popCOsGR_pooled <- NULL

# Define and create new directory and subdirectory
# to contain plotted Fisher exact test results
outDir <- paste0("./fisher_exact_tests/")
plotDir <- paste0(outDir, FDRname, "/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Load crossover interval counts in proportionally scaled windows
# and sum counts across chromosome arms
# pop1
pop1TelCenMatrix <- read.table(paste0("./", pop1Name, "_pooled_CO_frequency_",
                                      propName, "_TelCenMatrix.txt"))
pop1TelCenProfile <- as.vector(rowSums(pop1TelCenMatrix))
print(paste0(pop1Name, " TEL-CEN crossover counts matrix"))
print(pop1TelCenMatrix)
print(paste0(pop1Name, " TEL-CEN summed crossover counts profile"))
print(pop1TelCenProfile)
# pop2
pop2TelCenMatrix <- read.table(paste0("./", pop2Name, "_pooled_CO_frequency_",
                                      propName, "_TelCenMatrix.txt"))
pop2TelCenProfile <- as.vector(rowSums(pop2TelCenMatrix))
print(paste0(pop2Name, " TEL-CEN crossover counts matrix"))
print(pop2TelCenMatrix)
print(paste0(pop2Name, " TEL-CEN summed crossover counts profile"))
print(pop2TelCenProfile)

# Function to create a 2x2 contingency table where pop1InWin and pop1OutWin are the
# number of crossovers in pop1 located within and outside of a given genomic window
contingencyTable <- function(pop1InWin, pop1OutWin, pop2InWin, pop2OutWin) {
  conTab <- matrix(c(pop1InWin, pop1OutWin, pop2InWin, pop2OutWin), ncol = 2)
  rownames(conTab) <- c("InWin", "OutWin")
  colnames(conTab) <- c(pop1Name, pop2Name)
  conTab
}

# Create and save list of contingency tables
contingencyTableList <- lapply(1:length(pop1TelCenProfile), function(x) {
  contingencyTable(pop1InWin = pop1TelCenProfile[x],
                   pop1OutWin = sum(pop1TelCenProfile)-pop1TelCenProfile[x],
                   pop2InWin = pop2TelCenProfile[x],
                   pop2OutWin = sum(pop2TelCenProfile)-pop2TelCenProfile[x]) 
})
save(contingencyTableList,
     file = paste0(outDir,
                   "contingencyTableList_TelCenProfiles_",
                   pop1Name, "_vs_", pop2Name, "_", propName, "_", FDRname, ".RData"))

# Perform Fisher's exact test on each contingency table and extract P-values
fisherPvals <- sapply(1:length(pop1TelCenProfile), function(x) {
  fisher.test(contingencyTable(pop1InWin = pop1TelCenProfile[x],
                               pop1OutWin = sum(pop1TelCenProfile)-pop1TelCenProfile[x],
                               pop2InWin = pop2TelCenProfile[x],
                               pop2OutWin = sum(pop2TelCenProfile)-pop2TelCenProfile[x]),
              alternative = "two.sided")$p.value
})

# Adjust P-values by applying the Benjamini-Hochberg multiple testing correction
fisherAdjPvals <- p.adjust(p = fisherPvals, method = "BH")

# Function to plot -log10-transformed P-values and adjusted P-values derived from Fisher's exact tests
minusLog10PvalPlot <- function(xplot,
                               Pvals, PvalsCol,
                               AdjPvals, AdjPvalsCol) {
  plot(x = xplot, y = -log10(Pvals), col = PvalsCol, type = "l", lwd = 1.5,
       ylim = c(0,
                pmax(-log10(0.05), max(-log10(Pvals)))),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = "")
  mtext(side = 2, line = 2.25, cex = 1, col = PvalsCol,
        text = bquote("-Log"[10]*"("*italic("P")*"-value)"))
  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
  abline(h = -log10(0.05), lty = 5, lwd = 1, col = PvalsCol)

  par(new = T)
  plot(x = xplot, y = -log10(AdjPvals), col = AdjPvalsCol, type = "l", lwd = 1.5,
       ylim = c(0,
                pmax(-log10(FDR), max(-log10(AdjPvals)))),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = bquote(.(pop1Name)*" versus "*.(pop2Name)*" crossover counts"),
       cex.main = 1)
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 1, adj = c(0.5, -3.0), xpd = NA, srt = -90, col = AdjPvalsCol,
       labels = bquote("-Log"[10]*"(BH-adjusted "*italic("P")*"-value)"))
  axis(side = 4, cex.axis = 1, lwd.tick = 1.5)
  if(prop <= 15) {
    axis(side = 1, cex.axis = 1, lwd.tick = 1.5,
         at = c(1, seq(2, prop, by = 1)),
         labels = c(expression(italic("TEL")),
                    seq(2, prop-1, by = 1),
                    expression(italic("CEN"))))
  } else {
    axis(side = 1, cex.axis = 1, lwd.tick = 1.5,
         at = c(1, pretty(1:length(xplot))[2:(length(pretty(1:length(xplot)))-1)], prop),
         labels = c(expression(italic("TEL")),
                    pretty(1:length(xplot))[2:(length(pretty(1:length(xplot)))-1)],
                    expression(italic("CEN"))))
  }
  mtext(side = 1, line = 2.25, cex = 1, text = paste0("Scaled windows (", propName, ")"))
  abline(h = -log10(FDR), lty = 5, lwd = 1, col = AdjPvalsCol)

  box(lwd = 1.5)
}

# Function to plot Tel to Cen profile of windowed pop1 vs pop2 crossovers (one Y-axis) 
pop1Vpop2GenomePlot <- function(xplot,
                                pop1, pop1Col,
                                pop2, pop2Col,
                                Ylabel,
                                legendLoc,
                                legendLabs) {
  plot(xplot, pop1, type = "l", lwd = 1.5, col = pop1Col,
       ylim = c(0,
                max(c(pop1, pop2))),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = bquote(italic("r"[s]) ~ " = " ~ .(round(cor(pop1, pop2, method = "spearman", use = "pairwise.complete.obs"), digits = 2))))
  lines(xplot, pop2, type = "l", lwd = 1.5, col = pop2Col)
  if(prop <= 15) {
    axis(side = 1, cex.axis = 1, lwd.tick = 1.5,
         at = c(1, seq(2, prop, by = 1)),
         labels = c(expression(italic("TEL")),
                    seq(2, prop-1, by = 1),
                    expression(italic("CEN"))))
  } else {
    axis(side = 1, cex.axis = 1, lwd.tick = 1.5,
         at = c(1, pretty(1:length(xplot))[2:(length(pretty(1:length(xplot)))-1)], prop),
         labels = c(expression(italic("TEL")),
                    pretty(1:length(xplot))[2:(length(pretty(1:length(xplot)))-1)],
                    expression(italic("CEN"))))
  }
  mtext(side = 1, line = 2.25, cex = 1, text = paste0("Scaled windows (", propName, ")"))
  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
  mtext(side = 2, line = 2.25, cex = 1, text = Ylabel, col = "black")
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(pop1Col, pop2Col),
         text.col = c(pop1Col, pop2Col),
         text.font = c(1, 1),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
}

pdf(file = paste0(plotDir, "Fisher_exact_test_Pvals_TelCenProfiles_", pop1Name, "_vs_", pop2Name, "_", propName, "_", FDRname, ".pdf"),
    height = 7, width = 7)
par(mfrow = c(2, 1))
par(mar = c(4.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))

# Plot P-values
minusLog10PvalPlot(xplot = 1:length(pop1TelCenProfile),
                   Pvals = fisherPvals,
                   AdjPvals = fisherAdjPvals,
                   PvalsCol = "dodgerblue3",
                   AdjPvalsCol = "red")
# Plot Tel to Cen profile of windowed pop1 vs pop2 crossovers
# normalised by total pop1 and pop2 crossovers, respectively 
pop1Vpop2GenomePlot(xplot = 1:length(pop1TelCenProfile),
                    pop1 = pop1TelCenProfile/(sum(pop1TelCenProfile)),
                    pop2 = pop2TelCenProfile/(sum(pop2TelCenProfile)),
                    Ylabel = "Crossovers/population crossovers",
                    legendLoc = "top",
                    legendLabs = c(pop1Name, pop2Name),
                    pop1Col = "grey50",
                    pop2Col = "forestgreen")
dev.off()
print(paste0("Fisher_exact_test_Pvals_TelCenProfiles_", pop1Name, "_vs_", pop2Name, "_", propName, "_", FDRname, ".pdf written to ", plotDir))

