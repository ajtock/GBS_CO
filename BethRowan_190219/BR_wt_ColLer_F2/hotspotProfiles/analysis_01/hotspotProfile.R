#!/applications/R/R-3.5.0/bin/Rscript

# Plot coverage profiles around crossover hotspots

# Usage:
# ./hotspotProfile.R 1000 1kb '1 kb'

args <- commandArgs(trailingOnly = T)
flankSize <- as.numeric(args[1])
flankName <- args[2]
flankNamePlot <- args[3]

library(parallel)

hotspots <- read.table("./hotspots_Matt.txt", header = T)
hotspots$chr <- paste0("Chr", hotspots$chr)

covListNames <- c(
                  "H3K4me3",
                  "H3K9me2",
                  "H3K27me1",
                  "H3K27me3",
                  "Nucleosomes",
                  "REC8-HA",
                  "SPO11-1-oligos"
                 )
covColours <- c(
                "forestgreen",
                "magenta3",
                "yellow4",
                "navy",
                "darkcyan",
                "green2",
                "red"
               )
makeTransparent <- function(thisColour, alpha = 150)
{
  newColour <- col2rgb(thisColour)
  apply(newColour, 2, function(x) {
    rgb(red = x[1], green = x[2], blue = x[3],
        alpha = alpha, maxColorValue = 255)
  })
}
covColours <- makeTransparent(covColours)

H3K4me3 <- read.table("/projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/WT_H3K4me3_ChIP14_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed",
                      colClasses = c(NA, NA, "NULL", NA))
H3K9me2 <- read.table("/projects/ajt200/BAM_masters/H3K9me2/WT/coverage/log2ChIPinput/log2_WT_H3K9me2_ChIP_WT_H3K9me2_input_norm_allchrs_coverage_coord_tab.bed",
                      colClass = c(NA, NA, "NULL", NA))
H3K27me1 <- read.table("/home/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K27me1/coverage/log2ChIPinput/WT_H3K27me1_Rep1_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed",
                       colClass = c(NA, NA, "NULL", NA))
H3K27me3 <- read.table("/home/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K27me3/coverage/log2ChIPinput/WT_H3K27me3_Rep1_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed",
                       colClass = c(NA, NA, "NULL", NA))
Nucleosomes <- read.table("/projects/ajt200/BAM_masters/nucleosomes/WT/coverage/nakedDNA_untrimmed_input/log2ChIPinput/log2wtNucNakedDNAuntrimmed_norm_allchrs_coverage_coord_tab.bed",
                          colClass = c(NA, NA, "NULL", NA))
REC8_HA <- read.table("/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/log2_REC8_HA_Rep1_ChIP_REC8_MYC_Rep2_input_norm_allchrs_coverage_coord_tab.bed",
                      colClass = c(NA, NA, "NULL", NA))
SPO11oligos <- read.table("/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/log2wtSPO11oligoRPI1NakedDNA_norm_allchrs_coverage_coord_tab.bed",
                          colClass = c(NA, NA, "NULL", NA))

covList <- list(
                H3K4me3,
                H3K9me2,
                H3K27me1,
                H3K27me3,
                Nucleosomes,
                REC8_HA,
                SPO11oligos
               )

lociCov <- lapply(seq_along(1:dim(hotspots)[1]), function(x) {
  mclapply(seq_along(covList), function(y) {
    covList[[y]][covList[[y]]$V1 == hotspots[x,]$chr
               & covList[[y]]$V2 >= hotspots[x,]$start-flankSize
               & covList[[y]]$V2 <= hotspots[x,]$end+flankSize,]
  }, mc.cores = length(covList))
})

save(lociCov,
     file = "./crossover_hotspot_profiles_list.RData")

# Function to overlay coverage profiles around a hotspot
locCovPlot <- function(xplot,
                       title,
                       dat1, dat2, dat3, dat4, dat5, dat6, dat7,
                       covColours,
                       flankSize,
                       flankLabL,
                       flankLabR,
                       featureStartLab,
                       featureEndLab,
                       Ylab,
                       legendLabs,
                       legendLoc) {
  plot(xplot, dat1, col = covColours[1], type = "l", lwd = 1.5,
       ylim = c(min(dat1, dat2, dat3, dat4, dat5, dat6, dat7),
                max(dat1, dat2, dat3, dat4, dat5, dat6, dat7)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = title, cex.main = 1)
  lines(xplot, dat2, col = covColours[2], type = "l", lwd = 1.5)
  lines(xplot, dat3, col = covColours[3], type = "l", lwd = 1.5)
  lines(xplot, dat4, col = covColours[4], type = "l", lwd = 1.5)
  lines(xplot, dat5, col = covColours[5], type = "l", lwd = 1.5)
  lines(xplot, dat6, col = covColours[6], type = "l", lwd = 1.5)
  lines(xplot, dat7, col = covColours[7], type = "l", lwd = 1.5)
  mtext(side = 2, line = 2.25, cex = 0.8, text = Ylab, col = "black")
  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
  axis(side = 1, cex.axis = 1, lwd.tick = 1.5,
       labels = c(flankLabL,
                  featureStartLab,
                  featureEndLab,
                  flankLabR),
       at = c(1,
              flankSize+1,
              length(dat1)-flankSize,
              length(dat1)))
  abline(v = c(flankSize+1,
               length(dat1)-flankSize),
         lty = 3, lwd = 3)
  legend(legendLoc,
         legend = legendLabs,
         col = covColours,
         text.col = covColours,
         text.font = c(1),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
  box(lwd = 1.5)
}
                       
pdf(paste0("Crossover_hotspot_profiles_", flankName, "flank.pdf"),
    height = 2.5*(dim(hotspots)[1]), width = 6)
par(mfrow = c(dim(hotspots)[1], 1))
par(mar = c(2.1, 4.2, 2.1, 2.1))
par(mgp = c(2.25, 1, 0))
for(x in 1:(dim(hotspots)[1])) {
locCovPlot(xplot = 1:length(lociCov[[x]][[1]]$V1),
           title = as.character(hotspots$sampleid)[x],
           dat1 = lociCov[[x]][[1]]$V4,
           dat2 = lociCov[[x]][[2]]$V4,
           dat3 = lociCov[[x]][[3]]$V4,
           dat4 = lociCov[[x]][[4]]$V4,
           dat5 = lociCov[[x]][[5]]$V4,
           dat6 = lociCov[[x]][[6]]$V4,
           dat7 = lociCov[[x]][[7]]$V4,
           covColours = covColours,
           flankSize = flankSize,
           flankLabL = paste0("-", flankNamePlot),
           flankLabR = paste0("+", flankNamePlot),
           featureStartLab = "Start", 
           featureEndLab = "End",
           Ylab = bquote("Z-standardized log"[2]*"(ChIP/control)"),
           legendLabs = covListNames,
           legendLoc = "topleft")
}
dev.off()
 
