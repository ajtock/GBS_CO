#!/applications/R/R-3.5.0/bin/Rscript

# Combine ZINB regression model summaries for different window sizes
# (output of regression_SNP_count_per_win_COs_midpoints.R) into one table

# Usage from within directory containing windowed CO counts files:
# ./regression_SNP_count_per_win_COs_midpoints_rbind.R coller.filtarb coller.filtmsh2 15000 15kb

#pop1Name <- "coller.filtarb"
#pop2Name <- "coller.filtmsh2"
#flankSize <- 15000
#flankName <- "15kb"

args <- commandArgs(trailingOnly = TRUE)
pop1Name <- args[1]
pop2Name <- args[2]
flankSize <- as.numeric(args[3])
flankName <- as.character(args[4])

winNames <- c("100bp", "200bp", "500bp", "1000bp", "1500bp", "2500bp", "5000bp")
sizeDir <- paste0(flankName, "_flank_", winNames, "_win/")
inDir <- paste0(sizeDir, "regression_models/")
outDir <- paste0(flankName, "_flank_all_wins/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
outDir <- paste0(outDir, "regression_models/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))

dfList <- lapply(seq_along(winNames), function(x) {
  read.table(paste0(inDir[x], pop1Name, "_v_", pop2Name, "_SNP_frequency_",
                    "at_crossover_midpoints_ZINB_", winNames[x], "_win.tsv"),
             header = T)
})
dfrbind <- do.call(rbind, dfList)
colnames(dfrbind) <- c("Window size", "ZINB", "Estimate", "2.5% CI", "97.5% CI", "Std. error", "Z-value", "P-value")
write.table(dfrbind,
            file = paste0(outDir, pop1Name, "_v_", pop2Name, "_SNP_frequency_",
                          "at_crossover_midpoints_ZINB_all_wins.tsv"),
            row.names = F, col.names = T, sep = "\t", quote = T)
