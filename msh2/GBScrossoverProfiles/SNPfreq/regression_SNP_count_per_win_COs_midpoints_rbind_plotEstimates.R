#!/applications/R/R-3.5.0/bin/Rscript

# Combine ZINB regression model summaries for different populations
# (output of regression_SNP_count_per_win_COs_midpoints.R) into one table
# and plot estimates and 95% confidence intervals,
# coloured by P-value (low = red, high = blue)

# Usage from within directory containing windowed CO counts files:
# ./regression_SNP_count_per_win_COs_midpoints_rbind_plotEstimates.R 15000 15kb

#flankSize <- 15000
#flankName <- "15kb"

args <- commandArgs(trailingOnly = TRUE)
flankSize <- as.numeric(args[1])
flankName <- as.character(args[2])

library(ggplot2)
library(gridExtra)

pop1Names <- c("coller.filtarb", "colclc.filtarb", "coller.filt", "coller.filt", "coller.filt")
pop2Names <- c("coller.filtmsh2", "colclc.filtmsh2", "coller.filtrecq4a4b", "coller.filtrecqHEI10mask", "coller.filtHEI10mask")
winNames <- c("100bp", "200bp", "500bp", "1000bp", "1500bp", "2500bp", "5000bp")
sizeDir <- paste0(flankName, "_flank_", winNames, "_win/")
inDir <- paste0(sizeDir, "regression_models/")
outDir <- paste0(flankName, "_flank_all_wins/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
outDir <- paste0(outDir, "regression_models/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))

# Load all tabulated model summaries
dfList <- lapply(seq_along(winNames), function(x) {
  lapply(seq_along(pop2Names), function(y) {
    read.table(paste0(inDir[x], pop1Names[y], "_v_", pop2Names[y], "_SNP_frequency_",
                      "at_crossover_midpoints_ZINB_", winNames[x], "_win.tsv"),
               header = T)
  })
})

# Get row 2, with estimated change in mutant
dfList <- lapply(seq_along(winNames), function(x) {
  lapply(seq_along(pop2Names), function(y) {
    dfList[[x]][[y]][2,]
  })
})

# For each window size, combine model summaries for different populations into one table
dfListrbind <- lapply(seq_along(winNames), function(x) {
  do.call(rbind, dfList[[x]])
})
for(x in seq_along(winNames)) {
  colnames(dfListrbind[[x]]) <- c("Window size", "Population", "Estimate", "2.5% CI", "97.5% CI", "Std. error", "Z-value", "P-value")
  dfListrbind[[x]]$Population <- c("msh2 (Col/Ler)", "msh2 (Col/CLC)", "recq4ab (Col/Ler)", "recq4ab HEI10 (Col/Ler)", "HEI10 (Col/Ler)")
}

# Plot estimates and 95% confidence intervals, coloured by P-value
estPlotFun <- function(dataFrame) {
  ggplot(data = dataFrame,
         mapping = aes(x = Population,
                       y = Estimate)) +
  geom_errorbar(mapping = aes(ymin = `2.5% CI`,
                              ymax = `97.5% CI`),
                width = 0.5, colour = "red") +
  geom_point(shape = 21, size = 10, fill = "red") +
  labs(x = "",
       y = "Estimate") +
  theme_bw() +
  theme(axis.line.y = element_line(size = 1.0, colour = "black"),
        axis.ticks.y = element_line(size = 1.0, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 18, colour = "black"),
        axis.text.x = element_text(size = 18, colour = "black", hjust = 0, vjust = 0.5, angle = 45),
        axis.title = element_text(size = 10, colour = "black"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3, 0.9, 0.0, 0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote("Estimated change in SNP frequency in central\n" ~
                 .(as.character(dataFrame$`Window size`[1])) ~
                 "overlapping crossover midpoints"))
}

estPlots <- lapply(seq_along(dfListrbind), function(x) {
  estPlotFun(dataFrame = dfListrbind[[x]])
})
estPlotsGA <- do.call(grid.arrange,
                      c(estPlots, nrow = length(estPlots)))
ggsave(estPlotsGA,
       file = paste0(outDir, "wt_v_mutant_population_SNP_frequency_",
                     "at_crossover_midpoints_ZINB_all_wins_estimates_95CIs.pdf"),
       width = 18, height = 18*length(dfListrbind), units = "cm")



dfrbind <- do.call(rbind, dfList)
colnames(dfrbind) <- c("Window size", "ZINB", "Estimate", "2.5% CI", "97.5% CI", "Std. error", "Z-value", "P-value")
write.table(dfrbind,
            file = paste0(outDir, pop1Name, "_v_", pop2Name, "_SNP_frequency_",
                          "at_crossover_midpoints_ZINB_all_wins.tsv"),
            row.names = F, col.names = T, sep = "\t", quote = T)
