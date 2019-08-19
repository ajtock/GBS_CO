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
library(grid)
library(gridExtra)
library(extrafont)

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
  dfListrbind[[x]]$Population <- factor(c("msh2", "msh2 (Col/CLC)", "recq4ab",
                                          "recq4ab HEI10", "HEI10"),
                                        levels = c("msh2", "msh2 (Col/CLC)", "recq4ab",
                                                   "recq4ab HEI10", "HEI10"))
}

# Plot estimates and 95% confidence intervals, coloured by P-value
estPlotFun <- function(dataFrame) {
  ggplot(data = dataFrame,
         mapping = aes(x = Population,
                       y = Estimate,
                       colour = `P-value`)) +
  labs(colour = "P-value") +
  geom_errorbar(mapping = aes(ymin = `2.5% CI`,
                              ymax = `97.5% CI`),
                width = 0.2, size = 2) +
  geom_point(shape = 19, size = 10) +
  scale_color_gradient2(low = "red",
                        mid = "yellow",
                        high = "blue",
                        midpoint = 0.05) +
  scale_x_discrete(breaks = as.vector(dataFrame$Population),
                   labels = c(bquote(italic(.(as.vector(dataFrame$Population)[1]))),
                              bquote(italic(.(as.vector(dataFrame$Population)[2]))),
                              bquote(italic(.(as.vector(dataFrame$Population)[3]))),
                              bquote(italic(.(as.vector(dataFrame$Population)[4]))),
                              bquote(italic(.(as.vector(dataFrame$Population)[5]))))) +
  labs(x = "",
       y = "Estimate") +
  theme_bw() +
  theme(axis.line.y = element_line(size = 1.5, colour = "black"),
        axis.ticks.length = unit(5, "pt"),
        axis.ticks.y = element_line(size = 1.5, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 18, colour = "black"),
        axis.text.x = element_text(size = 18, colour = "black", hjust = 1.0, vjust = 1.0, angle = 45),
        axis.title = element_text(size = 18, colour = "black"),
        legend.text = element_text(size = 18, colour = "black"),
        legend.title = element_text(size = 18, colour = "black"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3, 0.9, 0.9, 0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  ggtitle(bquote(.(as.character(dataFrame$`Window size`[1])) ~
                 "centred on crossover midpoints"))
}
estPlots <- lapply(seq_along(dfListrbind), function(x) {
  estPlotFun(dataFrame = dfListrbind[[x]])
})
estPlotsGA <- do.call(grid.arrange,
                      c(estPlots, nrow = length(estPlots)))
ggsave(estPlotsGA,
       file = paste0(outDir, "wt_v_mutant_population_SNP_frequency_",
                     "at_crossover_midpoints_ZINB_all_wins_estimates_95CIs.pdf"),
       width = 20, height = 18*length(dfListrbind), units = "cm")
