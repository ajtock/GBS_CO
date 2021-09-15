#!/applications/R/R-3.5.0/bin/Rscript

# Plot bar graph of log2(observed:expected) peaks overlapping other features

# Usage:
# ./bargraph_permTest_COs_vs_ovule_hypermeth.R "Serra/Underwood COs" " " SerraUnderwood_COs SerraUnderwood_COs 10000

library(ggplot2)
library(ggthemes)

#plotTitle <- "Serra/Underwood COs"
#pt1Name <- " "
#pt1LibName <- "SerraUnderwood_COs" 
#ptOrderLibName <- "SerraUnderwood_COs"
## Number of permutations (randomisations) performed
#perms <- "10000"

args <- commandArgs(trailingOnly = T)
plotTitle <- args[1]
pt1Name <- args[2]
pt1LibName <- args[3]
ptOrderLibName <- args[4]
# Number of permutations (randomisations) performed
perms <- as.character(args[5])

ptOrderDir <- "./"
pt1Dir <- "./"

plotDir <- "./bar_graphs/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))


otherNames <- c(
                "loci",
                "lociUp500",
                "lociDown500",
                "lociUp1000",
                "lociDown1000",
                "lociUp2000",
                "lociDown2000"
               )
otherNamesPlot <- c(
                    "Hyper-CHH loci",
                    "500 bp upstream",
                    "500 bp downstream",
                    "1000 bp upstream",
                    "1000 bp downstream",
                    "2000 bp upstream",
                    "2000 bp downstream"
                   ) 

# Load permutation test results for peak set to be used for ordering
# of other features in bar graph
load(paste0(ptOrderDir, "permTest_", perms, "perms_",
            ptOrderLibName, "_vs_ovule_hypermeth_loci.RData"))
ptOrder <- ptLociOtherPerChrom
ptLociOtherPerChrom <- NULL

# Load permutation test results to be used for plotting
load(paste0(pt1Dir, "permTest_", perms, "perms_",
            pt1LibName, "_vs_ovule_hypermeth_loci.RData"))
pt1 <- ptLociOtherPerChrom
ptLociOtherPerChrom <- NULL

# ptOrder
ptOrder_Pval <- lapply(seq_along(ptOrder), function(x) {
  ptOrder[[x]]$numOverlaps$pval
})
ptOrder_Obs <- lapply(seq_along(ptOrder), function(x) {
  ptOrder[[x]]$numOverlaps$observed
})
ptOrder_Perm <- lapply(seq_along(ptOrder), function(x) {
  ptOrder[[x]]$numOverlaps$permuted
})
ptOrder_Exp <- lapply(seq_along(ptOrder), function(x) {
  mean(ptOrder[[x]]$numOverlaps$permuted)
})
ptOrder_log2ObsExp <- lapply(seq_along(ptOrder_Obs), function(x) {
  log2((ptOrder_Obs[[x]]+1)/(ptOrder_Exp[[x]]+1))
})
ptOrder_Zscore <- lapply(seq_along(ptOrder), function(x) {
  ptOrder[[x]]$numOverlaps$zscore
})
ptOrder_AltHyp <- lapply(seq_along(ptOrder), function(x) {
  ptOrder[[x]]$numOverlaps$alternative
})
ptOrder_alpha0.05 <- lapply(seq_along(ptOrder_Perm), function(x) {
  if(ptOrder_AltHyp[[x]] == "greater") {
    quantile(ptOrder_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(ptOrder_Perm[[x]], 0.05)[[1]]
  }
})
ptOrder_log2alpha0.05 <- lapply(seq_along(ptOrder_alpha0.05), function(x) {
  log2((ptOrder_alpha0.05[[x]]+1)/(ptOrder_Exp[[x]]+1))
})

# Order permutation test statistics and feature names by
# decreasing log2(observed/expected) overlaps of wt treatment peaks
ptOrder_log2ObsExp_sorted <- sort.int(unlist(ptOrder_log2ObsExp),
                                      decreasing = T)
ptOrder_log2alpha0.05_sorted <- unlist(ptOrder_log2alpha0.05[sort.int(unlist(ptOrder_log2ObsExp),
                                                                      decreasing = T,
                                                                      index.return = T)$ix])
ptOrder_otherNames_sorted <- otherNames[sort.int(unlist(ptOrder_log2ObsExp),
                                             decreasing = T,
                                             index.return = T)$ix]
ptOrder_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(ptOrder_log2ObsExp),
                                                     decreasing = T,
                                                     index.return = T)$ix]

# pt1
pt1_Pval <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$pval
})
pt1_Obs <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$observed
})
pt1_Perm <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$permuted
})
pt1_Exp <- lapply(seq_along(pt1), function(x) {
  mean(pt1[[x]]$numOverlaps$permuted)
})
pt1_log2ObsExp <- lapply(seq_along(pt1_Obs), function(x) {
  log2((pt1_Obs[[x]]+1)/(pt1_Exp[[x]]+1))
})
pt1_Zscore <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$zscore
})
pt1_AltHyp <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$alternative
})
pt1_alpha0.05 <- lapply(seq_along(pt1_Perm), function(x) {
  if(pt1_AltHyp[[x]] == "greater") {
    quantile(pt1_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt1_Perm[[x]], 0.05)[[1]]
  }
})
pt1_log2alpha0.05 <- lapply(seq_along(pt1_alpha0.05), function(x) {
  log2((pt1_alpha0.05[[x]]+1)/(pt1_Exp[[x]]+1))
})

# Order permutation test statistics and feature names by
# decreasing log2(observed/expected) overlaps of wt treatment peaks
pt1_log2ObsExp_sorted <- unlist(pt1_log2ObsExp[sort.int(unlist(ptOrder_log2ObsExp),
                                                        decreasing = T,
                                                        index.return = T)$ix])
pt1_log2alpha0.05_sorted <- unlist(pt1_log2alpha0.05[sort.int(unlist(ptOrder_log2ObsExp),
                                                              decreasing = T,
                                                              index.return = T)$ix])
pt1_otherNames_sorted <- otherNames[sort.int(unlist(ptOrder_log2ObsExp),
                                         decreasing = T,
                                         index.return = T)$ix]
pt1_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(ptOrder_log2ObsExp),
                                                 decreasing = T,
                                                 index.return = T)$ix]

# Combine in data.frame
df <- data.frame(Sample = rep(c(pt1Name),
                              each = length(ptOrder_log2ObsExp_sorted)),
                 Feature = rep(pt1_otherNamesPlot_sorted, 1),
                 log2ObsExp = c(pt1_log2ObsExp_sorted),
                 log2alpha0.05 = c(pt1_log2alpha0.05_sorted))

df$Feature <- factor(df$Feature,
                               levels = pt1_otherNamesPlot_sorted)
df$Sample <- factor(df$Sample,
                    levels = c(pt1Name))

bp <- ggplot(data = df,
             mapping = aes(x = Feature,
                           y = log2ObsExp,
                           fill = Sample)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
  scale_fill_manual(name = "",
                    values = c("dodgerblue2"),
                    labels = c(pt1Name)) +
  geom_point(mapping = aes(Feature, log2alpha0.05),
             position = position_dodge(0.9),
             shape = "-", colour  = "grey70", size = 20) +
  labs(y = expression("Log"[2]*"(observed/expected) CO overlap")) +
  scale_y_continuous(limits = c(-0.9, 0.9)) +
  scale_x_discrete(position = "top") +
  guides(fill = guide_legend(direction = "horizontal",
                             label.position = "top",
                             label.theme = element_text(size = 20, hjust = 0, vjust = 0.5, angle = 90),
                             nrow = 1,
                             byrow = TRUE)) +
  theme_bw() +
  theme(axis.line.y = element_line(size = 1, colour = "black"),
        axis.ticks.y = element_line(size = 1, colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black", hjust = 0.5, vjust = 0.5, angle = 90),
        axis.title.y = element_text(size = 20, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 20, colour = "black", hjust = 0, vjust = 0.5, angle = 90),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        #legend.position = c(0.05, 0.30),
        legend.position = "none",
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(colour = "transparent",
                                  fill = "transparent"),
        plot.margin = unit(c(5.5, 5.5, 10.5, 5.5), "pt"),
        plot.title = element_text(size = 17, colour = "black", hjust = 0.5)) +
  ggtitle(paste0(plotTitle, " (n = 17,077; ", prettyNum(as.character(perms),
                                                        big.mark = ",", trim = "T"),
                 " permutations)"))
ggsave(paste0(plotDir, "barplot_overlap_permTestResults_",
              as.character(perms), "perms_",
              "log2_Observed_Expected_",
              pt1LibName, ".pdf"),
       plot = bp,
       height = 10, width = 7)
save(bp,
     file = paste0(plotDir, "barplot_overlap_permTestResults_",
                   as.character(perms), "perms_",
                   "log2_Observed_Expected_",
                   pt1LibName, ".RData"))
