#!/applications/R/R-3.3.2/bin/Rscript

# Plot bar chart of log2(observed:expected) GBS-derived crossover intervals
# overlapping other features

# Usage:
# /applications/R/R-3.3.2/bin/Rscript other_features_wt_vs_msh2_ColLer_ColCvi_GBS_COs.R "Col/Ler and Col/Cvi GBS crossover intervals" HS_wt_ColLer_F2_2018 SB_msh2_ColLer_F2_2018 SB_wt_ColCvi_F2_2018 SB_msh2_ColCvi_F2_2018 "Col/Ler" "Col/Cvi" wt msh2 10000 

library(ggplot2)
library(ggthemes)

#dataName <- "Col/Ler and Col/Cvi GBS crossover intervals"
#pop1_geno1 <- "HS_wt_ColLer_F2_2018"
#pop1_geno2 <- "SB_msh2_ColLer_F2_2018"
#pop2_geno1 <- "SB_wt_ColCvi_F2_2018"
#pop2_geno2 <- "SB_msh2_ColCvi_F2_2018"
#pop1 <- "Col/Ler"
#pop2 <- "Col/Cvi"
#geno1 <- "wt"
#geno2 <- "msh2"
## Number of permutations (randomisations) performed
#perms <- 10000

args <- commandArgs(trailingOnly = T)
dataName <- as.character(args[1])
pop1_geno1 <- as.character(args[2])
pop1_geno2 <- as.character(args[3])
pop2_geno1 <- as.character(args[4])
pop2_geno2 <- as.character(args[5])
pop1 <- as.character(args[6])
pop2 <- as.character(args[7])
geno1 <- as.character(args[8])
geno2 <- as.character(args[9])
# Number of permutations (randomisations) performed
perms <- as.numeric(args[10])

plotDir <- "./plots/"

otherNames <- c("MSH4_Rep1GR", "REC8_HA_Rep1GR", "nucleRnucsGR",
                "SPO11GR", "SPO11_ChIP4GR",
                "H3K4me3GR", "H3K9me2GR",
                "genesGR", "promotersGR", "terminatorsGR",
                "TSSdownstream500GR", "TTSupstream500GR",
                "exonsGR", "intronsGR", "TEsGR")

otherNamesPlot <- c("MSH4 peaks", "REC8 peaks", "Nucleosomes",
                    "SPO11-1-oligo hotspots", "SPO11-1 ChIP peaks",
                    "H3K4me3 peaks", "H3K9me2 peaks",
                    "Genes", "Gene promoters", "Gene terminators",
                    "Gene 5' ends", "Gene 3' ends",
                    "Gene exons", "Gene introns", "Transposons")

load(paste0("permTest_", pop1_geno1, "_COs_vs_others.RData"))
pop1_geno1_pt <- ptCOsOtherPerChrom
ptCOsOtherPerChrom <- NULL
 
load(paste0("permTest_", pop1_geno2, "_COs_vs_others.RData"))
pop1_geno2_pt <- ptCOsOtherPerChrom
ptCOsOtherPerChrom <- NULL

load(paste0("permTest_", pop2_geno1, "_COs_vs_others.RData"))
pop2_geno1_pt <- ptCOsOtherPerChrom
ptCOsOtherPerChrom <- NULL

load(paste0("permTest_", pop2_geno2, "_COs_vs_others.RData"))
pop2_geno2_pt <- ptCOsOtherPerChrom
ptCOsOtherPerChrom <- NULL

# pop1_geno1_pt
pop1_geno1_pt_Pval <- lapply(seq_along(pop1_geno1_pt), function(x) {
  pop1_geno1_pt[[x]]$numOverlaps$pval
})
pop1_geno1_pt_Obs <- lapply(seq_along(pop1_geno1_pt), function(x) {
  pop1_geno1_pt[[x]]$numOverlaps$observed
})
pop1_geno1_pt_Perm <- lapply(seq_along(pop1_geno1_pt), function(x) {
  pop1_geno1_pt[[x]]$numOverlaps$permuted
})
pop1_geno1_pt_Exp <- lapply(seq_along(pop1_geno1_pt), function(x) {
  mean(pop1_geno1_pt[[x]]$numOverlaps$permuted)
})
pop1_geno1_pt_log2ObsExp <- lapply(seq_along(pop1_geno1_pt_Obs), function(x) {
  log2(pop1_geno1_pt_Obs[[x]]/pop1_geno1_pt_Exp[[x]])
})
pop1_geno1_pt_Zscore <- lapply(seq_along(pop1_geno1_pt), function(x) {
  pop1_geno1_pt[[x]]$numOverlaps$zscore
})
pop1_geno1_pt_AltHyp <- lapply(seq_along(pop1_geno1_pt), function(x) {
  pop1_geno1_pt[[x]]$numOverlaps$alternative
})
pop1_geno1_pt_alpha0.05 <- lapply(seq_along(pop1_geno1_pt_Perm), function(x) {
  if(pop1_geno1_pt_AltHyp[[x]] == "greater") {
    quantile(pop1_geno1_pt_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pop1_geno1_pt_Perm[[x]], 0.05)[[1]]
  }
})
pop1_geno1_pt_log2alpha0.05 <- lapply(seq_along(pop1_geno1_pt_alpha0.05), function(x) {
  log2(pop1_geno1_pt_alpha0.05[[x]]/pop1_geno1_pt_Exp[[x]])
})

pop1_geno1_pt_log2ObsExp_sorted <-                                     sort.int(unlist(pop1_geno1_pt_log2ObsExp), decreasing = T)
pop1_geno1_pt_log2ObsExp_sortedTest <- unlist(pop1_geno1_pt_log2ObsExp[sort.int(unlist(pop1_geno1_pt_log2ObsExp), decreasing = T, index.return = T)$ix])
try( if( identical(pop1_geno1_pt_log2ObsExp_sorted, pop1_geno1_pt_log2ObsExp_sortedTest) != TRUE )
      { stop("pop1_geno1_pt_log2ObsExp_sorted and pop1_geno1_pt_log2ObsExp_sortedTest are not identical!") }
)
pop1_geno1_pt_log2alpha0.05_sorted <- unlist(pop1_geno1_pt_log2alpha0.05[sort.int(unlist(pop1_geno1_pt_log2ObsExp), decreasing = T, index.return = T)$ix])
pop1_geno1_pt_otherNames_sorted <- otherNames[sort.int(unlist(pop1_geno1_pt_log2ObsExp), decreasing = T, index.return = T)$ix]
pop1_geno1_pt_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(pop1_geno1_pt_log2ObsExp), decreasing = T, index.return = T)$ix]

# pop1_geno2_pt
pop1_geno2_pt_Pval <- lapply(seq_along(pop1_geno2_pt), function(x) {
  pop1_geno2_pt[[x]]$numOverlaps$pval
})
pop1_geno2_pt_Obs <- lapply(seq_along(pop1_geno2_pt), function(x) {
  pop1_geno2_pt[[x]]$numOverlaps$observed
})
pop1_geno2_pt_Perm <- lapply(seq_along(pop1_geno2_pt), function(x) {
  pop1_geno2_pt[[x]]$numOverlaps$permuted
})
pop1_geno2_pt_Exp <- lapply(seq_along(pop1_geno2_pt), function(x) {
  mean(pop1_geno2_pt[[x]]$numOverlaps$permuted)
})
pop1_geno2_pt_log2ObsExp <- lapply(seq_along(pop1_geno2_pt_Obs), function(x) {
  log2(pop1_geno2_pt_Obs[[x]]/pop1_geno2_pt_Exp[[x]])
})
pop1_geno2_pt_Zscore <- lapply(seq_along(pop1_geno2_pt), function(x) {
  pop1_geno2_pt[[x]]$numOverlaps$zscore
})
pop1_geno2_pt_AltHyp <- lapply(seq_along(pop1_geno2_pt), function(x) {
  pop1_geno2_pt[[x]]$numOverlaps$alternative
})
pop1_geno2_pt_alpha0.05 <- lapply(seq_along(pop1_geno2_pt_Perm), function(x) {
  if(pop1_geno2_pt_AltHyp[[x]] == "greater") {
    quantile(pop1_geno2_pt_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pop1_geno2_pt_Perm[[x]], 0.05)[[1]]
  }
})
pop1_geno2_pt_log2alpha0.05 <- lapply(seq_along(pop1_geno2_pt_alpha0.05), function(x) {
  log2(pop1_geno2_pt_alpha0.05[[x]]/pop1_geno2_pt_Exp[[x]])
})

pop1_geno2_pt_log2ObsExp_sorted <- unlist(pop1_geno2_pt_log2ObsExp[sort.int(unlist(pop1_geno1_pt_log2ObsExp), decreasing = T, index.return = T)$ix])
pop1_geno2_pt_log2alpha0.05_sorted <- unlist(pop1_geno2_pt_log2alpha0.05[sort.int(unlist(pop1_geno1_pt_log2ObsExp), decreasing = T, index.return = T)$ix])
pop1_geno2_pt_otherNames_sorted <- otherNames[sort.int(unlist(pop1_geno1_pt_log2ObsExp), decreasing = T, index.return = T)$ix]
pop1_geno2_pt_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(pop1_geno1_pt_log2ObsExp), decreasing = T, index.return = T)$ix]

# pop2_geno1_pt
pop2_geno1_pt_Pval <- lapply(seq_along(pop2_geno1_pt), function(x) {
  pop2_geno1_pt[[x]]$numOverlaps$pval
})
pop2_geno1_pt_Obs <- lapply(seq_along(pop2_geno1_pt), function(x) {
  pop2_geno1_pt[[x]]$numOverlaps$observed
})
pop2_geno1_pt_Perm <- lapply(seq_along(pop2_geno1_pt), function(x) {
  pop2_geno1_pt[[x]]$numOverlaps$permuted
})
pop2_geno1_pt_Exp <- lapply(seq_along(pop2_geno1_pt), function(x) {
  mean(pop2_geno1_pt[[x]]$numOverlaps$permuted)
})
pop2_geno1_pt_log2ObsExp <- lapply(seq_along(pop2_geno1_pt_Obs), function(x) {
  log2(pop2_geno1_pt_Obs[[x]]/pop2_geno1_pt_Exp[[x]])
})
pop2_geno1_pt_Zscore <- lapply(seq_along(pop2_geno1_pt), function(x) {
  pop2_geno1_pt[[x]]$numOverlaps$zscore
})
pop2_geno1_pt_AltHyp <- lapply(seq_along(pop2_geno1_pt), function(x) {
  pop2_geno1_pt[[x]]$numOverlaps$alternative
})
pop2_geno1_pt_alpha0.05 <- lapply(seq_along(pop2_geno1_pt_Perm), function(x) {
  if(pop2_geno1_pt_AltHyp[[x]] == "greater") {
    quantile(pop2_geno1_pt_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pop2_geno1_pt_Perm[[x]], 0.05)[[1]]
  }
})
pop2_geno1_pt_log2alpha0.05 <- lapply(seq_along(pop2_geno1_pt_alpha0.05), function(x) {
  log2(pop2_geno1_pt_alpha0.05[[x]]/pop2_geno1_pt_Exp[[x]])
})

pop2_geno1_pt_log2ObsExp_sorted <- unlist(pop2_geno1_pt_log2ObsExp[sort.int(unlist(pop1_geno1_pt_log2ObsExp), decreasing = T, index.return = T)$ix])
pop2_geno1_pt_log2alpha0.05_sorted <- unlist(pop2_geno1_pt_log2alpha0.05[sort.int(unlist(pop1_geno1_pt_log2ObsExp), decreasing = T, index.return = T)$ix])
pop2_geno1_pt_otherNames_sorted <- otherNames[sort.int(unlist(pop1_geno1_pt_log2ObsExp), decreasing = T, index.return = T)$ix]
pop2_geno1_pt_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(pop1_geno1_pt_log2ObsExp), decreasing = T, index.return = T)$ix]

# pop2_geno2_pt
pop2_geno2_pt_Pval <- lapply(seq_along(pop2_geno2_pt), function(x) {
  pop2_geno2_pt[[x]]$numOverlaps$pval
})
pop2_geno2_pt_Obs <- lapply(seq_along(pop2_geno2_pt), function(x) {
  pop2_geno2_pt[[x]]$numOverlaps$observed
})
pop2_geno2_pt_Perm <- lapply(seq_along(pop2_geno2_pt), function(x) {
  pop2_geno2_pt[[x]]$numOverlaps$permuted
})
pop2_geno2_pt_Exp <- lapply(seq_along(pop2_geno2_pt), function(x) {
  mean(pop2_geno2_pt[[x]]$numOverlaps$permuted)
})
pop2_geno2_pt_log2ObsExp <- lapply(seq_along(pop2_geno2_pt_Obs), function(x) {
  log2(pop2_geno2_pt_Obs[[x]]/pop2_geno2_pt_Exp[[x]])
})
pop2_geno2_pt_Zscore <- lapply(seq_along(pop2_geno2_pt), function(x) {
  pop2_geno2_pt[[x]]$numOverlaps$zscore
})
pop2_geno2_pt_AltHyp <- lapply(seq_along(pop2_geno2_pt), function(x) {
  pop2_geno2_pt[[x]]$numOverlaps$alternative
})
pop2_geno2_pt_alpha0.05 <- lapply(seq_along(pop2_geno2_pt_Perm), function(x) {
  if(pop2_geno2_pt_AltHyp[[x]] == "greater") {
    quantile(pop2_geno2_pt_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pop2_geno2_pt_Perm[[x]], 0.05)[[1]]
  }
})
pop2_geno2_pt_log2alpha0.05 <- lapply(seq_along(pop2_geno2_pt_alpha0.05), function(x) {
  log2(pop2_geno2_pt_alpha0.05[[x]]/pop2_geno2_pt_Exp[[x]])
})

pop2_geno2_pt_log2ObsExp_sorted <- unlist(pop2_geno2_pt_log2ObsExp[sort.int(unlist(pop1_geno1_pt_log2ObsExp), decreasing = T, index.return = T)$ix])
pop2_geno2_pt_log2alpha0.05_sorted <- unlist(pop2_geno2_pt_log2alpha0.05[sort.int(unlist(pop1_geno1_pt_log2ObsExp), decreasing = T, index.return = T)$ix])
pop2_geno2_pt_otherNames_sorted <- otherNames[sort.int(unlist(pop1_geno1_pt_log2ObsExp), decreasing = T, index.return = T)$ix]
pop2_geno2_pt_otherNamesPlot_sorted <- otherNamesPlot[sort.int(unlist(pop1_geno1_pt_log2ObsExp), decreasing = T, index.return = T)$ix]

df <- data.frame(Sample = rep(c(pop1_geno1, pop1_geno2, pop2_geno1, pop2_geno2),
                              each = length(pop1_geno1_pt_log2ObsExp_sorted)),
                 Annotation_feature = rep(pop1_geno1_pt_otherNamesPlot_sorted, 2),
                 log2ObsExp = c(pop1_geno1_pt_log2ObsExp_sorted,
                                pop1_geno2_pt_log2ObsExp_sorted,
                                pop2_geno1_pt_log2ObsExp_sorted,
                                pop2_geno2_pt_log2ObsExp_sorted),
                 log2alpha0.05 = c(pop1_geno1_pt_log2alpha0.05_sorted,
                                   pop1_geno2_pt_log2alpha0.05_sorted,
                                   pop2_geno1_pt_log2alpha0.05_sorted,
                                   pop2_geno2_pt_log2alpha0.05_sorted))

df$Annotation_feature <- factor(df$Annotation_feature,
                                levels = pop1_geno1_pt_otherNamesPlot_sorted)
df$Sample <- factor(df$Sample,
                    levels = c(pop1_geno1, pop1_geno2,
                               pop2_geno1, pop2_geno2))

bp <- ggplot(data = df,
             mapping = aes(x = Annotation_feature,
                           y = log2ObsExp,
                           fill = Sample)) +
      geom_bar(stat = "identity",
               position = position_dodge()) +
      scale_fill_manual(name = "Population",
                        values = c("red",
                                   "red4",
                                   "blue",
                                   "navy"),
                        labels = c(bquote(.(geno1) ~ .(pop1)),
                                   bquote(italic(.(geno2)) ~ .(pop1)),
                                   bquote(.(geno1) ~ .(pop2)),
                                   bquote(italic(.(geno2)) ~ .(pop2)))) +
      geom_point(mapping = aes(Annotation_feature, log2alpha0.05),
                 position = position_dodge(0.9),
                 shape = "-", colour  = "grey80", size = 3.5) +
      labs(x = "Annotation feature",
           y = expression("Log"[2]*"(observed:expected) crossover overlap")) +
      theme_bw() +
      theme(axis.line.y = element_line(size = 0.5, colour = "black"),
            axis.ticks.y = element_line(size = 0.25, colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 7),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste0(dataName, " (", as.character(perms), " permutations)"))
ggsave(paste0(plotDir, "barplot_other_features_permTestResults_",
              as.character(perms), "perms_",
              "log2_Observed_Expected_",
              pop1_geno1, "_", pop1_geno2, "_peaks_",
              pop2_geno1, "_", pop2_geno2, "_peaks.pdf"),
       plot = bp,
       height = 5, width = 8)
save(bp,
     file = paste0(plotDir, "barplot_other_features_permTestResults_",
                   as.character(perms), "perms_",
                   "log2_Observed_Expected_",
                   pop1_geno1, "_", pop1_geno2, "_peaks_",
                   pop2_geno1, "_", pop2_geno2, "_peaks.RData"))
