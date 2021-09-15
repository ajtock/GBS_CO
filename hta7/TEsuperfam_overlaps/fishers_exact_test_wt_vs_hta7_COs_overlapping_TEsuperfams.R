#!/applications/R/R-4.0.0/bin/Rscript

# Compare proportions of wt and hta7 crossovers overlapping features using Fisher's exact tests

# Usage: 
# /applications/R/R-4.0.0/bin/Rscript fishers_exact_test_wt_vs_hta7_COs_overlapping_TEsuperfams.R coller.filt coller.hta7 wt hta7 "Pericentromeric crossovers" peri

#pop1Name <- "coller.filt"
#pop2Name <- "coller.hta7"
#pop1NamePlot <- "wt"
#pop2NamePlot <- "hta7"
#plotTitle <- "Pericentromeric crossovers"
#region <- "peri"

args <- commandArgs(trailingOnly = T)
pop1Name <- args[1]
pop2Name <- args[2]
pop1NamePlot <- args[3]
pop2NamePlot <- args[4]
plotTitle <- args[5]
region <- args[6]

library(GenomicRanges)
library(regioneR)
library(ggplot2)
library(ggthemes)
library(ggsignif)

# Genomic definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrStart <- c(rep(1, 5))
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# Pericentromeric regions are as defined in Supplemental Table S26
# of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))
if(region == "peri") {
  mask <- toGRanges(data.frame(rep(chrs, 2),
                               c(chrStart, pericenEnd),
                               c(pericenStart, chrLens)))
} else if(region == "arm") {
  mask <- toGRanges(data.frame(chrs,
                               pericenStart,
                               pericenEnd))
} else if(region == "genomewide") {
  mask <- GRanges()
}

plotDir <- "bar_graphs/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Load pop1COs and convert into GRanges object
pop1COs <- read.table(paste0("../COs/", pop1Name, "cos.txt"),
                      header = T)
pop1COsGR_pooled <- GRanges(seqnames = paste0("Chr",
                                              pop1COs$chr),
                            ranges = IRanges(start = pop1COs$start,
                                             end = pop1COs$end),
                            strand = "*",
                            midpoint = pop1COs$cos,
                            lib = paste0(as.character(pop1COs$lane), ".",
                                         as.character(pop1COs$lib)))
# Convert library (individual) names (e.g., from "1.1" to "1.01") to enable
# sorting individuals by increasing number
pop1COsGR_pooled$lib[grep("\\d\\.\\d$",
                          pop1COsGR_pooled$lib)] <- paste0(substr(x = pop1COsGR_pooled$lib[grep("\\d\\.\\d$",
                                                                                                pop1COsGR_pooled$lib)],
                                                                  start = 1, stop = 1),
                                                           ".0",
                                                           substr(x = pop1COsGR_pooled$lib[grep("\\d\\.\\d$",
                                                                                                pop1COsGR_pooled$lib)],
                                                                  start = 3, stop = 3))
pop1libNames <- sort(unique(pop1COsGR_pooled$lib))
print("Individuals with crossovers:")
print(length(unique(pop1COsGR_pooled$lib)))

pop1COs_mask_overlaps <- findOverlaps(query = mask,
                                      subject = pop1COsGR_pooled,
                                      ignore.strand = TRUE,
                                      select = "all")
if(length(subjectHits(pop1COs_mask_overlaps)) > 0) {
  pop1COsGR <- pop1COsGR_pooled[-subjectHits(pop1COs_mask_overlaps)]
} else {
  pop1COsGR <- pop1COsGR_pooled
}
print(paste0(pop1NamePlot, " crossovers  ", region, ":"))
print(pop1COsGR)
print(length(pop1COsGR))

# Load pop2COs and convert into GRanges object
pop2COs <- read.table(paste0("../COs/", pop2Name, "cos.txt"),
                      header = T)
pop2COsGR_pooled <- GRanges(seqnames = paste0("Chr",
                                              pop2COs$chr),
                            ranges = IRanges(start = pop2COs$start,
                                             end = pop2COs$end),
                            strand = "*",
                            midpoint = pop2COs$cos,
                            lib = paste0(as.character(pop2COs$lane), ".",
                                         as.character(pop2COs$lib)))
# Convert library (individual) names (e.g., from "1.1" to "1.01") to enable
# sorting individuals by increasing number
pop2COsGR_pooled$lib[grep("\\d\\.\\d$",
                          pop2COsGR_pooled$lib)] <- paste0(substr(x = pop2COsGR_pooled$lib[grep("\\d\\.\\d$",
                                                                                                pop2COsGR_pooled$lib)],
                                                                  start = 1, stop = 1),
                                                           ".0",
                                                           substr(x = pop2COsGR_pooled$lib[grep("\\d\\.\\d$",
                                                                                                pop2COsGR_pooled$lib)],
                                                                  start = 3, stop = 3))
pop2libNames <- sort(unique(pop2COsGR_pooled$lib))
print("Individuals with crossovers:")
print(length(unique(pop2COsGR_pooled$lib)))

pop2COs_mask_overlaps <- findOverlaps(query = mask,
                                      subject = pop2COsGR_pooled,
                                      ignore.strand = TRUE,
                                      select = "all")
if(length(subjectHits(pop2COs_mask_overlaps)) > 0) {
  pop2COsGR <- pop2COsGR_pooled[-subjectHits(pop2COs_mask_overlaps)]
} else {
  pop2COsGR <- pop2COsGR_pooled
}
print(paste0(pop2NamePlot, " crossovers  ", region, ":"))
print(pop2COsGR)
print(length(pop2COsGR))

# TEs
DNAfamNames <- c("dna", "heli", "ptmari", "mudr", "enspm", "hat", "harbinger")
RNAfamNames <- c("rna", "gypsy", "copia", "linel1", "sine")
DNAfamNamesPlot <- c("DNA", "Helitron", "Pogo/Tc1/Mariner", "MuDR", "EnSpm/CACTA", "hAT", "Harbinger")
RNAfamNamesPlot <- c("RNA", "Gypsy LTR", "Copia LTR", "LINE-1", "SINE")

DNAdir <- "/projects/ajt200/TAIR10/TE_classes/DNA/"
RNAdir <- "/projects/ajt200/TAIR10/TE_classes/RNA/"

### DNA TEs
TEsDNAGR <- lapply(seq_along(DNAfamNames), function(x) {
  TEsDNA <- read.table(paste0(DNAdir,
                              "TAIR10_Buisine_TEs_strand_tab_ann_",
                              DNAfamNames[x], ".txt"),
                       header = T)
  TEsDNAGR_all <- GRanges(seqnames = TEsDNA$chr,
                          ranges = IRanges(start = TEsDNA$start,
                                           end = TEsDNA$end),
                          strand = "*")
  maskTEsDNAoverlaps <- findOverlaps(query = mask,
                                     subject = TEsDNAGR_all,
                                     ignore.strand = TRUE,
                                     select = "all")
  if(length(subjectHits(maskTEsDNAoverlaps)) > 0) {
    TEsDNAGR_all[-subjectHits(maskTEsDNAoverlaps)]
  } else {
    TEsDNAGR_all
  }
})
### RNA TEs
TEsRNAGR <- lapply(seq_along(RNAfamNames), function(x) {
  TEsRNA <- read.table(paste0(RNAdir,
                              "TAIR10_Buisine_TEs_strand_tab_ann_",
                              RNAfamNames[x], ".txt"),
                       header = T)
  TEsRNAGR_all <- GRanges(seqnames = TEsRNA$chr,
                          ranges = IRanges(start = TEsRNA$start,
                                           end = TEsRNA$end),
                          strand = "*")
  maskTEsRNAoverlaps <- findOverlaps(query = mask,
                                     subject = TEsRNAGR_all,
                                     ignore.strand = TRUE,
                                     select = "all")
  if(length(subjectHits(maskTEsRNAoverlaps)) > 0) {
    TEsRNAGR_all[-subjectHits(maskTEsRNAoverlaps)]
  } else {
    TEsRNAGR_all
  }
})

otherNames <- c(
                DNAfamNames,
                RNAfamNames
               )
otherNamesPlot <- c(
                    DNAfamNamesPlot,
                    RNAfamNamesPlot
                   )
grl <- c(
         TEsDNAGR,
         TEsRNAGR
        )

# Function to create 2x2 contingency table where pop1COsIn and pop1COsOut are the
# number of crossovers overlapping and not overlapping a given set of features
contingencyTable <- function(pop1COsIn, pop1COsOut, pop2COsIn, pop2COsOut) {
  conTab <- matrix(c(pop1COsIn, pop1COsOut, pop2COsIn, pop2COsOut), ncol = 2)
  rownames(conTab) <- c("inFeature", "outFeature")
  colnames(conTab) <- c(pop1NamePlot, pop2NamePlot)
  conTab
}

# Count number of crossovers overlapping given sets of features
pop1COsoverlaps <- sapply(seq_along(grl), function(x) {
  sum(countOverlaps(query = pop1COsGR,
                    subject = grl[[x]],
                    ignore.strand = TRUE) > 0)
})
pop2COsoverlaps <- sapply(seq_along(grl), function(x) {
  sum(countOverlaps(query = pop2COsGR,
                    subject = grl[[x]],
                    ignore.strand = TRUE) > 0)
})

# Calculate proportion of COs overlapping given sets of features
pop1COsprop <- sapply(seq_along(pop1COsoverlaps), function(x) {
  (pop1COsoverlaps[x]/(length(pop1COsGR)-pop1COsoverlaps[x]))*100
})
pop2COsprop <- sapply(seq_along(pop2COsoverlaps), function(x) {
  (pop2COsoverlaps[x]/(length(pop2COsGR)-pop2COsoverlaps[x]))*100
})

# Perform chi-square test on each contingency table and extract P-values
fisherPvalsCOs <- sapply(1:length(grl), function(x) {
  fisher.test(x = contingencyTable(pop1COsIn = pop1COsoverlaps[x],
                                   pop1COsOut = length(pop1COsGR)-pop1COsoverlaps[x],
                                   pop2COsIn = pop2COsoverlaps[[x]], 
                                   pop2COsOut = length(pop2COsGR)-pop2COsoverlaps[x]),
              alternative = "two.sided")$p.value
})

fisherAdjPvalsCOs <- p.adjust(fisherPvalsCOs, method = "BH")

# Disable scientific notation
options(scipen = 999)

# Simplify P-values
fisherPvalsCOschar <- sapply(seq_along(fisherPvalsCOs),
  function(x) {
  if(fisherPvalsCOs[x] < 0.0001) {
    "< 0.0001"
  } else {
    paste0("= ", as.character(round(fisherPvalsCOs[x], digits = 4)))
  }
})
fisherAdjPvalsCOschar <- sapply(seq_along(fisherAdjPvalsCOs),
  function(x) {
  if(fisherAdjPvalsCOs[x] < 0.0001) {
    "< 0.0001"
  } else {
    paste0("= ", as.character(round(fisherAdjPvalsCOs[x], digits = 4)))
  }
})

# Combine in data.frame
df <- data.frame(Sample = rep(c(pop1NamePlot, pop2NamePlot),
                              each = length(grl)),
                 Feature = rep(otherNamesPlot, 2),
                 Overlap_proportion = c(pop1COsprop,
                                        pop2COsprop))

df$Feature <- factor(df$Feature,
                     levels = otherNamesPlot)
df$Sample <- factor(df$Sample,
                    levels = c(pop1NamePlot, pop2NamePlot))

bp <- ggplot(data = df,
             mapping = aes(x = Feature,
                           y = Overlap_proportion,
                           fill = Sample)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
  scale_fill_manual(name = "",
                    values = c("dodgerblue2",
                               "darkorange2"),
                    labels = c(pop1NamePlot,
                               bquote(italic(.(pop2NamePlot))))) +
  labs(y =  expression("% COs overlapping features"[ ])) +
  scale_y_continuous(limits = c(0, max(c(pop1COsprop, pop2COsprop))+10)) +
  scale_x_discrete(position = "top") +
  guides(fill = guide_legend(direction = "horizontal",
                             label.position = "top",
                             label.theme = element_text(size = 20, hjust = 0, vjust = 0.5, angle = 90),
                             reverse = TRUE,
                             nrow = 2,
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
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(colour = "transparent",
                                  fill = "transparent"),
        legend.key.size = unit(1, "cm"),
        plot.margin = unit(c(5.5, 5.5, 10.5, 5.5), "pt"),
        plot.title = element_text(size = 17, colour = "black", hjust = 0.5)) +
  ggtitle(plotTitle) +
  annotate(geom = "text",
           size = 6,
           x = c(0.8, 1.2,
                 1.8, 2.2,
                 2.8, 3.2,
                 3.8, 4.2,
                 4.8, 5.2,
                 5.8, 6.2,
                 6.8, 7.2, 
                 7.8, 8.2,
                 8.8, 9.2,
                 9.8, 10.2,
                 10.8, 11.2,
                 11.8, 12.2),
#           x = c(1.0,
#                 2.0,
#                 3.0,
#                 4.0,
#                 5.0,
#                 6.0,
#                 7.0,
#                 8.0,
#                 9.0,
#                 10.0,
#                 11.0,
#                 12.0),
           y = c(max(c(pop1COsprop, pop2COsprop))+5), angle = 90,
           label = c(paste0("P ", fisherPvalsCOschar[1]),
                     paste0("Q ", fisherAdjPvalsCOschar[1]),
                     paste0("P ", fisherPvalsCOschar[2]),
                     paste0("Q ", fisherAdjPvalsCOschar[2]),
                     paste0("P ", fisherPvalsCOschar[3]),
                     paste0("Q ", fisherAdjPvalsCOschar[3]),
                     paste0("P ", fisherPvalsCOschar[4]),
                     paste0("Q ", fisherAdjPvalsCOschar[4]),
                     paste0("P ", fisherPvalsCOschar[5]),
                     paste0("Q ", fisherAdjPvalsCOschar[5]),
                     paste0("P ", fisherPvalsCOschar[6]),
                     paste0("Q ", fisherAdjPvalsCOschar[6]),
                     paste0("P ", fisherPvalsCOschar[7]),
                     paste0("Q ", fisherAdjPvalsCOschar[7]),
                     paste0("P ", fisherPvalsCOschar[8]),
                     paste0("Q ", fisherAdjPvalsCOschar[8]),
                     paste0("P ", fisherPvalsCOschar[9]),
                     paste0("Q ", fisherAdjPvalsCOschar[9]),
                     paste0("P ", fisherPvalsCOschar[10]),
                     paste0("Q ", fisherAdjPvalsCOschar[10]),
                     paste0("P ", fisherPvalsCOschar[11]),
                     paste0("Q ", fisherAdjPvalsCOschar[11]),
                     paste0("P ", fisherPvalsCOschar[12]),
                     paste0("Q ", fisherAdjPvalsCOschar[12])))
#  geom_signif(y_position = c(70), xmin = c(0.4), xmax = c(0.8),
#              annotation = c(paste0("P = ", fisherPvalsCOs[1])),
#              tip_length = 0) 
ggsave(paste0(plotDir, "barplot_proportion_of_COs_in_",
              pop1NamePlot, "_and_", pop2NamePlot, "_overlapping_TEsuperfams",
              "_BH_adjusted_P_",
              region, ".pdf"),
       plot = bp,
       height = 8, width = 10)
