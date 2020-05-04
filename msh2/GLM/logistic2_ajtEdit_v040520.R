#/applications/R/R-3.5.0/bin/Rscript

#
# author: Andy Tock (adapted from original script logistic2.R by Tom Hardcastle)
# contact: ajt200@cam.ac.uk
# date: 04.05.2020
# 

# Build GLM with GBS-derived Col/Ler crossovers as the response variable and
# SPO11-1-oligos, MNase-seq-derived nucleosome occupancy, H3K4me3 ChIP-seq,
# BS-seq-derived DNA methylation, Col/Ler SNPs, and interval width as pedictor variables

library(segmentSeq)
library(glm2)

# Genomic definitions
genomeFAI <- read.table("/projects/ajt200/TAIR10/TAIR10_chr_all.fa.fai",
                        colClasses = c(rep(NA, 2), rep("NULL", 3)), header = F)[1:5,]
chrs <- paste0("Chr", genomeFAI$V1)
chrLens <- genomeFAI$V2

# Load SNPs (481252 SNPs in BC.complete.tiger.txt)
# and create a GRanges object of 481257 intervals within which a crossover can be detected
# There are 5 more intervals than SNPs due to the interval
# between the last SNP in each chromosome and the end of that chromosome
SNPs <- read.table("BC.complete.tiger.txt",
                   header = T)
print(dim(SNPs))
#[1] 481252      6
SNPsGR <- do.call("c", lapply(1:5, function(chr) {
  chrSNPs <- SNPs[SNPs[,1] == chr,]
  GRanges(seqnames = chrs[chr],
          ranges = IRanges(start = c(1, chrSNPs[,2]),
                           end = c(chrSNPs[,2], chrLens[chr])))
}))
print(SNPsGR)

# Load 3320 wild type crossovers
COsBED <- read.table("/projects/ajt200/GBS_CO/HS_CU_080617/wt/GBS_crossovers.bed",
                     header = F)
# Note addition of 1 to start coordinates as these are in BED format
# (i.e., 0-based start coordinates)
# NOTE: Tom merged overlapping COs (i.e., reduce(COsGR), which gives 3115 crossovers) 
COsGR <- GRanges(seqnames = COsBED$V1,
                 ranges = IRanges(start = COsBED$V2+1,
                                  end = COsBED$V3))
COsGR <- reduce(COsGR)

# Identify which SNP intervals overlap at least one crossover interval
# Tom set overlapType to "within", meaning that called overlaps correspond to
# SNP intervals that are completely contained by at least one crossover interval
COevents <- getOverlaps(coordinates = SNPsGR,
                        segments = COsGR,
                        overlapType = "within",
                        whichOverlaps = FALSE,
                        ignoreStrand = TRUE)
# Number of SNP intervals completely contained by at least one crossover interval
print(sum(COevents))
#[1] 4484
# Number of SNP intervals not completely contained by at least one crossover interval
print(sum(!COevents))
#[1] 476773

# Define regionsGR object as SNP intervals that are
# not completely contained by at least one crossover interval
# + all crossover intervals

regionsGR <- c(SNPsGR[!COevents], COsGR)
regionsGR <- sort(regionsGR)
print(length(regionsGR))
#[1] 479888

# Load SPO11-1-oligos, nuclesome occupancy, H3K4me3 ChIP-seq and DNA methylation data
# and calculate sum and mean in each region

SPO11 <- read.table(paste0("/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/noZscore/",
                           "log2wtSPO11oligoRPI1NakedDNA_norm_allchrs_coverage_coord_tab_noZscore.bed"),
                    colClasses = c(NA, rep("NULL", 2), NA),
                    header = F)
splitSPO11 <- list()
for(x in 1:5) {
  chr_SPO11 <- SPO11[SPO11[,1] == chrs[x],]$V4
  chr_regionsGR <- regionsGR[seqnames(regionsGR) == chrs[x]]
  splitSPO11 <- c(splitSPO11,
                  split(x = chr_SPO11,
                        f = c(rep(1:length(chr_regionsGR),
                                  times = width(chr_regionsGR)-1),
                              length(chr_regionsGR))))
}
meanSPO11 <- sapply(splitSPO11, mean)
sumSPO11 <- sapply(splitSPO11, sum)

MNase <- read.table(paste0("/projects/ajt200/BAM_masters/nucleosomes/WT/coverage/nakedDNA_untrimmed_input/",
                           "log2ChIPinput/noZscore/log2wtNucNakedDNAuntrimmed_norm_allchrs_coverage_coord_tab_noZscore.bed"),
                    colClasses = c(NA, rep("NULL", 2), NA),
                    header = F)
splitMNase <- list()
for(x in 1:5) {
  chr_MNase <- MNase[MNase[,1] == chrs[x],]$V4
  chr_regionsGR <- regionsGR[seqnames(regionsGR) == chrs[x]]
  splitMNase <- c(splitMNase,
                  split(x = chr_MNase,
                        f = c(rep(1:length(chr_regionsGR),
                                  times = width(chr_regionsGR)-1),
                              length(chr_regionsGR))))
}
meanMNase <- sapply(splitMNase, mean)
sumMNase <- sapply(splitMNase, sum)

H3K4me3 <- read.table(paste0("/projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/noZscore/",
                           "WT_H3K4me3_ChIP14_log2ChIPinput_norm_allchrs_coverage_coord_tab_noZscore.bed"),
                    colClasses = c(NA, rep("NULL", 2), NA),
                    header = F)
splitH3K4me3 <- list()
for(x in 1:5) {
  chr_H3K4me3 <- H3K4me3[H3K4me3[,1] == chrs[x],]$V4
  chr_regionsGR <- regionsGR[seqnames(regionsGR) == chrs[x]]
  splitH3K4me3 <- c(splitH3K4me3,
                  split(x = chr_H3K4me3,
                        f = c(rep(1:length(chr_regionsGR),
                                  times = width(chr_regionsGR)-1),
                              length(chr_regionsGR))))
}
meanH3K4me3 <- sapply(splitH3K4me3, mean)
sumH3K4me3 <- sapply(splitH3K4me3, sum)

# DNA methylation
CG <- read.table(paste0("/home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/",
                        "GSM980986_WT_rep2_CG.wig.bed.gr.tab.bed"),
                 colClasses = c(rep(NA, 2), "NULL", NA),
                 header = F)
CHG <- read.table(paste0("/home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/",
                         "GSM980986_WT_rep2_CHG.wig.bed.gr.tab.bed"),
                  colClasses = c(rep(NA, 2), "NULL", NA),
                  header = F)
CHH <- read.table(paste0("/home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/",
                         "GSM980986_WT_rep2_CHH.wig.bed.gr.tab.bed"),
                  colClasses = c(rep(NA, 2), "NULL", NA),
                  header = F)
meanCG <- NULL
meanCHG <- NULL
meanCHH <- NULL
meanDNAmeth <- NULL
for(x in 1:5) {
  windowsGR <- regionsGR[seqnames(regionsGR) == chrs[x]]
  windowsGR <- GRanges(seqnames = chrs[x],
                       ranges = IRanges(start = start(windowsGR),
                                        end = end(windowsGR)-1))
  end(windowsGR[length(windowsGR)]) <- chrLens[x]
  # Define DNA methylation coordinates as GRanges objects
  # and calculate mean methylation proportions in each window
  # CG
  chr_CG <- CG[CG[,1] == paste0("chr", x),]
  chr_CG_GR <- GRanges(seqnames = chrs[x],
                       ranges = IRanges(start = chr_CG[,2],
                                        width = 1),
                       strand = "*")
  CGoverlaps <- getOverlaps(coordinates = windowsGR,
                            segments = chr_CG_GR,
                            overlapType = "overlapping",
                            whichOverlaps = TRUE,
                            ignoreStrand = TRUE)
  CGwinVals <- sapply(CGoverlaps, function(x) mean(as.numeric(chr_CG[,3][x])))
  # CHG
  chr_CHG <- CHG[CHG[,1] == paste0("chr", x),]
  chr_CHG_GR <- GRanges(seqnames = chrs[x],
                        ranges = IRanges(start = chr_CHG[,2],
                                         width = 1),
                        strand = "*")
  CHGoverlaps <- getOverlaps(coordinates = windowsGR,
                             segments = chr_CHG_GR,
                             overlapType = "overlapping",
                             whichOverlaps = TRUE,
                             ignoreStrand = TRUE)
  CHGwinVals <- sapply(CHGoverlaps, function(x) mean(as.numeric(chr_CHG[,3][x])))
  # CHH
  chr_CHH <- CHH[CHH[,1] == paste0("chr", x),]
  chr_CHH_GR <- GRanges(seqnames = chrs[x],
                        ranges = IRanges(start = chr_CHH[,2],
                                         width = 1),
                        strand = "*")
  CHHoverlaps <- getOverlaps(coordinates = windowsGR,
                             segments = chr_CHH_GR,
                             overlapType = "overlapping",
                             whichOverlaps = TRUE,
                             ignoreStrand = TRUE)
  CHHwinVals <- sapply(CHHoverlaps, function(x) mean(as.numeric(chr_CHH[,3][x])))
  # Mean of all 3 contexts
  CwinVals <- sapply(seq_along(CGwinVals), function(x) {
                mean(c(CGwinVals[x], CHGwinVals[x], CHHwinVals[x]))
              })
  # Combine values from all chromosomes
  meanCG <- c(meanCG, CGwinVals)
  meanCHG <- c(meanCHG, CHGwinVals)
  meanCHH <- c(meanCHH, CHHwinVals)
  meanDNAmeth <- c(meanDNAmeth, CwinVals)
}

# Load larger Col/Ler SNPs data set (519834 SNPs) for calculating SNP frequency
# within each interval in regionsGR
newSNPs <- read.table("/projects/ajt200/GBS_CO/msh2/SNPs/collerF2.complete.tiger.txt",
                      header = T)[,1:2]
colnames(newSNPs) <- c("chr", "pos")
print(dim(newSNPs))
#[1] 519834      6
SNPs <- SNPs[,1:2]
colnames(SNPs) <- c("chr", "pos")
print(dim(SNPs))
#[1] 481252      2
allSNPs <- rbind(SNPs[,1:2], newSNPs[,1:2]) 
print(dim(allSNPs))
#[1] 1001086       2
allSNPs <- unique(allSNPs)
print(dim(allSNPs))
#[1] 534775      2
allSNPsGR <- GRanges(seqnames = paste0("Chr", allSNPs[,1]),
                     ranges = IRanges(start = allSNPs[,2],
                                      width = 1))
meanSNPs <- NULL
sumSNPs <- NULL
for(x in 1:5) {
  windowsGR <- regionsGR[seqnames(regionsGR) == chrs[x]]
  chr_allSNPsGR <- allSNPsGR[seqnames(allSNPsGR) == chrs[x]]
  winSNPs <- countOverlaps(query = windowsGR,
                           subject = chr_allSNPsGR,
                           type = "any",
                           ignore.strand = TRUE)
  meanSNPs <- c(meanSNPs, winSNPs/width(windowsGR))
  sumSNPs <- c(sumSNPs, winSNPs)
} 
 
# Create data object for model
dat <- cbind.data.frame(CO = regionsGR %in% COsGR,
                        meanSPO11 = meanSPO11,
                        meanMNase = meanMNase,
                        meanH3K4me3 = meanH3K4me3,
                        meanDNAmeth = meanDNAmeth,
                        meanSNPs = meanSNPs,
                        sumSPO11 = sumSPO11,
                        sumMNase = sumMNase,
                        sumH3K4me3 = sumH3K4me3,
                        sumSNPs = sumSNPs,
                        width = width(regionsGR))
save(dat, file = "df_for_GLM.RData")

# Build binomial GLM with "logit" link function
#############################################################
# run binomial GLM model, link function is logit, plot data #
#############################################################
# predict function produces predicted values 
# Tom uses sumspo - why not meanspo - because width is a variable?

glmCO <- glm2(formula = CO ~ sumSPO11 * sumMNase * sumH3K4me3 * meanDNAmeth * width,
              family = binomial(link="logit"),
              data = dat)
#Warning messages:
#1: step size truncated due to increasing deviance
#2: glm.fit2: fitted probabilities numerically 0 or 1 occurred

 
glmCO <- glm(co~band+sumREC8_ChIP*(gene+promoter+terminator+RC_Helitron+DNA_Pogo+DNA_Tc1+DNA_Mariner+DNA_MuDR+DNA_Harbinger+DNA_HAT+LTR_Copia+DNA_En.Spm+SINE+LINE_L1+LTR_Gypsy+width), family = binomial(link = "logit"), data = dat, control = list(maxit = 50))
predict1 <- predict(glmCO, type = "response")
glmSum <- summary(glmCO)
coeff <- glmSum$coefficients
save(glmSum, file = paste0(outDir, "GLM_REC8_ChIP_summary.RData"))
write.csv(coeff, file = paste0(outDir, "GLM_REC8_ChIP_coeffs.csv"))
save(predict1, file = paste0(outDir, "GLM_REC8_ChIP_predict.RData"))

glmCO <- glm(co~band+sumnuc*(gene+promoter+terminator+RC_Helitron+DNA_Pogo+DNA_Tc1+DNA_Mariner+DNA_MuDR+DNA_Harbinger+DNA_HAT+LTR_Copia+DNA_En.Spm+SINE+LINE_L1+LTR_Gypsy+width), family = binomial(link = "logit"), data = dat, control = list(maxit = 50))
predict2 <- predict(glmCO, type = "response")
glmSum <- summary(glmCO)
coeff <- glmSum$coefficients
save(glmSum, file = paste0(outDir, "GLM_nuc_summary.RData"))
write.csv(coeff, file = paste0(outDir, "GLM_nuc_coeffs.csv"))
save(predict2, file = paste0(outDir, "GLM_nuc_predict.RData"))

glmCO <- glm(co~band+sumspo*(gene+promoter+terminator+RC_Helitron+DNA_Pogo+DNA_Tc1+DNA_Mariner+DNA_MuDR+DNA_Harbinger+DNA_HAT+LTR_Copia+DNA_En.Spm+SINE+LINE_L1+LTR_Gypsy+width), family = binomial(link = "logit"), data = dat, control = list(maxit = 50))
predict3 <- predict(glmCO, type = "response")
glmSum <- summary(glmCO)
coeff <- glmSum$coefficients
save(glmSum, file = paste0(outDir, "GLM_SPO11-oligos_summary.RData"))
write.csv(coeff, file = paste0(outDir, "GLM_SPO11-oligos_coeffs.csv"))
save(predict3, file = paste0(outDir, "GLM_SPO11-oligos_predict.RData"))

glmCO <- glm(co~band+sumMSH4_ChIP*(gene+promoter+terminator+RC_Helitron+DNA_Pogo+DNA_Tc1+DNA_Mariner+DNA_MuDR+DNA_Harbinger+DNA_HAT+LTR_Copia+DNA_En.Spm+SINE+LINE_L1+LTR_Gypsy+width), family = binomial(link = "logit"), data = dat, control = list(maxit = 50))
predict4 <- predict(glmCO, type = "response")
glmSum <- summary(glmCO)
coeff <- glmSum$coefficients
save(glmSum, file = paste0(outDir, "GLM_MSH4_ChIP_summary.RData"))
write.csv(coeff, file = paste0(outDir, "GLM_MSH4_ChIP_coeffs.csv"))
save(predict4, file = paste0(outDir, "GLM_MSH4_ChIP_predict.RData"))

glmCO <- glm(co~band+(sumREC8_ChIP+sumMSH4_ChIP+sumnuc+sumspo)*(gene+promoter+terminator+RC_Helitron+DNA_Pogo+DNA_Tc1+DNA_Mariner+DNA_MuDR+DNA_Harbinger+DNA_HAT+LTR_Copia+DNA_En.Spm+SINE+LINE_L1+LTR_Gypsy+width), family = binomial(link = "logit"), data = dat, control = list(maxit = 50))
predict5 <- predict(glmCO, type = "response")
glmSum <- summary(glmCO)
coeff <- glmSum$coefficients
save(glmSum, file = paste0(outDir, "GLM_REC8_MSH4_nuc_SPO11_summary.RData"))
write.csv(coeff, file = paste0(outDir, "GLM_REC8_MSH4_nuc_SPO11_coeffs.csv"))
save(predict5, file = paste0(outDir, "GLM_REC8_MSH4_nuc_SPO11_predict.RData"))

###################################################
# plot predicted CO overlaps for annotation class #
###################################################

#regid <- c(list(geneOverlap, promoterOverlap, terminatorOverlap), lapply(colnames(tesOverlap), function(tesn){tesOverlap[,tesn]}))
#names(regid) <- c("gene", "promoter", "terminator", colnames(tesOverlap))
#pdf(file = paste0(outDir, "GLM_REC8_ChIP_annotation_CO_boxplot.pdf"))
#boxplot(lapply(regid, function(x) sapply(1:100, function(ii) sum(rbinom(sum(x), 1, prob = predict1[x]))/sum(dat$width[x])*1e6)), main = "predicted CO/Mb by annotation")
#tco <- sapply(regid, function(x) sum(dat$co[x])/sum(dat$width[x])*1e6)
#points(x = 1:length(regid), y = tco, col = "red", pch = 19)
#dev.off()

#####################################################
# plot predicted CO for regionsGR by REC8_ChIP hexile #
#####################################################

#sort(dat$meanREC8_ChIP)[round(0:7*(nrow(dat)/7))]
#[1]    4.384697    6.054059    7.394007    8.744343   10.205637   12.277201
#[7] 3242.308239
#levels(cut(dat$meanREC8_ChIP, breaks = sort(dat$meanREC8_ChIP)[round(0:7*(nrow(dat)/7))]))
#[1] "(4.38,6.05]"     "(6.05,7.39]"     "(7.39,8.74]"     "(8.74,10.2]"    
#[5] "(10.2,12.3]"     "(12.3,3.24e+03]"
#length(dat$meanREC8_ChIP[dat$meanREC8_ChIP < 4.384697])
#[1] 68613
# omits 68613 of 480348 regionsGR for which meanREC8_ChIP < 4.384697
# so changed below from breaks = sort(dat$meanREC8_ChIP)[round(0:7*(nrow(dat)/7))] to breaks = c(0, sort(dat$meanREC8_ChIP)[round(0:6*(nrow(dat)/6))]) 
# so that...
#levels(cut(dat$meanREC8_ChIP, breaks = c(0, sort(dat$meanREC8_ChIP)[round(0:6*(nrow(dat)/6))])))
#[1] "(0,4.75]"        "(4.75,6.52]"     "(6.52,8.07]"     "(8.07,9.68]"    
#[5] "(9.68,11.8]"     "(11.8,3.24e+03]"
ssID <- split(1:nrow(dat), cut(dat$meanREC8_ChIP, breaks = c(0, sort(dat$meanREC8_ChIP)[round(0:6*(nrow(dat)/6))])))
pdf(file = paste0(outDir, "GLM_REC8_ChIP_hexiles_CO_boxplot.pdf"))
boxplot(lapply(ssID, function(x) sapply(1:100, function(ii) sum(rbinom(length(x), 1, prob = predict1[x]))/sum(dat$width[x])*1e6)), main = "GLM: CO~band+REC8 ChIP*(annotation+width)", xlab = "REC8 ChIP hexiles", ylab = "Number of crossovers per Mb", names = c("1", "2", "3", "4", "5", "6"))
tco <- sapply(ssID, function(x) sum(dat$co[x])/sum(dat$width[x])*1e6)
points(x = 1:length(ssID), y = tco, col = "red", pch = 19)
legend("top", legend = c("Observed", "Predicted"), text.col = c("red", "black"), bty = "n")
dev.off()

ssID <- split(1:nrow(dat), cut(dat$meanREC8_ChIP, breaks = sort(dat$meanREC8_ChIP)[round(0:7*(nrow(dat)/7))]))
pdf(file = paste0(outDir, "GLM_REC8_ChIP_hexiles_CO_boxplot_lessnoise?.pdf"))
boxplot(lapply(ssID, function(x) sapply(1:100, function(ii) sum(rbinom(length(x), 1, prob = predict1[x]))/sum(dat$width[x])*1e6)), main = "GLM: CO~band+REC8 ChIP*(annotation+width)", xlab = "REC8 ChIP hexiles", ylab = "Number of crossovers per Mb", names = c("1", "2", "3", "4", "5", "6"))
tco <- sapply(ssID, function(x) sum(dat$co[x])/sum(dat$width[x])*1e6)
points(x = 1:length(ssID), y = tco, col = "red", pch = 19)
legend("top", legend = c("Observed", "Predicted"), text.col = c("red", "black"), bty = "n")
dev.off()

#######################################################
# plot predicted CO for regionsGR by nucleosomes hexile #
#######################################################

ssID <- split(1:nrow(dat), cut(dat$meannuc, breaks = c(0, sort(dat$meannuc)[round(0:6*(nrow(dat)/6))])))
pdf(file = paste0(outDir, "GLM_nuc_hexiles_CO_boxplot.pdf"))
tco <- sapply(ssID, function(x) sum(dat$co[x])/sum(dat$width[x])*1e6)
boxplot(lapply(ssID, function(x) sapply(1:100, function(ii) sum(rbinom(length(x), 1, prob = predict2[x]))/sum(dat$width[x])*1e6)), main = "GLM: CO~band+nucleosomes*(annotation+width)", xlab = "Nucleosomes hexiles", ylab = "Number of crossovers per Mb", names = c("1", "2", "3", "4", "5", "6"), ylim = c(min(tco), 40))
points(x = 1:length(ssID), y = tco, col = "red", pch = 19)
legend("top", legend = c("Observed", "Predicted"), text.col = c("red", "black"), bty = "n")
dev.off()

ssID <- split(1:nrow(dat), cut(dat$meannuc, breaks = sort(dat$meannuc)[round(0:7*(nrow(dat)/7))]))
pdf(file = paste0(outDir, "GLM_nuc_hexiles_CO_boxplot_lessnoise?.pdf"))
tco <- sapply(ssID, function(x) sum(dat$co[x])/sum(dat$width[x])*1e6)
boxplot(lapply(ssID, function(x) sapply(1:100, function(ii) sum(rbinom(length(x), 1, prob = predict2[x]))/sum(dat$width[x])*1e6)), main = "GLM: CO~band+nucleosomes*(annotation+width)", xlab = "Nucleosomes hexiles", ylab = "Number of crossovers per Mb", names = c("1", "2", "3", "4", "5", "6"), ylim = c(min(tco), max(tco)))
points(x = 1:length(ssID), y = tco, col = "red", pch = 19)
legend("top", legend = c("Observed", "Predicted"), text.col = c("red", "black"), bty = "n")
dev.off()

##########################################################
# plot predicted CO for regionsGR by SPO11-1-oligos hexile #
##########################################################

ssID <- split(1:nrow(dat), cut(dat$meanspo, breaks = c(0, sort(dat$meanspo)[round(0:6*(nrow(dat)/6))])))
pdf(file = paste0(outDir, "GLM_SPO11-1-oligos_hexiles_CO_boxplot.pdf"))
tco <- sapply(ssID, function(x) sum(dat$co[x])/sum(dat$width[x])*1e6)
boxplot(lapply(ssID, function(x) sapply(1:100, function(ii) sum(rbinom(length(x), 1, prob = predict3[x]))/sum(dat$width[x])*1e6)), main = "GLM: CO~band+SPO11-1-oligos*(annotation+width)", xlab = "SPO11-1-oligos hexiles", ylab = "Number of crossovers per Mb", names = c("1", "2", "3", "4", "5", "6"), ylim = c(min(tco), 53))
points(x = 1:length(ssID), y = tco, col = "red", pch = 19)
legend("top", legend = c("Observed", "Predicted"), text.col = c("red", "black"), bty = "n")
dev.off()

ssID <- split(1:nrow(dat), cut(dat$meanspo, breaks = sort(dat$meanspo)[round(0:7*(nrow(dat)/7))]))
pdf(file = paste0(outDir, "GLM_SPO11-1-oligos_hexiles_CO_boxplot_lessnoise?.pdf"))
tco <- sapply(ssID, function(x) sum(dat$co[x])/sum(dat$width[x])*1e6)
boxplot(lapply(ssID, function(x) sapply(1:100, function(ii) sum(rbinom(length(x), 1, prob = predict3[x]))/sum(dat$width[x])*1e6)), main = "GLM: CO~band+SPO11-1-oligos*(annotation+width)", xlab = "SPO11-1-oligos hexiles", ylab = "Number of crossovers per Mb", names = c("1", "2", "3", "4", "5", "6"))
points(x = 1:length(ssID), y = tco, col = "red", pch = 19)
legend("top", legend = c("Observed", "Predicted"), text.col = c("red", "black"), bty = "n")
dev.off()

#####################################################
# plot predicted CO for regionsGR by MSH4_ChIP hexile #
#####################################################

ssID <- split(1:nrow(dat), cut(dat$meanMSH4_ChIP, breaks = c(0, sort(dat$meanMSH4_ChIP)[round(0:6*(nrow(dat)/6))])))
pdf(file = paste0(outDir, "GLM_MSH4_ChIP_hexiles_CO_boxplot.pdf"))
boxplot(lapply(ssID, function(x) sapply(1:100, function(ii) sum(rbinom(length(x), 1, prob = predict4[x]))/sum(dat$width[x])*1e6)), main = "GLM: CO~band+MSH4 ChIP*(annotation+width)", xlab = "MSH4 ChIP hexiles", ylab = "Number of crossovers per Mb", names = c("1", "2", "3", "4", "5", "6"))
tco <- sapply(ssID, function(x) sum(dat$co[x])/sum(dat$width[x])*1e6)
points(x = 1:length(ssID), y = tco, col = "red", pch = 19)
legend("top", legend = c("Observed", "Predicted"), text.col = c("red", "black"), bty = "n")
dev.off()

ssID <- split(1:nrow(dat), cut(dat$meanMSH4_ChIP, breaks = sort(dat$meanMSH4_ChIP)[round(0:7*(nrow(dat)/7))]))
pdf(file = paste0(outDir, "GLM_MSH4_ChIP_hexiles_CO_boxplot_lessnoise?.pdf"))
boxplot(lapply(ssID, function(x) sapply(1:100, function(ii) sum(rbinom(length(x), 1, prob = predict4[x]))/sum(dat$width[x])*1e6)), main = "GLM: CO~band+MSH4 ChIP*(annotation+width)", xlab = "MSH4 ChIP hexiles", ylab = "Number of crossovers per Mb", names = c("1", "2", "3", "4", "5", "6"))
tco <- sapply(ssID, function(x) sum(dat$co[x])/sum(dat$width[x])*1e6)
points(x = 1:length(ssID), y = tco, col = "red", pch = 19)
legend("top", legend = c("Observed", "Predicted"), text.col = c("red", "black"), bty = "n")
dev.off()

##########################################
# plot predicted CO for chromosome bands #
##########################################

pdf(file = paste0(outDir, "GLM_REC8_ChIP_chrbands_CO_boxplot.pdf"), height = 6, width = 9)
boxplot(lapply(split(1:nrow(dat), dat$band), function(x) sapply(1:100, function(ii) sum(rbinom(length(x), 1, prob = predict1[x]))/sum(dat$width[x])*1e6)), main = "GLM: CO~band+REC8 ChIP*(annotation+width)", xlab = "2-Mb chromosome bands", ylab = "Number of crossovers per Mb")
ssb <- split(1:nrow(dat), dat$band)
points(x = unique(dat$band), y = sapply(ssb, function(x) sum(dat$co[x])/sum(dat$width[x])*1e6), col = "red", pch = 19)
abline(v = c(8, 17, 31, 37, 50), lty = 2, lwd = 2.5)
legend("topleft", legend = c("Observed", "Predicted"), text.col = c("red", "black"), bty = "n")
dev.off()

pdf(file = paste0(outDir, "GLM_REC8_MSH4_nuc_SPO11_chrbands_CO_boxplot.pdf"), height = 6, width = 9)
boxplot(lapply(split(1:nrow(dat), dat$band), function(x) sapply(1:100, function(ii) sum(rbinom(length(x), 1, prob = predict5[x]))/sum(dat$width[x])*1e6)), main = "GLM: CO~band+(REC8 ChIP+MSH4 ChIP+nucleosomes+SPO11-1-oligos)*(annotation+width)", xlab = "2-Mb chromosome bands", ylab = "Number of crossovers per Mb")
ssb <- split(1:nrow(dat), dat$band)
points(x = unique(dat$band), y = sapply(ssb, function(x) sum(dat$co[x])/sum(dat$width[x])*1e6), col = "red", pch = 19)
abline(v = c(8, 17, 31, 37, 50), lty = 2, lwd = 2.5)
legend("topleft", legend = c("Observed", "Predicted"), text.col = c("red", "black"), bty = "n")
dev.off()

#########################################
# plot predicted CO for each SNP region #
#########################################

#pdf(file = paste0(outDir, "GLM_REC8_ChIP_snpreg_CO.pdf"))
#xy <- cbind(x = cumsum(width(regionsGR)), y = (predict1))
#xy <- cbind(x = cumsum(width(regionsGR)), y = (predict1))
#xy <- cbind(x = cumsum(width(regionsGR)), y = (predict1))
#xy <- cbind(x = cumsum(width(regionsGR)), y = (predict1))
#xy <- cbind(x = cumsum(width(regionsGR)), y = (predict1))
#xy <- cbind(x = cumsum(width(regionsGR)), y = (predict1))
#plot(xy[xy[,2]>0.01,], pch = ".", main = "Likelihood of CO at each SNP window");
#abline(v = cumsum(chrlens), col = "red", lty = 3)
#abline(v = c(0, cumsum(chrlens)[-5])+start(centGR), col = "blue", lty = 3)
#abline(v = c(0, cumsum(chrlens)[-5])+end(centGR), col = "blue", lty = 3)
#dev.off()

sessionInfo()

