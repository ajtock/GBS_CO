##################################################################################
# Generate normalised coverage files for unique-both data                        #
##################################################################################

gbs_dir <- "/projects/ajt200/GBS_Chris/wtcu_co/" 
plot_dir <- "/projects/ajt200/GBS_Chris/wtcu_co/plots/"
bed_dir <- "/home/meiosis/ajt200/BS_Seq/Stroud_2013/wig/bed/"
library(Biostrings)
library(GenomicAlignments)
library(segmentSeq)
library(parallel)

######################
# Define chromosomes #
######################

chr_ends <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
cen_start <- c(13999493, 2930298, 12782751, 3400400, 10934103)
cen_end <- c(15983278, 4804568, 14750881, 4169512, 13194594)
pericen_start <- c(13999493, 2930298, 12782751, 3400400, 10934103)
pericen_end <- c(15983278, 4804568, 14750881, 4169512, 13194594)
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")

################################
# make cumulative genomes      #
################################

sumchr <- cumsum(c(0, chr_ends))
print(sumchr)
sumchr_tot <- sumchr[length(sumchr)]
print(sumchr_tot)

centromeres <- sapply(seq_along(centromeres), function(x) {
  centromeres[x] + sumchr[x]
})
print(centromeres)
pericen_start <- sapply(seq_along(pericen_start), function(x) {
  pericen_start[x] + sumchr[x]
})
print(pericen_start)
pericen_end <- sapply(seq_along(pericen_end), function(x) {
  pericen_end[x] + sumchr[x]
})
print(pericen_end)

####################################################################
# Make cumulative COs using GRanges object containing CO intervals #
# Use interval midpoints                                           #
####################################################################

load(file = paste0(gbs_dir, "cos.gr.coords.RData"))
wtco <- cos.gr.coords
GSM980986_WT_rep2_CG <- read.table(file = paste0(bed_dir, "GSM980986_WT_rep2_CG.wig.bed.gr.tab.bed"))
GSM980986_WT_rep2_CHG <- read.table(file = paste0(bed_dir, "GSM980986_WT_rep2_CHG.wig.bed.gr.tab.bed"))
GSM980986_WT_rep2_CHH <- read.table(file = paste0(bed_dir, "GSM980986_WT_rep2_CHH.wig.bed.gr.tab.bed"))

cum_co_meth_dat <- NULL
for(i in 1:5) {
  # define windows as GRanges object
  windows <- seq(1, chr_ends[i], by = 10000)
  windows <- c(windows, chr_ends[i])
  cum_windows <- windows + sumchr[i]
  windows_iranges <- IRanges(start = windows, width = 10000)
  windows_granges <- GRanges(seqnames = chrs[i], strand = "+", ranges = windows_iranges)
  
  # define and count wt COs within windows
  wtco_chr <- wtco[seqnames(wtco) == chrs[i]]
  win_wtco <- countOverlaps(windows_granges, wtco_chr)
  norm_win_wtco <- win_wtco / length(wtco)
  
  # calculate mean methylation levels (all contexts) within windows
  wtCG <- GSM980986_WT_rep2_CG[which(GSM980986_WT_rep2_CG$V1 == paste0("chr", i)),]
  wtCG_ir_coords <- IRanges(start = wtCG$V2, width = 1)
  wtCG_gr_coords <- GRanges(seqnames = chrs[i], strand = "+", ranges = wtCG_ir_coords)
  overlapswtCG <- getOverlaps(windows_granges, wtCG_gr_coords, whichOverlaps = TRUE)
  wtCG_win_vals <- sapply(overlapswtCG, function(x) mean(as.numeric(wtCG$V4[x])))

  wtCHG <- GSM980986_WT_rep2_CHG[which(GSM980986_WT_rep2_CHG$V1 == paste0("chr", i)),]
  wtCHG_ir_coords <- IRanges(start = wtCHG$V2, width = 1)
  wtCHG_gr_coords <- GRanges(seqnames = chrs[i], strand = "+", ranges = wtCHG_ir_coords)
  overlapswtCHG <- getOverlaps(windows_granges, wtCHG_gr_coords, whichOverlaps = TRUE)
  wtCHG_win_vals <- sapply(overlapswtCHG, function(x) mean(as.numeric(wtCHG$V4[x])))

  wtCHH <- GSM980986_WT_rep2_CHH[which(GSM980986_WT_rep2_CHH$V1 == paste0("chr", i)),]
  wtCHH_ir_coords <- IRanges(start = wtCHH$V2, width = 1)
  wtCHH_gr_coords <- GRanges(seqnames = chrs[i], strand = "+", ranges = wtCHH_ir_coords)
  overlapswtCHH <- getOverlaps(windows_granges, wtCHH_gr_coords, whichOverlaps = TRUE)
  wtCHH_win_vals <- sapply(overlapswtCHH, function(x) mean(as.numeric(wtCHH$V4[x])))

  dat <- cbind(cum_windows, norm_win_wtco, wtCG_win_vals, wtCHG_win_vals, wtCHH_win_vals)
  cum_co_meth_dat <- rbind(cum_co_meth_dat, dat)
}
write.table(cum_co_meth_dat, file = paste0(gbs_dir, "wtcuco_wtmeth_genome_10kb.txt"))

test <- seq(1, 1000, by = 1)
j = 200
ma <- rep(1, test[j])/test[j]
filt_wtco <- stats::filter(cum_co_meth_dat[,2], ma)
which_na <- which(is.na(filt_wtco) == TRUE)
left_na <- which_na[which(which_na < 100)]
left_val <- filt_wtco[left_na[length(left_na)]+1]
filt_wtco[left_na] <- left_val
right_na <- which_na[which(which_na > 100)]
right_val <- filt_wtco[right_na[1]-1]
filt_wtco[right_na] <- right_val
filt_wtco_n <- filt_wtco[!is.na(filt_wtco)]
filt <- cbind(cum_co_meth_dat[,1], filt_wtco)
write.table(filt, file = paste0(gbs_dir, "filt_wtcuco_genome_10kb.txt"))

filt_wtCG <- stats::filter(cum_co_meth_dat[,3], ma); filt_wtCG_n <- filt_wtCG[!is.na(filt_wtCG)]
filt_wtCHG <- stats::filter(cum_co_meth_dat[,4], ma); filt_wtCHG_n <- filt_wtCHG[!is.na(filt_wtCHG)]
filt_wtCHH <- stats::filter(cum_co_meth_dat[,5], ma); filt_wtCHH_n <- filt_wtCHH[!is.na(filt_wtCHH)]

filt_wtMeth <- cbind(cum_co_meth_dat[,1], filt_wtCG, filt_wtCHG, filt_wtCHH)
write.table(filt_wtMeth, file = paste0(gbs_dir, "filt_wtMeth_genome_10kb.txt"))

pdf(paste0(plot_dir, "wtco_wtmeth_genome_10kb_v260417.pdf"), height = 4, width = 10)
par(mfrow = c(1, 1))
par(mar = c(5.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
ymin <- min(filt_wtco_n)
ymax <- max(filt_wtco_n)
ymin2 <- min(filt_wtCG_n, filt_wtCHG_n, filt_wtCHH_n)
ymax2 <- max(filt_wtCG_n, filt_wtCHG_n, filt_wtCHH_n)
plot(filt_wtMeth[,1], filt_wtMeth[,2], type = "l", col = "red",
     ylim = c(ymin2, ymax2), ylab = "",
     ann = F, xaxt = "n", yaxt = "n")
#lines(filt_wtMeth[,1], filt_wtMeth[,3], type = "l", col = "purple")
#lines(filt_wtMeth[,1], filt_wtMeth[,4], type = "l", col = "blue")
par(new = TRUE)
plot(filt_wtMeth[,1], filt_wtMeth[,3], type = "l", col = "purple",
     ylim = c(ymin2, ymax2), ylab = "",
     ann = F, xaxt = "n", yaxt = "n")
par(new = TRUE)
plot(filt_wtMeth[,1], filt_wtMeth[,4], type = "l", col = "blue",
     ylim = c(ymin2, ymax2), ylab = "",
     ann = F, xaxt = "n", yaxt = "n")
axis(side = 4)
mtext(side = 4, line = 2, cex = 1, text = "Methylation proportion")
par(new = TRUE)
plot(filt[,1], filt[,2], type = "l", col = "black",
     ylim = c(ymin, ymax),
     main = "10-kb windows",
     xlab = "Coordinates (bp)",
     ylab = "Normalised CO frequency")
#abline(h = mean(filt_wtco_n), lty = 2, lwd = 1, col = "red")
abline(v = sumchr, lty = 1, lwd = 1)
abline(v = centromeres, lty = 2, lwd = 1, col = "red")
abline(v = pericen_start, lty = 3, lwd = 2)
abline(v = pericen_end, lty = 3, lwd = 2)
legend("topright",
       legend = c("wt COs", "wt CG", "wt CHG", "wt CHH"),
       col = c("black", "red", "purple", "blue"),
       text.col = c("black", "red", "purple", "blue"),
       ncol = 1, cex = 0.8, lwd = 1.2, bty = "n")
dev.off()


## alternative approach to CO counting within windows - results are comparable
#cum_co_dat2 <- NULL
#for(i in 1:5) {
#  windows <- seq(1, chr_ends[i], by = 10000)
#  windows <- c(windows, chr_ends[i])
#  wtco_chr <- wtco[seqnames(wtco) == chrs[i]]
#  wtco_chr_ranges <- ranges(wtco_chr)
#  mid_interval <- width(wtco_chr_ranges) / 2
#  co_coords <- as.numeric(start(wtco_chr_ranges) + mid_interval)
#  norm_win_wtco2 <- NULL
#  for(k in 1:length(windows)-1) {
#    print(k)
#    which_co <- length(which(co_coords >= windows[k] & co_coords < windows[k+1]))
#    norm_which_co <- which_co / length(wtco)
#    norm_win_wtco2 <- c(norm_win_wtco2, norm_which_co)
#  }
#  cum_windows <- windows + sumchr[i]
#  dat <- cbind(cum_windows, norm_win_wtco2)
#  cum_co_dat2 <- rbind(cum_co_dat2, dat)
#}
#write.table(cum_co_dat2, file = paste0(gbs_dir, "wtcu_co_genome_10kb_2.txt"))
#
#test <- seq(1, 1000, by = 1)
#j = 200
#ma <- rep(1, test[j])/test[j]
#filt_wtco2 <- stats::filter(cum_co_dat2[,2], ma)
#which_na <- which(is.na(filt_wtco2) == TRUE)
#left_na <- which_na[which(which_na < 100)]
#left_val <- filt_wtco2[left_na[length(left_na)]+1]
#filt_wtco2[left_na] <- left_val
#right_na <- which_na[which(which_na > 100)]
#right_val <- filt_wtco2[right_na[1]-1]
#filt_wtco2[right_na] <- right_val
#filt_wtco2_n <- filt_wtco2[!is.na(filt_wtco2)]
#filt2 <- cbind(cum_co_dat2[,1], filt_wtco2)
#write.table(filt2, file = paste0(gbs_dir, "filt_wtcu_co_genome_10kb_2.txt"))




##################################################################################
# Generate mean methylation values for windows per BS library and context  #
##################################################################################

#BS_id <- c("GSM980986_WT_rep2_CG", "GSM980986_WT_rep2_CHG", "GSM980986_WT_rep2_CHH")

#GSM980986_WT_rep2_CG <- read.table(file = paste0(bed_dir, "GSM980986_WT_rep2_CG.wig.bed.gr.tab.bed"))
#GSM980986_WT_rep2_CHG <- read.table(file = paste0(bed_dir, "GSM980986_WT_rep2_CHG.wig.bed.gr.tab.bed"))
#GSM980986_WT_rep2_CHH <- read.table(file = paste0(bed_dir, "GSM980986_WT_rep2_CHH.wig.bed.gr.tab.bed"))

#ids <- c("GSM980986", "GSM981031")

#plotMethBed(ids, splitContexts = TRUE)

#BSGenomePlotData <- function() {
#  read.table(file = paste0(bed_dir, "GSM980986_WT_rep2_CG.wig.bed.gr.tab.bed"

sessionInfo()

