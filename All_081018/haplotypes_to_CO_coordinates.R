#!/applications/R/R-3.4.0/bin/Rscript

# Convert GBS-derived haplotype files into crossover interval files

# Usage:
# ./haplotypes_to_CO_coordinates.R "n96.wtF2" ".bchqsnvmask.smooth.co.txt" "CU_wt_ColLer_F2_2016_lane1"

args <- commandArgs(trailingOnly = T)
#prefix <- "n96.wtF2"
#suffix <- ".bchqsnvmask.smooth.co.txt"
#outName <- "CU_wt_ColLer_F2_2016_lane1"
prefix <- args[1]
suffix <- args[2]
outName <- args[3]

library(GenomicRanges)

inDir <- "/projects/ajt200/GBS_CO/All_081018/cos_files/"
outDir <- "/projects/ajt200/GBS_CO/All_081018/coordinates/"
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))

filePath <- system(paste0("ls ", inDir, prefix, "*", suffix),
                   intern = T)
popCOs <- NULL
for(h in 1:length(filePath)) {
  haplo <- read.table(filePath[h])
  allChrs <- NULL
  for(i in 1:5) {
    chr_haplo <- haplo[haplo[,2] == i,]
    if((length(chr_haplo[,1]) > 1) == T) {
      start <- chr_haplo[,4][-length(chr_haplo[,4])]
      end <- chr_haplo[,3][-1]
      midpoint <- start+( round( (end-start)/2 ) )
      chr <- rep(i, length(start))
      lib <- rep(chr_haplo[,1][1], length(start))
      width <- end-start+1
      geno1 <- as.character(chr_haplo[,5][-length(chr_haplo[,5])])
      geno2 <- as.character(chr_haplo[,5][-1])
      chrTable <- cbind(lib, chr, start, end, midpoint, width,
                        geno1, geno2)
      allChrs <- rbind(allChrs, chrTable)
      print(dim(allChrs))
    }
  }
  popCOs <- rbind(popCOs, allChrs)
}
popCOs <- as.data.frame(popCOs)
popCOs <- popCOs[order(popCOs$lib, decreasing = F),]
print(paste0(outName, " rows (COs) and columns:"))
print(dim(popCOs))
print(paste0(outName, " individuals:"))
print(length(unique(popCOs$lib)))
write.table(popCOs,
            paste0(outDir, outName, ".txt"),
            sep = "\t", col.names = T, row.names = F, quote = F)

# Convert dataframe into GRanges object
popCOs <- read.table(paste0(outDir, outName, ".txt"), header = T)
popCOsGR <- GRanges(seqnames = paste0("Chr", popCOs$chr),
                    ranges = IRanges(start = popCOs$start,
                                     end = popCOs$end),
                    strand = "*",
                    midpoint = popCOs$midpoint,
                    lib = popCOs$lib,
                    geno1 = as.character(popCOs$geno1),
                    geno2 = as.character(popCOs$geno2))
popCOsGR <- sortSeqlevels(popCOsGR)
popCOsGR <- sort(popCOsGR, ignore.strand = T)
save(popCOsGR,
     file = paste0(outDir, outName, "_GRanges.RData"))

