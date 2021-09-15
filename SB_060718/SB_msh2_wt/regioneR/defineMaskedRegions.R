#!/applications/R/R-3.3.2/bin/Rscript

# Create and save regions to mask in overlap analyses
library(regioneR)

outDir <- "/projects/ajt200/GBS_CO/SB_060718/SB_msh2_wt/regioneR/"

# Col/Ler mask
maskGR <- toGRanges(data.frame(chr   = c("Chr3",
                                         rep("Chr4", times = 3)),
                               start = c(3000000,
                                         131000, 170000, 7800500),
                               end   = c(9000000,
                                         140000, 179000, 8030000)))
save(maskGR,
     file = paste0(outDir, "SB_ColLer_F2_2018_maskGR.RData"))
maskGR <- NULL

# Col/Cvi mask
maskGR <- toGRanges(data.frame(chr   = c(rep("Chr1", times = 3),
                                         rep("Chr2", times = 5),
                                         "Chr3",
                                         rep("Chr4", times = 2)),
                               start = c(3460000, 4250000, 11303000,
                                         446500, 527100, 6500000, 16190000, 18870000,
                                         4000000,
                                         3125000, 4171000),
                               end   = c(3750000, 4650000, 11330000,
                                         457300, 538500, 10300000, 16205000, 18880000,
                                         9500000,
                                         3400000, 4802000)))
save(maskGR,
     file = paste0(outDir, "SB_ColCvi_F2_2018_maskGR.RData"))
maskGR <- NULL
