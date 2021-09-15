######################################################
# Count number of COs within pericentromeric regions #
######################################################

library(GenomicRanges)
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
# pericentromeric regions are as defined in Supplemental Table S26 of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

pericenGR <- do.call("c", lapply(1:5, function(chr) {
  GRanges(seqnames = chrs[chr], IRanges(start = pericenStart[chr], end = pericenEnd[chr]))
}))
print(pericenGR)

coDir <- "/home/meiosis/cal66/GBS_mega_recombination/"
cos.wtcu.gr.coords <- load(file = paste0(coDir, "cos.wtcu.gr.coords"))
cos.wtcu.gr.coords <- sort(cos.gr.coords)
cos.cmt3.gr.coords <- load(file = paste0(coDir, "cos.cmt3.gr.coords"))
cos.cmt3.gr.coords <- sort(cos.gr.coords)
cos.hei10.gr.coords <- load(file = paste0(coDir, "cos.hei10.gr.coords"))
cos.hei10.gr.coords <- sort(cos.gr.coords)
cos.hei10_cmt3.gr.coords <- load(file = paste0(coDir, "cos.hei10_cmt3.gr.coords"))
cos.hei10_cmt3.gr.coords <- sort(cos.gr.coords)
cos.recq4ab_cmt3.gr.coords <- load(file = paste0(coDir, "cos.recq4ab_cmt3.gr.coords"))
cos.recq4ab_cmt3.gr.coords <- sort(cos.gr.coords)
cos.recq4ab.gr.coords <- load(file = paste0(coDir, "cos.recq4ab.gr.coords"))
cos.recq4ab.gr.coords <- sort(cos.gr.coords)

coGRList <- GRangesList(cos.wtcu.gr.coords, cos.cmt3.gr.coords, cos.hei10.gr.coords, cos.hei10_cmt3.gr.coords, cos.recq4ab_cmt3.gr.coords, cos.recq4ab.gr.coords)
coNames <- c("cos.wtcu.gr.coords", "cos.cmt3.gr.coords", "cos.hei10.gr.coords", "cos.hei10_cmt3.gr.coords", "cos.recq4ab_cmt3.gr.coords", "cos.recq4ab.gr.coords")

coCountCen <- lapply(seq_along(coGRList), function(x) {
                countOverlaps(pericenGR, coGRList[[x]])
              })
names(coCountCen) <- coNames
print(coCountCen)
#$cos.wtcu.gr.coords
#[1] 54 71 47 48 56
#
#$cos.cmt3.gr.coords
#[1] 92 81 92 82 85
#
#$cos.hei10.gr.coords
#[1] 125 105  84  36  95 # alternative approach below counts 94 COs within the pericentromeric region of Chr05, 
#                        # omitting one CO interval whose start coordinate is < pericenStart[5]
#$cos.hei10_cmt3.gr.coords
#[1] 128 141  96  48 108
#
#$cos.recq4ab_cmt3.gr.coords
#[1] 127 104  75  66 112

## A more convoluted sanity check to confirm above CO counts 
## and return CO intervals within pericentromeric regions as GRanges objects
coCountCenAlt <- lapply(seq_along(coGRList), function(x) {
                   lapply(1:5, function(y) {
                     coChr <- coGRList[[x]][seqnames(coGRList[[x]]) == chrs[y],]
                     coChrPericen <- coChr[start(coChr) >= pericenStart[y] &
                                        end(coChr) <= pericenEnd[y],]
                   })
                 })
names(coCountCenAlt) <- coNames
print(coCountCenAlt)


sessionInfo()
