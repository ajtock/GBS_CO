#######################################################################
# Define crossover intervals using 96 GBS-derived haplotype files     #
#######################################################################

library(dplyr)
sessionInfo()

in.dir1 <- "/projects/ajt200/GBS_Chris/HEI10_cmt3/"
in.dir2 <- "/projects/ajt200/GBS_Chris/recq4ab_cmt3/"

HEI10_cmt3 <- paste0(in.dir1, "n96.CL.HEI10.cmt3.lane1.",
                     c(1:96), ".bchqsnvmask.smooth.co.txt")
recq4ab_cmt3 <- paste0(in.dir2, "n96.CL.recq.cmt3.lane1.",
                       c(1:96), ".bchqsnvmask.smooth.co.txt")

HEI10_cmt3_haplo <- lapply(seq_along(HEI10_cmt3), function(x) {
  read.table(HEI10_cmt3[[x]])
})

pop_HEI10_cmt3_co_intervals <- NULL
for(i in 1:96) {
  print(i)
  ind_HEI10_cmt3_co_all <- NULL
  for(j in 1:length(HEI10_cmt3_haplo[[i]]$V3)+1) {
    ind_HEI10_cmt3_co <- cbind(HEI10_cmt3_haplo[[i]]$V1[j],
                               HEI10_cmt3_haplo[[i]]$V2[j],
                               HEI10_cmt3_haplo[[i]]$V4[j-1],
                               HEI10_cmt3_haplo[[i]]$V3[j])
    ind_HEI10_cmt3_co_all <- rbind(ind_HEI10_cmt3_co_all, ind_HEI10_cmt3_co)
    ind_HEI10_cmt3_co_all <- ind_HEI10_cmt3_co_all[which(ind_HEI10_cmt3_co_all[,3] < ind_HEI10_cmt3_co_all[,4]),]
  }
  pop_HEI10_cmt3_co_intervals <- rbind(pop_HEI10_cmt3_co_intervals, ind_HEI10_cmt3_co_all)
  pop_HEI10_cmt3_co_intervals <- tbl_df(pop_HEI10_cmt3_co_intervals)
}
co_interval_size <- pop_HEI10_cmt3_co_intervals$V4-pop_HEI10_cmt3_co_intervals$V3
pop_HEI10_cmt3_co_intervals <- cbind(pop_HEI10_cmt3_co_intervals, co_interval_size)
pop_HEI10_cmt3_co_intervals <- tbl_df(pop_HEI10_cmt3_co_intervals)
pop_HEI10_cmt3_co_intervals_sort <- arrange(pop_HEI10_cmt3_co_intervals, V2, V3)
write.table(pop_HEI10_cmt3_co_intervals_sort, file=paste0(in.dir1, "HEI10_cmt3_co_intervals.txt"))
HEI10_cmt3_gff <- cbind(pop_HEI10_cmt3_co_intervals_sort[,2],
                        rep(".", length(pop_HEI10_cmt3_co_intervals_sort[,2])),
                        rep("CO_interval", length(pop_HEI10_cmt3_co_intervals_sort[,2])),
                        pop_HEI10_cmt3_co_intervals_sort[,3],
                        pop_HEI10_cmt3_co_intervals_sort[,4],
                        pop_HEI10_cmt3_co_intervals_sort[,5],
                        rep("+", length(pop_HEI10_cmt3_co_intervals_sort[,2])),
                        rep(".", length(pop_HEI10_cmt3_co_intervals_sort[,2])),
                        rep(".", length(pop_HEI10_cmt3_co_intervals_sort[,2])))
colnames(HEI10_cmt3_gff) <- c("chr", "source", "feature", "start", "end", "size", "strand", "frame", "ann")
write.table(HEI10_cmt3_gff, file=paste0(in.dir1, "HEI10_cmt3_co_intervals.gff"), row.names=F, sep="\t")


recq4ab_cmt3_haplo <- lapply(seq_along(recq4ab_cmt3), function(x) {
  read.table(recq4ab_cmt3[[x]])
})

pop_recq4ab_cmt3_co_intervals <- NULL
for(i in 1:96) {
  print(i)
  ind_recq4ab_cmt3_co_all <- NULL
  for(j in 1:length(recq4ab_cmt3_haplo[[i]]$V3)+1) {
    ind_recq4ab_cmt3_co <- cbind(recq4ab_cmt3_haplo[[i]]$V1[j],
                               recq4ab_cmt3_haplo[[i]]$V2[j],
                               recq4ab_cmt3_haplo[[i]]$V4[j-1],
                               recq4ab_cmt3_haplo[[i]]$V3[j])
    ind_recq4ab_cmt3_co_all <- rbind(ind_recq4ab_cmt3_co_all, ind_recq4ab_cmt3_co)
    ind_recq4ab_cmt3_co_all <- ind_recq4ab_cmt3_co_all[which(ind_recq4ab_cmt3_co_all[,3] < ind_recq4ab_cmt3_co_all[,4]),]
  }
  pop_recq4ab_cmt3_co_intervals <- rbind(pop_recq4ab_cmt3_co_intervals, ind_recq4ab_cmt3_co_all)
  pop_recq4ab_cmt3_co_intervals <- tbl_df(pop_recq4ab_cmt3_co_intervals)
}
co_interval_size <- pop_recq4ab_cmt3_co_intervals$V4-pop_recq4ab_cmt3_co_intervals$V3
pop_recq4ab_cmt3_co_intervals <- cbind(pop_recq4ab_cmt3_co_intervals, co_interval_size)
pop_recq4ab_cmt3_co_intervals <- tbl_df(pop_recq4ab_cmt3_co_intervals)
pop_recq4ab_cmt3_co_intervals_sort <- arrange(pop_recq4ab_cmt3_co_intervals, V2, V3)
write.table(pop_recq4ab_cmt3_co_intervals_sort, file=paste0(in.dir2, "recq4ab_cmt3_co_intervals.txt"))
recq4ab_cmt3_gff <- cbind(pop_recq4ab_cmt3_co_intervals_sort[,2],
                        rep(".", length(pop_recq4ab_cmt3_co_intervals_sort[,2])),
                        rep("CO_interval", length(pop_recq4ab_cmt3_co_intervals_sort[,2])),
                        pop_recq4ab_cmt3_co_intervals_sort[,3],
                        pop_recq4ab_cmt3_co_intervals_sort[,4],
                        pop_recq4ab_cmt3_co_intervals_sort[,5],
                        rep("+", length(pop_recq4ab_cmt3_co_intervals_sort[,2])),
                        rep(".", length(pop_recq4ab_cmt3_co_intervals_sort[,2])),
                        rep(".", length(pop_recq4ab_cmt3_co_intervals_sort[,2])))
colnames(recq4ab_cmt3_gff) <- c("chr", "source", "feature", "start", "end", "size", "strand", "frame", "ann")
write.table(recq4ab_cmt3_gff, file=paste0(in.dir2, "recq4ab_cmt3_co_intervals.gff"), row.names=F, sep="\t")


