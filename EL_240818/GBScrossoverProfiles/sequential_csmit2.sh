#!/bin/bash

csmit -m 20G -c 1 "/applications/R/R-3.3.2/bin/Rscript crossover_Profiles_commandArgs.R 5000 5kb 20 /projects/ajt200/BAM_masters/RNAseq_meiocyte_Feng/snakemake_RNAseq_STAR/mapped/both/pb/Col_RNAseq_meiocyte_Rep2_MappedOn_TAIR10_chr_all_both_sort_norm.perbase Col_RNAseq_meiocyte_Rep2" & sleep 5;
csmit -m 20G -c 1 "/applications/R/R-3.3.2/bin/Rscript crossover_Profiles_commandArgs.R 5000 5kb 20 /projects/ajt200/BAM_masters/RNAseq_meiocyte_Feng/snakemake_RNAseq_STAR/mapped/both/pb/taf4b_RNAseq_meiocyte_Rep2_MappedOn_TAIR10_chr_all_both_sort_norm.perbase taf4b_RNAseq_meiocyte_Rep2" & sleep 5;

csmit -m 20G -c 1 "/applications/R/R-3.3.2/bin/Rscript crossover_Profiles_commandArgs.R 5000 5kb 20 /projects/ajt200/BAM_masters/RNAseq_meiocyte_Feng/snakemake_RNAseq_STAR/mapped/unique/pb/Col_RNAseq_meiocyte_Rep2_MappedOn_TAIR10_chr_all_unique_sort_norm.perbase Col_RNAseq_meiocyte_Rep2_unique" & sleep 5;
csmit -m 20G -c 1 "/applications/R/R-3.3.2/bin/Rscript crossover_Profiles_commandArgs.R 5000 5kb 20 /projects/ajt200/BAM_masters/RNAseq_meiocyte_Feng/snakemake_RNAseq_STAR/mapped/unique/pb/taf4b_RNAseq_meiocyte_Rep2_MappedOn_TAIR10_chr_all_unique_sort_norm.perbase taf4b_RNAseq_meiocyte_Rep2_unique" & sleep 5;
wait
