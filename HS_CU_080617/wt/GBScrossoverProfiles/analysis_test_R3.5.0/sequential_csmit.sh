#!/bin/bash

csmit -m 30G -c 1 "/applications/R/R-3.5.0/bin/Rscript ./crossover_Profiles_commandArgs.R 5000 5kb 50 /home/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/WT/coverage/log2ChIPinput/REC8_MYC_Rep1_input/noZscore/log2_REC8_HA_Rep2_ChIP_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab_noZscore.bed REC8_HA_Rep2" & sleep 10;
csmit -m 30G -c 1 "/applications/R/R-3.5.0/bin/Rscript ./crossover_Profiles_commandArgs.R 5000 5kb 50 /projects/ajt200/BAM_masters/nucleosomes/WT/coverage/nakedDNA_untrimmed_input/log2ChIPinput/noZscore/log2wtNucNakedDNAuntrimmed_norm_allchrs_coverage_coord_tab_noZscore.bed MNase" & sleep 10;
csmit -m 30G -c 1 "/applications/R/R-3.5.0/bin/Rscript ./crossover_Profiles_commandArgs.R 5000 5kb 50 /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/noZscore/log2wtSPO11oligoRPI1NakedDNA_norm_allchrs_coverage_coord_tab_noZscore.bed SPO11_1_oligos_RPI1" & sleep 10;
