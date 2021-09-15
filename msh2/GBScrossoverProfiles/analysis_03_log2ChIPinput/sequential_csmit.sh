#!/bin/bash

csmit -m 150G -c 11 "/applications/R/R-3.4.0/bin/Rscript ./crossover_Profiles_commandArgs.R 5000 5kb 20 /projects/ajt200/BAM_masters/nucleosomes/WT/coverage/nakedDNA_untrimmed_input/log2ChIPinput/noZscore/log2wtNucNakedDNAuntrimmed_norm_allchrs_coverage_coord_tab_noZscore.bed MNase" & sleep 10;
csmit -m 150G -c 11 "/applications/R/R-3.4.0/bin/Rscript ./crossover_Profiles_commandArgs.R 5000 5kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/noZscore/log2wtSPO11oligoRPI1NakedDNA_norm_allchrs_coverage_coord_tab_noZscore.bed SPO11_1_oligos_RPI1" & sleep 10;
wait
