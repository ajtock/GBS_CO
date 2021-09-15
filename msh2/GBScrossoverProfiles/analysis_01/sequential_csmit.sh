#!/bin/bash

csmit -m 50G -c 10 "./crossover_Profiles_commandArgs.R 5000 5kb 20 /projects/ajt200/BAM_masters/nucleosomes/WT/coverage/nakedDNA_untrimmed_input/log2ChIPinput/log2wtNucNakedDNAuntrimmed_norm_allchrs_coverage_coord_tab.bed MNase" & sleep 10;
csmit -m 50G -c 10 "./crossover_Profiles_commandArgs.R 5000 5kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/log2wtSPO11oligoRPI1NakedDNA_norm_allchrs_coverage_coord_tab.bed SPO11_1_oligos_RPI1" & sleep 10;
wait

