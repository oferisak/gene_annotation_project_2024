library(glue)
library(dplyr)
library(tidyr)
library(stringr)

source('./lib/clinvar_functions.R')
source('./lib/bedtools_functions.R')

# target options
target_files<-list('twist_v2'='./data/target_files/twist_hg19_exome_comp_spikein_v2.0.2_targets_sorted.re_annotated_0.pad50.bed',
                   'cegat_exome_v5'='./data/target_files/S000054_cexome5_hg19_snvs_targets.pad50.bed',
                   'idt_exome_v2'='./data/target_files/xgen-exome-research-panel-v2-targets-hg19.bed')


# genes files
gene_files<-list('refseq'='./data/accessory_files/GCF_000001405.25_hg19.p13.refseq.bed')