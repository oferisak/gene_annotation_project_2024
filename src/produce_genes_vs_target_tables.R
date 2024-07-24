setwd('/media/SSD/Bioinformatics/Projects/gene_annotation_project_2024/')
source('./lib/initialization.R')

# target options
target_files<-list('twist_v2'='./data/target_files/twist_hg19_exome_comp_spikein_v2.0.2_targets_sorted.re_annotated_0.pad50.bed',
                   'cegat_exome_v5'='./data/target_files/S000054_cexome5_hg19_snvs_targets.pad50.bed')

# genes files
gene_files<-list('refseq'='./data/accessory_files/GCF_000001405.25_hg19.p13.refseq.bed')

# selected target file
target_file_name<-'cegat_exome_v5'
gene_file_name<-'refseq'

gene_file_cols<-c('chr','start','end','info','strand','exon_num','gene','transcript','exon_size','type')

for (target_file_name in names(target_files)){
  message(glue('checking the coverage of {gene_file_name} genes in target file {target_file_name}'))
  coverage_command<-glue('bedtools coverage -a {gene_files[[gene_file_name]]} -b {target_files[[target_file_name]]}')
  coverage_output<-read.table(text=system(coverage_command,intern = T))
  colnames(coverage_output)<-c(gene_file_cols,
                                  'num_o_overlapping_regions','overlap_size','total_size','overlap_percent')
  
  incomplete_coverage<-coverage_output%>%filter(overlap_percent!=1)
  
  # write the full coverage output
  write.table(coverage_output,file=glue('./output/gene_vs_target/{target_file_name}.{gene_file_name}.all_regions.tsv'),row.names = F,quote = F,sep = '\t')
  
  # write the disclaimer coverage output
  write.table(incomplete_coverage,
              file=glue('./output/gene_vs_target/{target_file_name}.{gene_file_name}.incomplete_coverage.tsv'),row.names = F,quote = F,sep = '\t')
  
  # per gene output
  incomplete_coverage_per_gene<-incomplete_coverage%>%
    group_by(gene,transcript)%>%
    summarize(exons_incomplete_coverage=paste0(glue('{exon_num}:{total_size-overlap_size} ({100*round(1-overlap_percent,3)}%)'),collapse='; '))
}
