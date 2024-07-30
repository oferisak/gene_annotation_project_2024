setwd('/media/SSD/Bioinformatics/Projects/gene_annotation_project_2024/')
source('./lib/initialization.R')

giab_strats_folder<-glue('./data/giab_stratification_files/hg19')
giab_strats_files<-list.files(giab_strats_folder,recursive = T,pattern = '\\.bed',full.names = T)

# genes files
gene_files<-list('refseq'='./data/accessory_files/GCF_000001405.25_hg19.p13.refseq.bed')

# selected gene file
gene_file_name<-'refseq'

gene_file_cols<-c('chr','start','end','info','strand','exon_num','gene','transcript','size','type')

for (target_file_name in names(target_files)){
  message(glue('Producing GIAB strat files for {target_file_name}'))
  gene_giab_intersection<-readr::read_delim(gene_files[[gene_file_name]],col_names = gene_file_cols,delim='\t')
  base_cols<-colnames(gene_giab_intersection)
  for (giab_strat_file in giab_strats_files){
    strat_name<-stringr::str_replace(basename(giab_strat_file),'(_slop.+)?.bed.gz','')%>%stringr::str_replace('GRCh37_','')
    message(glue('intersecting with {strat_name}'))
    intersection_command<-glue('bedtools intersect -a {target_files[[target_file_name]]} -b {giab_strat_file} | bedtools coverage -a {gene_files[[gene_file_name]]} -b  -')
    intersection_output<-read.table(text=system(intersection_command,intern = T))
    colnames(intersection_output)<-c(gene_file_cols,
                                     'num_o_overlapping_regions','overlap_size','total_size','overlap_percent')
    gene_giab_intersection<-gene_giab_intersection%>%left_join(
      intersection_output%>%rename(!!sym(strat_name):=overlap_percent)%>%select(gene_file_cols,strat_name)
    )
  }
  strat_cols<-setdiff(colnames(gene_giab_intersection),base_cols)
  write.table(gene_giab_intersection,file=glue('./output/gene_giab_stratifications/{target_file_name}.gene_giab_stratifications.csv'),row.names = F,quote = F,sep = '\t')
  gene_giab_intersection_filtered<-gene_giab_intersection%>%
    filter(if_any(all_of(strat_cols), ~ . > 0))%>%
    mutate(non_zero_stratifications = apply(across(all_of(strat_cols)), 1, function(row) {
      non_zero_columns <- names(row)[row != 0]
      non_zero_values <- row[row != 0]
      paste(paste(non_zero_columns, "(", 100*round(non_zero_values,3), "%)", sep = ""), collapse = "; ")
    }))
  write.table(gene_giab_intersection_filtered,
              file=glue('./output/gene_giab_stratifications/{target_file_name}.gene_giab_stratifications.non_zero.csv'),
              row.names = F,quote = F,sep = '\t')
}


  
