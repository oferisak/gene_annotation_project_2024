setwd('/media/SSD/Bioinformatics/Projects/gene_annotation_project_2024/')
source('./lib/initialization.R')

# accessory files
alias_file<-'./data/accessory_files/Homo_sapiens.gene_info'

# definitions
# should the input clinvar file be filtered for P/LP variants?
is_filter_clinvar_vcf_for_plp<-TRUE
# should the clinvar file be converted into hg19 format?
is_convert_to_hg19<-TRUE

# clinvar options
clinvar_files<-list('clinvar_plp_20240708'='./data/clinvar_files/clinvar_20240708.vcf.gz')

# selected target and clinvar 
for (clinvar_file_name in names(clinvar_files)){
  clinvar_file<-clinvar_files[[clinvar_file_name]]
  
  if (is_filter_clinvar_vcf_for_plp){
    clinvar_file<-filter_clinvar_vcf_for_plp(clinvar_file)
  }
  
  if (is_convert_to_hg19){
    clinvar_file<-convert_vcf_to_hg19(clinvar_file,output_folder = './data/clinvar_files')
  }
  for (target_file_name in names(target_files)){
    message(glue('Generating ClinVar coverage files for {target_file_name} and {clinvar_file_name}'))
    # grab corresponding target and clinvar files
    target_file<-target_files[[target_file_name]]
    # grab gene aliases for output file
    gene_alias<-readr::read_delim(alias_file)%>%as.data.frame()
    # remove ambiguous genes
    gene_alias<-gene_alias%>%filter(!duplicated(Symbol))
    # now intersect the VCF with the target BED
    coverage_command<-glue('bedtools coverage -a {clinvar_file} -b {target_files[[target_file_name]]}')
    coverage_output<-read.table(text=system(coverage_command,intern = T),comment.char = '')
    colnames(coverage_output)<-c('chr','pos','id','ref','alt','X1','X2','info','num_o_overlapping_regions','num_o_covered_bases','variant_size','variant_fraction_covered')
    
    # parse INFO file into columns
    coverage_output <- coverage_output %>%
      mutate(
        allele_id = stringr::str_extract(info, 'ALLELEID=(\\d+)', group = 1),
        clnsig = stringr::str_extract(info, 'CLNSIG=([^;]+)', group = 1),
        vartype = stringr::str_extract(info, 'MC=([^;]+)', group = 1),
        vartype = ifelse(is.na(vartype),str_extract(info, 'CLNVC=([^;]+)', group = 1),vartype),
        vartype = ifelse(grepl(',', vartype),
                         stringr::str_split(vartype, ',') %>% unlist() %>% head(1),
                         vartype),
        gene=stringr::str_extract(info, 'GENEINFO=([^:]+)', group = 1),
        var=glue('{chr}:{pos}{ref}>{alt}'),
        in_target=ifelse(variant_fraction_covered>0,'in_target','not_in_target')
      )
    if (is_filter_clinvar_vcf_for_plp){
      coverage_output<-coverage_output%>%filter(grepl('pathogenic',clnsig,ignore.case=T))
    }
    # clinvar coverage table
    clinvar_coverage_table<-coverage_output%>%
      select(gene,id,chr,pos,ref,alt,var,clnsig,vartype,in_target)%>%
      left_join(gene_alias%>%select(gene=Symbol,synonyms=Synonyms))
    
    # Full table
    full_clinvar_coverage_file_name<-glue('{target_file_name}.{clinvar_file_name}.clinvar_coverage.all_variants.tsv')
    write.table(clinvar_coverage_table%>%select(-synonyms),file=glue('./output/clinvar_coverage/{full_clinvar_coverage_file_name}'),
                sep='\t',row.names = F,quote = F)
    # Missing table
    missing_clinvar_coverage_file_name<-glue('{target_file_name}.{clinvar_file_name}.clinvar_coverage.missing_variants.tsv')
    write.table(clinvar_coverage_table%>%select(-synonyms)%>%filter(in_target=='not_in_target'),file=glue('./output/clinvar_coverage/{missing_clinvar_coverage_file_name}'),
                sep='\t',row.names = F,quote = F)
    
    
    per_gene_table<-clinvar_coverage_table%>%
      filter(in_target=='not_in_target')%>%
      filter(chr!='chrM')%>%
      group_by(gene)%>%
      summarize(missing_variants_text=paste0(var,collapse='; '),
                missing_variants_id_text=paste0(id,collapse='; '))%>%
      select(gene,missing_variants_text,missing_variants_id_text)
    
    # add the aliases for every gene
    per_gene_table<-per_gene_table%>%left_join(gene_alias%>%select(gene=Symbol,synonyms=Synonyms))
    
    per_gene_name<-glue('{target_file_name}.{clinvar_file_name}.clinvar_coverage.per_gene')
    write.table(per_gene_table,file=glue('./output/clinvar_coverage/{per_gene_name}.tsv'),sep='\t',row.names = F,quote = F)
  }
}

