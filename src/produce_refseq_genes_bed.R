setwd('/media/SSD/Bioinformatics/Projects/gene_annotation_project_2024/')
source('./lib/initialization.R')

# chromosome conversion table
refseq_features<-readr::read_delim('/media/SSD/Bioinformatics/Databases/refseq/GCF_000001405.25_GRCh37.p13_feature_table.txt.gz')
chromosome_conversion_table<-refseq_features%>%distinct(chromosome,genomic_accession)
# gff files
refseq_gff_raw<-readr::read_delim('/media/SSD/Bioinformatics/Databases/refseq/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz',comment = '#',
                              col_names = c('genomic_accession','refseq','region','start','end','X1','strand','X2','INFO'))
refseq_gff<-refseq_gff_raw%>%
  mutate(gene=stringr::str_extract(INFO,'gene_id ([^;]+)',group = 1)%>%stringr::str_replace_all('\"',''),
         transcript=stringr::str_extract(INFO,'transcript_id ([^;]+)',group = 1)%>%stringr::str_replace_all('\"',''),
         exon_num=stringr::str_extract(INFO,'exon_number ([^;]+)',group = 1)%>%stringr::str_replace_all('\"',''))%>%
  select(genomic_accession,region,start,end,gene,transcript,exon_num,strand)%>%
  mutate(size=end-start)%>%
  left_join(chromosome_conversion_table)%>%
  mutate(info=glue('{gene}|{transcript}|{exon_num}'))

# produce gene name transcript file
basic_gene_transcript<-refseq_gff%>%filter(region=='exon')%>%distinct(gene,transcript)
write.table(basic_gene_transcript,'./output/refseq_basic_gene_transcript.tsv',row.names = F,quote = F,sep='\t')
options(scipen = 999)

refseq_exons_bed<-
  refseq_gff%>%
  filter(region=='exon')%>%
  select(chr=chromosome,start,end,info,strand,exon_num,gene,transcript,size)
refseq_exons_bed_file<-'./data/accessory_files/GCF_000001405.25_GRCh37.p13.exons.bed'
write.table(refseq_exons_bed,file = refseq_exons_bed_file,
            row.names = F,quote = F,sep = '\t',col.names = F)

refseq_cds_bed<-
  refseq_gff%>%
  filter(region=='CDS')%>%
  select(chr=chromosome,start,end,info,strand,exon_num,gene,transcript,size)
refseq_cds_bed_file<-'./data/accessory_files/GCF_000001405.25_GRCh37.p13.cds.bed'
write.table(refseq_cds_bed,file = refseq_cds_bed_file,
            row.names = F,quote = F,sep = '\t',col.names = F)

# generate non coding regions
intersection_command<-glue('bedtools subtract -a {refseq_exons_bed_file} -b {refseq_cds_bed_file}')
refseq_non_coding_bed<-read.table(text = system(intersection_command,intern = T))
colnames(refseq_non_coding_bed)<-colnames(refseq_exons_bed)
  
# generate final cds + non coding regions bed file
write.table(
  refseq_cds_bed%>%mutate(type='cds')%>%bind_rows(
    refseq_non_coding_bed%>%mutate(type='non-coding',exon_num=as.character(exon_num))
  ),file = './data/accessory_files/GCF_000001405.25_GRCh37.p13.refseq.bed',
  row.names = F,quote = F,sep = '\t',col.names = F
)

# generate hg19 cds + non coding regions bed file
write.table(
  refseq_cds_bed%>%
    mutate(chr=glue('chr{chr}'),
           type='cds')%>%
    bind_rows(
    refseq_non_coding_bed%>%
      mutate(chr=glue('chr{chr}'),type='non-coding',exon_num=as.character(exon_num))
  ),file = './data/accessory_files/GCF_000001405.25_hg19.p13.refseq.bed',
  row.names = F,quote = F,sep = '\t',col.names = F
)

