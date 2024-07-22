intersect_bed_with_vcf<-function(vcf,bed){
  if (grepl('.gz$',vcf)){
    bedtools_intersect_command<-glue('zcat {vcf} | bedtools intersect -c -a - -b {bed}')
  }else{
    bedtools_intersect_command<-glue('cat {vcf} | bedtools intersect -c -a - -b {bed}')
  }
  # intersect clinvar with the target file
  message(glue('Intersecting bed file: {bed} with VCF: {vcf}.\ncommand:{bedtools_intersect_command}'))
  intersection_output<-read.table(text=system(bedtools_intersect_command,intern = T),sep = '\t',comment.char = '')
  # change the column names
  colnames(intersection_output)<-c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','OVERLAP')
  return(intersection_output)
}

intersect_bed_files<-function(bed1,bed2){
  
}