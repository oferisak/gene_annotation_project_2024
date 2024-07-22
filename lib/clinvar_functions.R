filter_clinvar_vcf_for_plp<-function(clinvar_vcf){
  plp_clinvar_name<-stringr::str_replace(clinvar_vcf,'.vcf','.plp.vcf')
  if(grepl('.gz$',clinvar_vcf)){
    filter_command<-glue('zcat {clinvar_vcf} | grep -i "^#\\|pathogenic" | grep -v Conflict | gzip > {plp_clinvar_name}')
  }else{
    filter_command<-glue('cat {clinvar_vcf} | grep -i "^#\\|pathogenic" | grep -v Conflict > {plp_clinvar_name}')
  }
  message(glue('filtering clinvar vcf {clinvar_vcf} for P/LP variants:\n{filter_command}'))
  system(filter_command)
  message(glue('done. file saved as {plp_clinvar_name}'))
  return(plp_clinvar_name)
}


convert_vcf_to_hg19 <- function(clinvar_file, output_folder) {
    # Extract the base name of the input file
    input_file_basename <- basename(clinvar_file)
    
    # Construct the output file path
    output_file <- file.path(output_folder, paste0("hg19_", sub("\\.gz$", "", input_file_basename)))
    
    # Determine if the input file is gzipped
    read_cmd <- ifelse(grepl("\\.gz$", clinvar_file), "zcat", "cat")
    
    # Open the input and output files
    output_con <- file(output_file, "w")
    clinvar_file_content<-system(glue('{read_cmd} {clinvar_file}'),intern = T)
    # Process the file line by line
    for (line in clinvar_file_content){
      #print(line)
      if (grepl("^##", line) || grepl("^#CHROM", line)) {
        # Write header or column definition lines as is
        writeLines(line, output_con)
      } else if (grepl("^MT", line)) {
        # Replace 'MT' with 'chrM'
        writeLines(stringr::str_replace(line,'^MT','chrM'),output_con)
      } else if (grepl("^[0-9XY]", line) || grepl("^chr[0-9XY]", line)) {
        # Prepend 'chr' to numeric chromosome lines if not already present
        if (!grepl("^chr", line)) {
          writeLines(paste0("chr", line), output_con)
        } else {
          writeLines(line, output_con)
        }
      }
    }
    
    # Close connections
    close(output_con)
    cat("Processing complete. Output saved to", output_file, "\n")
    return(output_file)
}
  
