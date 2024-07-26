library(shinydashboard)
library(shinythemes)
library(DT)
library(shiny)
library(glue)
library(dplyr)
# TARGET: Target and Gene Region Evaluation Tool

# genes_files
genes<-readr::read_delim('./data/refseq_basic_gene_transcript.tsv')
all_gene_symbols<-genes%>%pull(gene)%>%unique()

# pseudogenes file
pseudogenes<-readr::read_delim('./data/pseudogenes.2024-05-17.csv')

# parse target specific data
target_files<-list.files('./data/target_data',full.names = T)
# target vs genes
target_vs_genes<-grep('incomplete_coverage',target_files,value=T)
target_vs_genes_targets<-do.call(rbind.data.frame, basename(target_vs_genes)%>%stringr::str_split('\\.'))[,1]
names(target_vs_genes)<-target_vs_genes_targets
# target vs clinvar missing
target_vs_clinvar_missing<-grep('missing_variants',target_files,value=T)
target_vs_clinvar_missing_targets<-do.call(rbind.data.frame, basename(target_vs_clinvar_missing)%>%stringr::str_split('\\.'))[,1]
names(target_vs_clinvar_missing)<-target_vs_clinvar_missing_targets
# target vs giab stratifications
target_vs_giab<-grep('giab_stratifications.non_zero',target_files,value=T)
target_vs_giab_targets<-do.call(rbind.data.frame, basename(target_vs_giab)%>%stringr::str_split('\\.'))[,1]
names(target_vs_giab)<-target_vs_giab_targets

# use only targets that have all the required information
target_choices <- Reduce(intersect, list(target_vs_genes_targets,target_vs_clinvar_missing_targets,target_vs_giab_targets))

# Define UI ####
ui <- navbarPage(
  theme=shinytheme('flatly'),
  title = 'TARGET',
  # title = tags$div(
  #   class = "navbar-brand",
  #   tags$img(src = "target_logo_1.png", height = "60px", alt = "Logo")  # Replace 'your_image.png' with your image file name
  # ),
  #title = 'TARGET',
  tabPanel("Gene List Upload",
           
             div(selectizeInput('target_selection', label='Select enrichment kit', multiple=F,choices=target_choices,width = 800,options= list(maxOptions = 5000)), 
                 style='font-size:150%;'),
             div(fileInput("gene_list_file", "Select gene list file",
                       multiple = FALSE,width = 800,
                       accept = c(".txt",".csv",".tsv",'.xls','.xlsx')),style='font-size:150%;'),
             selectizeInput('gene_list_selection',choices = NULL,selected=NULL,multiple=TRUE,label='Gene List',width = 800),
             div(verbatimTextOutput(outputId = 'num_of_parsed_genes'),style='font-size:125%;'),
             div(DT::dataTableOutput(outputId = 'gene_list_table'), 
                 style='font-size:150%;')
           
  ),
  
  tabPanel("Target Coverage",
           h3('Gene regions that are not inside the target region'),
           br(),
           div(sliderInput('min_exon_covered_size',label='Show exons in which coverage is lower than this rate:',
                           value=0.8,min = 0,max=1,step = 0.05,width = '500px')),
           div(DT::dataTableOutput(outputId='target_gene_cov_table'),
               style='font-size:100%;')
  ),
  tabPanel("ClinVar P/LP",
           h3('ClinVar P/LP variants that are not inside the target region'),
           br(),
           div(DT::dataTableOutput(outputId='missing_clinvar_table'),
               style='font-size:100%;')
  ),
  tabPanel("Difficult Regions",
           div(sliderInput('min_difficult_region_size',label='Show exons with at least this rate of difficult regions',
                           value=0.75,min = 0,max=1,step = 0.05,width = '500px')),
           div(DT::dataTableOutput(outputId='difficult_regions_table'),
               style='font-size:100%;')
  ),
  
  tabPanel("Pseudogenes",
           h3('Genes with Pseudogene'),
           br(),
           div(DT::dataTableOutput(outputId='pseudogenes_table'),
               style='font-size:100%;')
  ),
  tabPanel("Summary table",
           div(DT::dataTableOutput(outputId='summary_table'),
               style='font-size:100%;')
  ),
  # tags$head(tags$style(HTML("
  #   .navbar-brand {
  #     padding: 0;
  #     margin: 0;
  #   }
  #   .navbar-nav {
  #     margin-left: auto;
  #     margin-right: auto;
  #   }
  #   .navbar {
  #     display: flex;
  #     align-items: center;
  #   }
  # ")))
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  updateSelectizeInput(inputId = 'gene_list_selection', choices = all_gene_symbols, server = TRUE)
  # upload gene list ####
  observeEvent(input$gene_list_file,{
    if (grepl('xls',input$gene_list_file$datapath)){
      gene_list<-readxl::read_excel(input$gene_list_file$datapath,col_names = F)
      colnames(gene_list)[1]<-'gene'
    }else{
      gene_list<-readr::read_delim(input$gene_list_file$datapath,col_names = c('gene'))
    }
    print(gene_list)
    original_genes<-gene_list%>%pull(gene)%>%unique()
    gene_list<-gene_list%>%filter(gene%in%all_gene_symbols)
    updateSelectizeInput(inputId = 'gene_list_selection', choices = all_gene_symbols,selected=gene_list$gene, server = TRUE)
    output$num_of_parsed_genes<-renderText({glue('Out of {length(original_genes)} input genes, identified {nrow(gene_list)} genes. could not find: {paste0(setdiff(original_genes,gene_list%>%pull(gene)),collapse=", ")}')})
  })
  
  gene_or_target_change <- reactive({
    list(input$target_selection,input$gene_list_selection)
  })
  # Target coverage ####
  target_coverage_table <- reactive({
    req(input$target_selection, input$gene_list_selection)
    readr::read_delim(target_vs_genes[[input$target_selection]]) %>%
      rename('rate_covered' = overlap_percent) %>%
      filter(gene %in% input$gene_list_selection) %>%
      filter(rate_covered <= input$min_exon_covered_size) %>%
      mutate(gene = factor(gene),
             type = factor(type),
             transcript = factor(transcript),
             bases_not_covered = exon_size - overlap_size) %>%
      relocate(gene, transcript, .before = chr) %>%
      select(-c(info, num_o_overlapping_regions, strand, exon_size, overlap_size, total_size))
  })
  
  # Difficult regions ####
  base_cols<-c('chr','start','end','info','strand','exon_num','gene','transcript','size','type','non_zero_stratifications')
  difficult_regions_table<-reactive({
    req(input$target_selection, input$gene_list_selection)
    drtable<-readr::read_delim(target_vs_giab[[input$target_selection]])
    dif_regions_cols<-setdiff(colnames(drtable),base_cols)
    drtable%>%
      filter(gene%in%input$gene_list_selection)%>%
      filter(if_any(dif_regions_cols,~.>=input$min_difficult_region_size))%>%
      mutate(gene=factor(gene),
             type=factor(type),
             transcript=factor(transcript))%>%
      relocate(gene,transcript,.before=chr)%>%
      select(-c(info,strand))
  })
  
  # Clinvar variants ####
  clinvar_table<-reactive({
    req(input$target_selection, input$gene_list_selection)
    cv_table<-
      readr::read_delim(target_vs_clinvar_missing[[input$target_selection]])%>%
      filter(gene%in%input$gene_list_selection)%>%
      mutate(gene=factor(gene),
             id=factor(id),
             ref=ifelse(nchar(ref)>30,glue('{substr(ref,1,30)}...'),ref),
             alt=ifelse(nchar(alt)>30,glue('{substr(alt,1,30)}...'),alt),
             var=ifelse(nchar(var)>40,glue('{substr(var,1,40)}...'),var))
    colnames(cv_table)<-tolower(colnames(cv_table))
    cv_table
  })
  
  # Pseudogenes ####
  pseudogenes_table<-reactive({
    req(input$target_selection, input$gene_list_selection)
    pseudogenes%>%filter(gene_symbol%in%input$gene_list_selection) %>%arrange(desc(number_of_pseudogenes))
  })

  
  # Summary table ####
  summary_table<-reactive({
    req(input$target_selection, input$gene_list_selection)
    ## pseudogenes summary table ####
    pseudogene_summary<-pseudogenes_table()%>%
      mutate(pseudogenes_text=ifelse(number_of_pseudogenes==1,
                                     glue('There is {number_of_pseudogenes} pseudogene ({pseudogenes}).'),
                                     glue('There are {number_of_pseudogenes} pseudogenes ({pseudogenes}).')))%>%
      select(gene=gene_symbol,pseudogenes_text)
    ## clinvar summary table ####
    clinvar_summary<-clinvar_table()%>%
      group_by(gene)%>%
      summarize(clinvar_text=ifelse(length(id)>0,
                                    glue('<b>ClinVar IDs not in target</b>: {paste0(id,collapse="; ")}'),''))%>%
      select(gene,clinvar_text)
    
    ## target coverage summary table ####
    exon_coverage_summary<-target_coverage_table()%>%
      group_by(gene,transcript)%>%
      summarize(exon_coverage_text=paste0(glue('exon {exon_num}:{type}:{100*round(rate_covered,3)}'),collapse='; '))%>%
      select(gene,transcript,exon_coverage_text)%>%
      ungroup()%>%group_by(gene)%>%
      summarize(transcript_exon_coverage_text=paste('<b>Exons with incomplete coverage: cds/non-coding: percent covered</b>:<br>',
                                                    paste0(glue('{transcript}: {exon_coverage_text}'),collapse='<br>')))
    
    ## difficult regions summary table
    difficult_regions_summary<-
      difficult_regions_table()%>%
      group_by(gene,transcript)%>%
      summarize(exon_difficult_regions_text=paste0(glue('exon {exon_num}: {non_zero_stratifications}'),collapse='; '))%>%
      ungroup()%>%group_by(gene)%>%
      summarize(difficult_regions_text=paste('<b>Regions difficult to sequence/map</b>:<br>',
                                                    paste0(glue('{transcript}: {exon_difficult_regions_text}'),collapse='<br>')))
    summary_table<-data.frame(gene=input$gene_list_selection)
    if (nrow(pseudogene_summary)>0){
      summary_table<-summary_table%>%left_join(pseudogene_summary)
    }
    if (nrow(clinvar_summary)>0){
      summary_table<-summary_table%>%left_join(clinvar_summary)
    }
    if (nrow(exon_coverage_summary)>0){
      summary_table<-summary_table%>%left_join(exon_coverage_summary)
    }
    if (nrow(difficult_regions_summary)>0){
      summary_table<-summary_table%>%left_join(difficult_regions_summary)
    }
    
    
    summary_text_cols<-grep('_text',colnames(summary_table),value=T)
    
    summary_table<-summary_table%>%
      rowwise() %>%
      filter(any(!is.na(across(all_of(summary_text_cols)))))%>%
      mutate(summary_text = apply(across(all_of(summary_text_cols)), 1, function(x) {
        paste(na.omit(x), collapse = '<br>')
      }))
  
  })
  
  # Target coverage ####
  output$target_gene_cov_table<-renderDT(target_coverage_table(),
                                         filter = list(
                                           position = 'top', clear = FALSE),
                                         options=list(pageLength = nrow(target_coverage_table())))
  # Difficult regions ####
  output$difficult_regions_table<-renderDT(difficult_regions_table()%>%select(-non_zero_stratifications),
                                           filter = list(
                                             position = 'top', clear = FALSE),
                                           options=list(pageLength = nrow(difficult_regions_table()),
                                                        scrollX=TRUE, scrollCollapse=TRUE))
  
  
  #   select(-c(info,non_zero_stratifications,strand))
  output$missing_clinvar_table<-renderDT(clinvar_table(),
                                         filter = list(
                                           position = 'top', clear = FALSE),
                                         options=list(pageLength = nrow(clinvar_table),
                                                      scrollX=TRUE, scrollCollapse=TRUE))
  
  output$pseudogenes_table<-renderDT(pseudogenes_table(),
                                     filter = list(position = 'top', clear = FALSE),
                                     options=list(pageLength = nrow(pseudogenes_table),
                                                  scrollX=TRUE, scrollCollapse=TRUE))
  output$summary_table<-renderDT(datatable(summary_table()%>%select(gene,summary_text),
                                           escape=FALSE,
                                           filter = list(position = 'top', clear = FALSE),
                                           options=list(pageLength = nrow(summary_table()),
                                                        scrollX=TRUE, scrollCollapse=TRUE)))
}
# Run the application 
shinyApp(ui = ui, server = server)