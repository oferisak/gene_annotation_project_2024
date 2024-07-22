library(shinydashboard)
library(shinythemes)
library(DT)
library(shiny)
library(glue)
library(dplyr)

target_info<-list()

# genes_files
genes<-readr::read_delim('./data/refseq_basic_gene_transcript.tsv')
all_gene_symbols<-genes%>%pull(gene)%>%unique()

# pseudogenes file
pseudogenes<-readr::read_delim('./data/pseudogenes.2024-05-17.csv')

# Twist target
target_info[['twist_v2']]<-list()
# clinvar disclaimer
target_info[['twist_v2']][['clinvar']]<-clinvar_data<-readr::read_delim('./data/clinvar_coverage.twist_v2.clinvar_plp_20240708.missing_variants.tsv',delim = '\t')
# target data
target_info[['twist_v2']][['gene_coverage']]<-readr::read_delim('./data/twist_v2_vs_refseq.disclaimer.csv',delim = '\t')
# difficult regions data
target_info[['twist_v2']][['difficult_regions']]<-readr::read_delim('./data/gene_giab_stratifications.disclaimer.csv',delim = '\t')

target_choices<-names(target_info)

# Define UI ####
ui <- navbarPage(
  theme=shinytheme('flatly'),
  title = 'Gene Annotation App',
  
  tabPanel("Gene List Upload",
           fluidRow(
             fileInput("gene_list_file", "Select gene list file",
                       multiple = FALSE,
                       accept = c(".txt",".csv",".tsv",'.xls','.xlsx')),
             selectizeInput('gene_list_selection',choices = NULL,selected=NULL,multiple=TRUE,label='Gene List'),
             div(selectizeInput('target_selection', label='Select enrichment kit', multiple=F,choices=target_choices,width = 800,options= list(maxOptions = 5000)), 
                 style='font-size:150%;'),
             div(verbatimTextOutput(outputId = 'num_of_parsed_genes'),style='font-size:125%;'),
             div(DT::dataTableOutput(outputId = 'gene_list_table'), 
                 style='font-size:125%;')
           )
  ),
  
  tabPanel("Target Coverage",
           div(DT::dataTableOutput(outputId='target_gene_cov_table'),
               style='font-size:100%;')
  ),
  tabPanel("ClinVar P/LP",
           div(DT::dataTableOutput(outputId='missing_clinvar_table'),
               style='font-size:100%;')
  ),
  tabPanel("Difficult Regions",
           div(DT::dataTableOutput(outputId='difficult_regions_table'),
               style='font-size:100%;')
  ),
  
  tabPanel("Pseudogenes",
           div(DT::dataTableOutput(outputId='pseudogenes_table'),
               style='font-size:100%;')
  )
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
  observeEvent(gene_or_target_change(),{
    # Target coverage ####
    target_coverage_table<-target_info[[input$target_selection]][['gene_coverage']]%>%
      filter(gene%in%input$gene_list_selection)%>%
      mutate(gene=factor(gene),
             type=factor(type),
             transcript=factor(transcript),
             bases_not_covered=exon_size-overlap_size)%>%
      rename('rate_covered'=overlap_percent)%>%
      relocate(gene,transcript,.before=chr)%>%
      select(-c(info,num_o_overlapping_regions,strand,exon_size,overlap_size,total_size))
    output$target_gene_cov_table<-renderDT(target_coverage_table,
                                           filter = list(
                                             position = 'top', clear = FALSE),
                                           options=list(pageLength = nrow(target_coverage_table)))
    
    # Difficult regions ####
    difficult_regions_table<-target_info[[input$target_selection]][['difficult_regions']]%>%
      filter(gene%in%input$gene_list_selection)%>%
      mutate(gene=factor(gene),
             type=factor(type),
             transcript=factor(transcript))%>%
      relocate(gene,transcript,.before=chr)%>%
      select(-c(info,non_zero_stratifications,strand))
    output$difficult_regions_table<-renderDT(difficult_regions_table,
                                           filter = list(
                                             position = 'top', clear = FALSE),
                                           options=list(pageLength = nrow(difficult_regions_table),
                                                        scrollX=TRUE, scrollCollapse=TRUE))
    
    # Clinvar variants ####
    clinvar_table<-target_info[[input$target_selection]][['clinvar']]%>%
      filter(gene%in%input$gene_list_selection)%>%
       mutate(gene=factor(gene),
              id=factor(id),
              ref=ifelse(nchar(ref)>30,glue('{substr(ref,1,30)}...'),ref),
              alt=ifelse(nchar(alt)>30,glue('{substr(alt,1,30)}...'),alt),
              var=ifelse(nchar(var)>40,glue('{substr(var,1,40)}...'),var))
    colnames(clinvar_table)<-tolower(colnames(clinvar_table))
    #   select(-c(info,non_zero_stratifications,strand))
    output$missing_clinvar_table<-renderDT(clinvar_table,
                                             filter = list(
                                               position = 'top', clear = FALSE),
                                             options=list(pageLength = nrow(clinvar_table),
                                                          scrollX=TRUE, scrollCollapse=TRUE))
    # Pseudogenes ####
    pseudogenes_table<-pseudogenes%>%filter(gene_symbol%in%input$gene_list_selection) %>%arrange(desc(number_of_pseudogenes))
    output$pseudogenes_table<-renderDT(pseudogenes_table,
                                             filter = list(
                                               position = 'top', clear = FALSE),
                                             options=list(pageLength = nrow(pseudogenes_table),
                                                          scrollX=TRUE, scrollCollapse=TRUE))
  })

}
# Run the application 
shinyApp(ui = ui, server = server)