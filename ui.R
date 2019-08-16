#loading all the R packages
library(shiny)
options(shiny.maxRequestSize = 1000*1024^2) #This command limits the size for file uploads to 1Gb. If unset, the maximum request size defaults to 5MB.
options(shiny.fullstacktrace = TRUE) #Controls whether full stack traces are dumped to the console when errors occur during Shiny app execution. The default is FALSE (pretty stack traces).
  library(tools)
  library(seqinr)
  library(stringr)
  library(tidyverse)
  require(Biostrings)
  library(gridExtra)
  library(shinyBS)
  library(DT)
  library(cowplot)
  library(grid)
  library(shinyWidgets)
  library(reshape)
  library(R.utils)
  library(promises)
  library(future)
  plan(multisession) 
  

ui <- navbarPage("Genome Assembly and Annotation Metrics",
      tabPanel(
          "Welcome", 
           img(src='gqc_white_background_transparent.jpg',  width="200", height="100"),
            helpText(HTML(
                                   '<style>
                             .greenText
                             {
                             color:green;
                             font-weight:bold;
                             }
                             .blackText
                             {
                             color:black;
                            
                             }
                              .blueText
                             {
                             color:blue;
                             
                             }
                              .boldblackText
                             {
                             color:black;
                             font-weight:bold;
                             
                             }
                             </style>
                           <p class="greenText">Welcome to GenomeQC website! </p>
                           <p class="blackText">GenomeQC is a user-friendly and interactive platform that generates descriptive summaries with intuitive graphics for genome assemblies and structural annotations. It also benchmarks user supplied assemblies and annotations against the publicly available reference genomes of their choice. </p>
                            <p class="blackText"> The web application is designed to compute assembly and annotation statistics for small to medium-sized genomes with an upper limit of 2.5 Gb (the approximate size of the maize genome).  </p>
                            <p class="blackText"> The tool’s analysis interface is organized into 3 sections: </p>
		            <p class="boldblackText"> 1. Compare reference genomes </p
                             class="blackText">This section displays various assembly and annotation metrics for the user-selected list of reference genomes. </p>
	                    <p class="boldblackText"> 2. Analyse your genome assembly </p
                             class="blackText">This section provides the user the option to perform analysis on their genome assembly as well as benchmark their analysis with pre-computed reference genomes. </p>
                            <p class="boldblackText"> 3. Analyse your genome annotation </p 
                             class="blackText">This section provides the user the option to perform analysis on their genome annotations as well as benchmark their analysis with the pre-computed reference annotations. </p>           
                              '))
                 ),             
      tabPanel("FAQs", helpText(HTML(
                   '<style>
                             .greenText
                             {
                             color:green;
                             font-weight:bold;
                             }
                             .blackText
                             {
                             color:black;
                            
                             }
                              .blueText
                             {
                             color:blue;
                             
                             }
                             </style>
                              <p class="blackText">Frequently Asked Questions (FAQs): </p>

                               <p class="greenText">Q: What should I do if the web-page gets disconnected? </p
                               class="blackText">A:  The Shiny server is sensitive to internet connectivity, and so, you may experience periodic disconnection of the page. Should this happen, the user can reload the page and resubmit a job.
Once the BUSCO job has been submitted, it is not dependent on an active internet connection.  Even if the webpage is disconnected or if the user closes the webpage/browser, the plots will still be emailed to the email address provided in the input box.</p>
                               <p class="greenText">Q: How long does it take for the BUSCO analysis to complete? </p
                               class="blackText">A: BUSCO analysis of genome assemblies and annotations is a computationally intensive job and the expected run time depends on the size of assemblies and annotation sets. The following lists the expected run time for different genomes: 
Genomes up to 200 Mb: up to 2 hours,
Genomes between 200Mb  and 400 Mb: 3-4 hours,
Genomes between 400Mb and 700 Mb: 4-8 hours,
Genomes between 700Mb and 1.5 Gb: 8-24 hours,
Genomes greater than 1.5 Gb: >1 day</p>
                               <p class="greenText">Q: What should I do if I did not receive BUSCO plots via email? </p
                               class="blackText">A: First, confirm that you clicked on the “Assembly Busco and Contamination Plots” tab.  This step is required in addition to clicking the “Click to Submit your Job” button.  Second, please check your spam folder if you do not find the plots in your inbox. Finally, if you still do not receive the BUSCO plots by email, please contact us with the details of your job submission including your job ID. </p>
                               <p class="blackText">Userguide for this web-application is available at this link: https://github.com/HuffordLab/GenomeQC </span></p>
                               <p class="blackText">Please send questions (after reading the User guide) to: <span class="blueText">john.portwood@ars.usda.gov </span></p>
                       
                                
                               '))
                 ),
                 
  tabPanel("Compare reference genomes",               
    sidebarLayout(
      sidebarPanel(
        HTML(
          "Please click the blue icon on the right of each input field to receive more info about the input field.", 
          "Click again on the icon after reading the info to close the pop-up box", sep = "<br/>"),
        
        tags$style(".popover{
                   max-width: 100%;}"),
        helpText(HTML('<p style="color:blue; fontsize: 9pt">Con: Contigs. Chr: Pseudomolecules</p>')),
        div(
          div(style="width:90%; display:inline-block; vertical-align: middle;", pickerInput(
            inputId = "references", 
            label = HTML('<p> Reference Genomes<span style="color: red; font-size:150%;">*</span></p>'),
            choices = c("MaizeB73_v4_con", "MaizeMo17_CAU_con", "MaizeW22_NRgene_con", "Sorghum_bicolor_NCBIv3_con", "Rice_japonica_Nipponbare4.0_chr", "Arabidopsis_thaliana_TAIR10.1_chr", "Brassica_rapa_v3_chr", "Brachypodium_distachyon_v3_con", "Setaria_italica_v2_chr", "Glycine_max_v2.1_con"),
            options = list(
              `actions-box` = TRUE, 
              size = 10,
              `selected-text-format` = "count > 3"
            ), 
            multiple = TRUE
          )),
          div(style="display:inline-block;vertical-align:top;", bsButton("q1", label = "", icon = icon("info"), style = "info", size = "extra-small"))
        ),
        bsPopover(id = "q1", title = "",
                  
                  content = paste("When the values in this field are selected, the app will",
                                  "provide you comparison of pre-computed reference/database values.", sep = "<br>"),
                  placement = "right",
                  trigger = "click",  options=list(container="body")),
        div(
          div(style="width:90%; display:inline-block; vertical-align: middle;",  textInput("id", label = HTML('<p> Enter your email address<span style="color: red; font-size:150%;">*</span></p>'), value = "")),
          div(style="display:inline-block;vertical-align:top;", bsButton("q2", label = "", icon = icon("info"), style = "info", size = "extra-small"))
        ),
        bsPopover(id = "q2", title = "",
                  content = paste("BUSCO analysis results", 
                                  "will be emailed to this address", sep = "<br>"),
                  placement = "right",
                  trigger = "click",  options=list(container="body")),
        
        actionButton("go", "Click To Submit Your Job", style = "color: blue; font-size:150%; font-weight: bold; height: 50px;
                    width: 300px"),
        
        helpText(HTML('<p style="color:black; font-size: 20px">Download your results</p>')),
        downloadButton('assemblyplot', 'Download NG(X)plot'),
        downloadButton('downloadassemblytable', 'Download asssembly metrics table'),
        downloadButton('downloadannotationtable', 'Download annotation metrics table')
      ),
      mainPanel(
        fluidRow(
          helpText(HTML('<p style="color:black; fontsize: 9pt">For the best results, please click on each tab from left to right one at a time, explore the tab completely, download its results and then move on to next tab. Click the BUSCO tab at the very end only.</p>')),
          helpText(HTML('<p style="color:green; fontsize: 9pt; font-weight: bold">Pop-up plots are available for each metric by clicking on the rows of Assembly and Annotation metrics tables!!!!</p>')),
          tabsetPanel(
            tabPanel("Assembly NG(x) Plot", helpText("NG(X) plots provide information on the contiguity of the assembled genome sequence. Higher the curve, better is the quality of the assembly in terms of contiguity."),  plotOutput("NGX_plot")),
            tabPanel("Assembly Metrics Table", helpText("This tab outputs different metrics including number of scaffolds, L50, N50, LG50, NG50, and gaps percentage values. A good quality assembly will have fewer total scaffolds, higher N50,NG50 values and lower L50, LG50, %N values."), DT::dataTableOutput("lengthmetrics")),
            tabPanel("Annotation Metrics Table", helpText("This tab provides a summary of different annotation metrics like number and average length of gene models, exons, transcripts, etc."), DT::dataTableOutput("annotationmetrics")),
            tabPanel("Assembly and Annotation Busco Plots", helpText("This tab outputs the % of BUSCO genes. A good quality assembly and annotation should have a high number of complete and single copy BUSCOs and a ow number of fragmented and missing BUSCO genes. Please be patient. All the plots will be emailed once the analysis is finished."), htmlOutput("reference_busco_plots"))
            
          )
        )
      )
    )
  ),
        
    tabPanel("Analyze your genome assembly", fluid = TRUE,
             
             sidebarLayout(
               sidebarPanel(
                 HTML(
                   "Please click the blue icon on the right of each input field to receive more info about the input field.", 
                   "Click again on the icon after reading the info to close the pop-up box", sep = "<br/>"),
                 
                 tags$style(".popover{
                   max-width: 100%;}"),
                 helpText(HTML('<p style="color:blue; fontsize: 9pt">Con: Contigs. Chr: Pseudomolecules</p>')),
                 div(
                   div(style="width:90%; display:inline-block; vertical-align: middle;", pickerInput(
                     inputId = "references1", 
                     label = HTML('<p> Reference Genomes</p>'),
                     choices = c("MaizeB73_v4_con", "MaizeMo17_CAU_con", "MaizeW22_NRgene_con", "Sorghum_bicolor_NCBIv3_con", "Rice_japonica_Nipponbare4.0_chr", "Arabidopsis_thaliana_TAIR10.1_chr", "Brassica_rapa_v3_chr", "Brachypodium_distachyon_v3_con", "Setaria_italica_v2_chr", "Glycine_max_v2.1_con"),         
                     options = list(
                       `actions-box` = TRUE, 
                       size = 10,
                       `selected-text-format` = "count > 3"
                     ), 
                     multiple = TRUE
                   )),
                   div(style="display:inline-block;vertical-align:top;", bsButton("q3", label = "", icon = icon("info"), style = "info", size = "extra-small"))
                 ),
                 bsPopover(id = "q3", title = "",
                           
                           content = paste("When the values in this field are selected, the app will",
                                           "provide you comparison of pre-computed reference/database values and also",
                                           "allows you to compare your genome results with these values", sep = "<br>"),
                           placement = "right",
                           trigger = "click",  options=list(container="body")),
                 div(
                   div(style="width:90%; display:inline-block; vertical-align: middle;",  textInput("id1", label = HTML('<p> Enter your email address<span style="color: red; font-size:150%;">*</span></p>'), value = "")),
                   div(style="display:inline-block;vertical-align:top;", bsButton("q4", label = "", icon = icon("info"), style = "info", size = "extra-small"))
                 ),
                 bsPopover(id = "q4", title = "",
                           content = paste("BUSCO and contamination analysis results", 
                                           "will be emailed to this address", sep = "<br>"),
                           placement = "right",
                           trigger = "click",  options=list(container="body")),
                 
        helpText(HTML('<p style="color:red; fontsize: 9pt">Provide just one word for the name</p>')),
        div(
          div(style="display:inline-block;vertical-align:top;",  textInput("user_input1", label = HTML('<p> Name of your genome assembly<span style="color: red; font-size:150%;">*</span></p>'), value = "")),
          div(style="display:inline-block;vertical-align:top;", bsButton("q5", label = "", icon = icon("info"), style = "info", size = "extra-small"))
         ),
        bsPopover(id = "q5", title = "",
                  content = paste("Please provide a one word name for your assembly.", 
                                  "The value in this input field will be used to label",
                                  "the resulting plots and table.",
                                  "", sep = "<br>"),
                  placement = "right",
                  trigger = "click",  options=list(container="body")),

        
        div(
          div(style="display:inline-block;vertical-align:top;", numericInput("estimated_genome_size", label = HTML('<p> Estimated genome size (Mb)<span style="color: red; font-size:150%;">*</span></p>'), value = "")),
          div(style="display:inline-block;vertical-align:top;", bsButton("q6", label = "", icon = icon("info"), style = "info", size = "extra-small"))
        ),
        bsPopover(id = "q6", title = "",
                  content = paste("You need to provide an estimate of the genome size. Please provide the size in megabases (Mb).",
                                  "Example: estimated genome size for maize is 2200 Mb. The value in this input field will be used to",
                                  "calculate NG metrics at different thresholds. Information on NG metrics can be found here:",
                                  "https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics", sep = "<br>"),
                  placement = "right",
                  trigger = "click",  options=list(container="body")),
        
        div(
          div(style="width:90%; display:inline-block; vertical-align: middle;", selectInput("busco_sets1", label = HTML('<p> BUSCO Datasets (select one)<span style="color: red; font-size:150%;">*</span></p>'), c("embryophyta_odb9 (Plants)","bacteria_odb9 (Bacteria)","proteobacteria_odb9 (Bacteria)","rhizobiales_odb9 (Bacteria)","betaproteobacteria_odb9 (Bacteria)","gammaproteobacteria_odb9 (Bacteria)","enterobacteriales_odb9 (Bacteria)","deltaepsilonsub_odb9 (Bacteria)","actinobacteria_odb9 (Bacteria)","cyanobacteria_odb9 (Bacteria)","firmicutes_odb9 (Bacteria)","clostridia_odb9 (Bacteria)","lactobacillales_odb9 (Bacteria)","bacillales_odb9 (Bacteria)","bacteroidetes_odb9 (Bacteria)","spirochaetes_odb9 (Bacteria)","tenericutes_odb9 (Bacteria)", "eukaryota_odb9 (Eukaryota)","fungi_odb9 (Fungi)","microsporidia_odb9 (Fungi)","dikarya_odb9 (Fungi)","ascomycota_odb9 (Fungi)","pezizomycotina_odb9 (Fungi)","eurotiomycetes_odb9 (Fungi)","sordariomyceta_odb9 (Fungi)","saccharomyceta_odb9 (Fungi)","saccharomycetales_odb9 (Fungi)","basidiomycota_odb9 (Fungi)","metazoa_odb9 (Metazoa)","nematoda_odb9 (Roundworms)","arthropoda_odb9 (Insects)","vertebrata_odb9 (Vertebrates)","mammalia_odb9 (Mammals)","protists_ensembl (Protist)"))),
          div(style="display:inline-block;vertical-align:top;", bsButton("q7", label = "", icon = icon("info"), style = "info", size = "extra-small"))
        ),
        bsPopover(id = "q7", title = "",
                  
                  content = paste("This field is required to calculate the percentage of Benchmark Universal Single-Copy Orthologs(BUSCO)",
                                  "in your genome assembly to provide intuitive quantitative measures of genomic data completeness in",
                                  "terms of expected gene content. BUSCO datasets are lineage-specific profiles. Genes in these datasets",
                                  "are selected from orthologous groups with genes present as single-copy orthologs in at least 90% of", 
                                  "the species. For more info visit this link: https://busco.ezlab.org/", sep = "<br>"),
                  placement = "right",
                  trigger = "click",  options=list(container="body")),
        div(
          div(style="width:90%; display:inline-block; vertical-align: middle;", selectInput("augustus_species", label = HTML('<p> AUGUSTUS Species (select one)<span style="color: red; font-size:150%;">*</span></p>'), c("maize","adorsata","aedes","amphimedon","ancylostoma_ceylanicum","anidulans","arabidopsis","aspergillus_fumigatus","aspergillus_nidulans","aspergillus_oryzae","aspergillus_terreus","bombus_impatiens1","bombus_terrestris2","botrytis_cinerea","b_pseudomallei","brugia","cacao","caenorhabditis","camponotus_floridanus","candida_albicans","candida_guilliermondii","candida_tropicalis","c_elegans_trsk","chaetomium_globosum","chicken","chlamy2011","chlamydomonas","chlorella","coccidioides_immitis","Conidiobolus_coronatus","coprinus","coprinus_cinereus","coyote_tobacco","cryptococcus","cryptococcus_neoformans_gattii","cryptococcus_neoformans_neoformans_B","cryptococcus_neoformans_neoformans_JEC21","culex","debaryomyces_hansenii","E_coli_K12","elephant_shark","encephalitozoon_cuniculi_GB","eremothecium_gossypii","fly","fusarium","fusarium_graminearum","galdieria","generic","heliconius_melpomene1","histoplasma","histoplasma_capsulatum","honeybee1","human","japaneselamprey","kluyveromyces_lactis","laccaria_bicolor","leishmania_tarentolae","lodderomyces_elongisporus","magnaporthe_grisea","maize5","nasonia","neurospora","neurospora_crassa","parasteatoda","pchrysosporium","pea_aphid","pfalciparum","phanerochaete_chrysosporium","pichia_stipitis","pneumocystis","rhizopus_oryzae","rhodnius","rice","saccharomyces","saccharomyces_cerevisiae_rm11-1a_1","saccharomyces_cerevisiae_S288C","s_aureus","schistosoma","schistosoma2","schizosaccharomyces_pombe","sealamprey","s_pneumoniae","sulfolobus_solfataricus","template_prokaryotic","tetrahymena","thermoanaerobacter_tengcongensis","tomato","toxoplasma","tribolium2012","trichinella","ustilago","ustilago_maydis","verticillium_albo_atrum1","verticillium_longisporum1","volvox","wheat","Xipophorus_maculatus","yarrowia_lipolytica","zebrafish"))),
          div(style="display:inline-block;vertical-align:top;", bsButton("q8", label = "", icon = icon("info"), style = "info", size = "extra-small"))),
        
        bsPopover(id = "q8", title = "",
                  
                  content = paste("This field is also required to calculate the BUSCO scores for your genome", 
                                  "assembly. Currently, AUGUSTUS (gene prediction program used in BUSCO pipeline)",
                                  "has been trained for predicting genes in this species list. If your species is not", 
                                  "listed in the drop-down menu, please select the most closely related species. More", 
                                  "info here: http://bioinf.uni-greifswald.de/augustus/", sep = "<br>"),
                  placement = "right",
                  trigger = "click",  options=list(container="body")),
        
        

        helpText(HTML('<p style="color:red; fontsize: 9pt">Maximum upload limit for genome fasta file is 1Gb (compressed file size)</p>')),
        
        
        div(
          div(style="width:90%; display:inline-block; vertical-align: middle;", fileInput("fasta_file1", label = HTML('<p> Upload Genome Fasta File (.gz format)<span style="color: red; font-size:150%;">*</span></p>'),  accept = c(".gz",".zip"), multiple = FALSE)),
          div(style="display:inline-block; vertical-align: top;", bsButton("q9", label = "", icon = icon("info"), style = "info", size = "extra-small", shape = "circular"))
        ),
        
        bsPopover(id = "q9", title = "",
                  
                  content = paste("Please compress your genome FASTA file and upload the",
                                  ".gz compression file only. Info on this format", 
                                  "here: https://en.wikipedia.org/wiki/Gzip", sep = "<br>"),
                  placement = "top",
                  trigger = "click",  options=list(container="body")),
        actionButton("go1", "Click To Submit Your Job", style = "color: blue; font-size:150%; font-weight: bold; height: 50px;
                        width: 300px"),
        
        helpText(HTML('<p style="color:black; font-size: 20px">Download your results</p>')),
    
        downloadButton('assemblyplot1', 'Download NG(X)plot'),
        downloadButton('downloadassemblytable1', 'Download asssembly metrics table')
        
               ),
        mainPanel(
          fluidRow(
            helpText(HTML('<p style="color:black; fontsize: 9pt">For the best results, please click on each tab from left to right one at a time, explore the tab completely, download its results and then move on to next tab. Click the BUSCO tab at the very end only.</p>')),
            helpText(HTML('<p style="color:green; fontsize: 9pt; font-weight: bold">Pop-up plots are available for each metric by clicking on the rows of Assembly and Annotation metrics tables!!!!</p>')),
            tabsetPanel(
              tabPanel("Assembly NG(x) Plot", helpText("NG(X) plots provide information on the contiguity of the assembled genome sequence. Higher the curve, better is the quality of the assembly in terms of contiguity. Please note: if you uploaded a large genome you might experience a short lag period between when you click this tab and when your job gets started."),  plotOutput("NGX_plot1")),
              tabPanel("Assembly Metrics Table", helpText("This tab outputs different metrics including number of scaffolds, L50, N50, LG50, NG50, and gaps percentage values. A high quality assembly will have fewer total scaffolds, higher N50, NG50 values and lower L50, LG50, %N values. Please note: if you uploaded a big genome you might experience a short lag period between when you click this tab and when your job gets started."), DT::dataTableOutput("lengthmetrics1")),
              tabPanel("Assembly Busco and Contamination Plots", helpText("This tab outputs the % of BUSCO genes.  A good quality assembly and annotation should have a higher number of complete and single copy BUSCOs and a lower number of fragmented and missing BUSCO genes. For contamination analysis, the assembled sequences are screened against the NCBI UniVec database to quickly identify sequences of vector origin or those of adaptors or linkers.  Please note: you might experience a short lag period between when you click this tab and when your BUSCO job gets submitted to the server. Please be patient. All the plots will be emailed once the analysis is finished. Please see the FAQ section for an estimate of expected run time based on the size of the submitted assembly.
"), htmlOutput("assembly_busco_plot"))
              
            )
          )
        )
      )
    ),
    
    tabPanel("Analyze your genome annotation", fluid = TRUE,
             
             sidebarLayout(
               sidebarPanel(
                 HTML(
                   "Please click the blue icon on the right of each input field to receive more info about the input field.", 
                   "Click again on the icon after reading the info to close the pop-up box", sep = "<br/>"),
                 tags$style(".popover{
                            max-width: 100%;}"),
                 helpText(HTML('<p style="color:blue; fontsize: 9pt">Con: Contigs. Chr: Pseudomolecules</p>')),
                 div(
                   div(style="width:90%; display:inline-block; vertical-align: middle;", pickerInput(
                     inputId = "references2", 
                     label = HTML('<p> Reference Genomes</p>'),
                     choices = c("MaizeB73_v4_con", "MaizeMo17_CAU_con", "MaizeW22_NRgene_con", "Sorghum_bicolor_NCBIv3_con", "Rice_japonica_Nipponbare4.0_chr", "Arabidopsis_thaliana_TAIR10.1_chr", "Brassica_rapa_v3_chr", "Brachypodium_distachyon_v3_con", "Setaria_italica_v2_chr", "Glycine_max_v2.1_con"),
                      options = list(
                       `actions-box` = TRUE, 
                       size = 10,
                       `selected-text-format` = "count > 3"
                     ), 
                     multiple = TRUE
                   )),
                   div(style="display:inline-block;vertical-align:top;", bsButton("q10", label = "", icon = icon("info"), style = "info", size = "extra-small"))
                 ),
                 bsPopover(id = "q10", title = "",
                           
                           content = paste("When the values in this field are selected, the app will",
                                           "provide you comparison of pre-computed reference/database values and also",
                                           "allows you to compare your genome results with these values", sep = "<br>"),
                           placement = "right",
                           trigger = "click",  options=list(container="body")),
                 div(
                   div(style="width:90%; display:inline-block; vertical-align: middle;",  textInput("id2", label = HTML('<p> Enter your email address<span style="color: red; font-size:150%;">*</span></p>'), value = "")),
                   div(style="display:inline-block;vertical-align:top;", bsButton("q11", label = "", icon = icon("info"), style = "info", size = "extra-small"))
                 ),
                 bsPopover(id = "q11", title = "",
                           content = paste("BUSCO analysis results will", 
                                           "be emailed to this address", sep = "<br>"),
                           placement = "right",
                           trigger = "click",  options=list(container="body")),

                  helpText(HTML('<p style="color:red; fontsize: 9pt">Provide just one word for the name</p>')),
                  div(
                   div(style="display:inline-block;vertical-align:top;",  textInput("user_input2", label = HTML('<p> Name of your genome annotation<span style="color: red; font-size:150%;">*</span></p>'), value = "")),
                   div(style="display:inline-block;vertical-align:top;", bsButton("q12", label = "", icon = icon("info"), style = "info", size = "extra-small"))
                 ),
                 bsPopover(id = "q12", title = "",
                           content = paste("Please provide a one word name for your assembly.", 
                                           "The value in this input field will be used to label",
                                           "the resulting plots and table.",
                                           "", sep = "<br>"),
                           placement = "right",
                           trigger = "click",  options=list(container="body")),
                 
        
                 div(
                   div(style="width:90%; display:inline-block; vertical-align: middle;", selectInput("busco_sets2", label = HTML('<p> BUSCO Datasets (select one)<span style="color: red; font-size:150%;">*</span></p>'), c("embryophyta_odb9 (Plants)","bacteria_odb9 (Bacteria)","proteobacteria_odb9 (Bacteria)","rhizobiales_odb9 (Bacteria)","betaproteobacteria_odb9 (Bacteria)","gammaproteobacteria_odb9 (Bacteria)","enterobacteriales_odb9 (Bacteria)","deltaepsilonsub_odb9 (Bacteria)","actinobacteria_odb9 (Bacteria)","cyanobacteria_odb9 (Bacteria)","firmicutes_odb9 (Bacteria)","clostridia_odb9 (Bacteria)","lactobacillales_odb9 (Bacteria)","bacillales_odb9 (Bacteria)","bacteroidetes_odb9 (Bacteria)","spirochaetes_odb9 (Bacteria)","tenericutes_odb9 (Bacteria)", "eukaryota_odb9 (Eukaryota)","fungi_odb9 (Fungi)","microsporidia_odb9 (Fungi)","dikarya_odb9 (Fungi)","ascomycota_odb9 (Fungi)","pezizomycotina_odb9 (Fungi)","eurotiomycetes_odb9 (Fungi)","sordariomyceta_odb9 (Fungi)","saccharomyceta_odb9 (Fungi)","saccharomycetales_odb9 (Fungi)","basidiomycota_odb9 (Fungi)","metazoa_odb9 (Metazoa)","nematoda_odb9 (Roundworms)","arthropoda_odb9 (Insects)","vertebrata_odb9 (Vertebrates)","mammalia_odb9 (Mammals)","protists_ensembl (Protist)"))),
                   div(style="display:inline-block;vertical-align:top;", bsButton("q13", label = "", icon = icon("info"), style = "info", size = "extra-small"))
                 ),
                 bsPopover(id = "q13", title = "",
                           
                           content = paste("This field is required to calculate the percentage of Benchmark Universal Single-Copy Orthologs(BUSCO)",
                                           "in your genome annotation set to provide intuitive quantitative measures of genomic data completeness in",
                                           "terms of expected gene content. BUSCO datasets are lineage-specific profiles. Genes in these datasets",
                                           "are selected from orthologous groups with genes present as single-copy orthologs in at least 90% of", 
                                           "the species. For more info visit this link: https://busco.ezlab.org/", sep = "<br>"),
                           placement = "right",
                           trigger = "click",  options=list(container="body")),

                 helpText(HTML('<p style="color:red; fontsize: 9pt">Maximum upload limit for genome fasta file is 1Gb (compressed file size)</p>')),
        
        
                 div(
                   div(style="width:90%; display:inline-block; vertical-align: middle;", fileInput("fasta_file2", label = HTML('<p> Upload Genome Fasta File (.gz format)<span style="color: red; font-size:150%;">*</span></p>'),  accept = c(".gz",".zip"), multiple = FALSE)),
                   div(style="display:inline-block; vertical-align: top;", bsButton("q14", label = "", icon = icon("info"), style = "info", size = "extra-small", shape = "circular"))
                 ),
                 
                 bsPopover(id = "q14", title = "",
                           
                           content = paste("Please compress your genome fasta file and upload the",
                                           ".gz compression file only. Info on this format", 
                                           "here: https://en.wikipedia.org/wiki/Gzip", sep = "<br>"),
                           placement = "top",
                           trigger = "click",  options=list(container="body")),
        

        helpText(HTML('<p style="color:red; fontsize: 1pt">Busco analysis of gene structural annotations requires that every contig or chromosome name found in the 1st column of the input GFF file must have a corresponding sequence entry in the fasta file as shown here: 
https://www.ncbi.nlm.nih.gov/genome/167?genome_assembly_id=161521 .As an alternative, you could upload a corresponding transcript fasta file in addition to the GFF file.</p>')),
        
        div(
          div(style="width:90%; display:inline-block; vertical-align: middle;", fileInput("gff_file", label = HTML('<p> Upload Structure Annotation File (gff, gff3 or gtf format in .gz compression)<span style="color: red; font-size:150%;">*</span></p>'),
                                                                                          accept = ".gz")),
          div(style="display:inline-block;vertical-align:top;", bsButton("q15", label = "", icon = icon("info"), style = "info", size = "extra-small"))),
          bsPopover(id = "q15", title = "",
                  
                  content = paste("Please compress your genome annotation file and upload the",
                                  ".gz compression file only.",
                                  "Visit this link for more info on structure annotation file",
                                  "formats: https://useast.ensembl.org/info/website/upload/gff3.html", sep = "<br>"),
                  placement = "top",
                  trigger = "click",  options=list(container="body")),
        div(
          div(style="width:90%; display:inline-block; vertical-align: middle;", fileInput("transcripts_file", "Upload Transcripts Fasta File (.gz compression format)",
                                                                                          accept = ".gz")),
          div(style="display:inline-block;vertical-align:top;", bsButton("q16", label = "", icon = icon("info"), style = "info", size = "extra-small"))),
          bsPopover(id = "q16", title = "",
                  
                  content = paste("This field is required to calculate annotation BUSCO scores.",
                                  "Currently the app is configured to use the information from transcripts",
                                  "fasta file if the user uploads it. if the user does not upload the",
                                  "transcripts file, the app will check if the information in the first",
                                  "column of the gff file corresponds to the headers in the fasta file as shown",
                                  "in the example. If there is a discrepancy, it will print an error message;",
                                  "otherwise, the job will be submitted.", sep = "<br>"),
                  placement = "top",
                  trigger = "click",  options=list(container="body")),
        
        actionButton("go2", "Click To Submit Your Job", style = "color: blue; font-size:150%; font-weight: bold; height: 50px;
                        width: 300px"),
        
        helpText(HTML('<p style="color:black; font-size: 20px">Download your results </p>')),
        downloadButton('downloadannotationtable1', 'Download annotation metrics table')
               ),
        mainPanel(
          fluidRow(
            helpText(HTML('<p style="color:black; fontsize: 9pt">For the best results, please click on each tab from left to right one at a time, explore the tab completely, download its results and then move on to next tab. Click the BUSCO tab at the very end only.</p>')),
            helpText(HTML('<p style="color:green; fontsize: 9pt; font-weight: bold">Pop-up plots are available for each metric by clicking on the rows of Assembly and Annotation metrics tables!!!!</p>')),
            tabsetPanel(
              tabPanel("Annotation Metrics Table", helpText("This tab provides a summary of different annotation metrics like number and average length of gene models, exons, transcripts, etc. Please note: if you uploaded a big genome you might experience a short lag period between when you click this tab and when your job gets started."), DT::dataTableOutput("annotationmetrics1")),
              tabPanel("Annotation Busco Plot", helpText("This tab outputs the % of BUSCO genes. A good quality assembly and annotation should have a higher number of complete and single copy BUSCOs and a lower number of fragmented and missing BUSCO genes.  Please note: you might experience a short lag period between when you click this tab and when your BUSCO job gets submitted to the server. Please be patient. All the plots will be emailed once the analysis is finished. Please see the FAQ section for an estimate of expected run time based on the size of the submitted annotation."), htmlOutput("annotation_busco_plot"))
            )
          )
        )
      )
    )
)

        
        
        
        
        
        
        
 
