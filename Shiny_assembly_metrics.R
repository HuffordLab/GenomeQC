Sys.setenv (PATH="/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/home/nancy/.local/bin:/home/nancy/bin:/home/nancy/assembly_app/augustus-3.2.3/bin")
Sys.setenv(AUGUSTUS_CONFIG_PATH="/home/nancy/assembly_app/augustus-3.2.3/config/")
#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
setwd("/Users/nancy/busco-master/scripts/") # SETTING THE WORKING DIRECTORY OF R TO THE BUSCO SCIRPTS DIRECTORY


library(shiny)
if (interactive()) {
  library(seqinr)
  library(stringr)
  library(tidyverse)
  require(Biostrings)
  library(gridExtra)
  
  
  ui <- fluidPage(
    titlePanel("Assembly Metrics"),
    sidebarLayout(
      sidebarPanel(
        numericInput("estimated_genome_size", label = h3("estimated genome size"), value = ""),
        selectInput("busco_sets", "Busco Datasets (select one)", c("bacteria_odb9","proteobacteria_odb9","rhizobiales_odb9","betaproteobacteria_odb9","gammaproteobacteria_odb9","enterobacteriales_odb9","deltaepsilonsub_odb9","actinobacteria_odb9","cyanobacteria_odb9","firmicutes_odb9","clostridia_odb9","lactobacillales_odb9","bacillales_odb9","bacteroidetes_odb9","spirochaetes_odb9","tenericutes_odb9","eukaryota_odb9","fungi_odb9","microsporidia_odb9","dikarya_odb9","ascomycota_odb9","pezizomycotina_odb9","eurotiomycetes_odb9","sordariomyceta_odb9","saccharomyceta_odb9","saccharomycetales_odb9","basidiomycota_odb9","metazoa_odb9","nematoda_odb9","arthropoda_odb9","insecta_odb9","endopterygota_odb9","hymenoptera_odb9","diptera_odb9","vertebrata_odb9","actinopterygii_odb9","tetrapoda_odb9","aves_odb9","mammalia_odb9","euarchontoglires_odb9","laurasiatheria_odb9","embryophyta_odb9","protists_ensembl","alveolata_stramenophiles_ensembl")),
        fileInput("fasta_file", "Choose a Fasta File",
                  accept = c(".fasta", 
                             ".fa",
                             ".fna"),
                  multiple = FALSE)
        
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("NG(x) Plot", plotOutput("NX_plot")), 
          tabPanel("Length Metrics Table", tableOutput("length_metrics")),
          tabPanel("Busco Plot", plotOutput("busco_plot"))
      )
    )
  )
)
  server <- function(input, output,session) {
    
    # THIS MODULE PRODUCES NG(X) PLOT GIVEN THE ESTIMATED GENOME SIZE AND COMPARE THE VALUES WITH PRE-COMPUTED REFERENCE VALUES
    
    output$NX_plot <- renderPlot({
      inFile <- input$fasta_file
      
      if (is.null(inFile))
        return(NULL)
      
      
      df <- read.fasta(inFile$datapath, as.string = TRUE)
      progress <- Progress$new(session, min=1, max=15)
      on.exit(progress$close())
      
      progress$set(message = 'Calculation in progress',
                   detail = 'This may take a while...')
      
      for (i in 1:15) {
        progress$set(value = i)
        Sys.sleep(0.5)
      }
      
      le <- getLength(df)
      
      
      sorted <- sort(le, decreasing = TRUE)
      NX_scaffoldlength <- vector()
      for (i in 1:100)
      {
        
        n <- sorted[cumsum(sorted) >= as.numeric(input$estimated_genome_size)*i/100][1]
        NX_scaffoldlength <- c(NX_scaffoldlength,n)
        n <- NULL
      }
      NX_threshold <- c(1:100)
      B73_length <- c(205271, 205271, 205271, 205271, 205271, 205271, 174301, 174301, 174301,
                      174301, 174301, 116448, 116448, 116448, 116448, 106448, 109864, 109864,
                      109864, 83770, 83770, 83770, 83221, 83221, 83221, 74831, 74831,
                      74831, 71404, 71404, 71404, 70286, 70286, 65062, 65062, 65062,
                      61621, 61621,  50068,  50068,  50062,  50062,  50062,  48144,  48144,
                      45430,  45430,  39295,  39295,  38994,  38994,  34516,  30911,  30911,
                      29382,  29382,  27184,  27181,  27181,  20348,  20218,  20218,  19880,
                      19687,  19338,  19338,  18133,  18433,  18118,  18535,  17731,  17731,
                      17116,  16894,  12734,  11050,  10236,  8618,  8596,  8454,  8023,
                      7842,  7646,  7534,  7441,  7154,  6445,  6245,  6112,  6023,
                      4857,  4384,  4253,  4115,   3984,   3564,   3334,   3275,   2008,
                      509)
      xdata <- NX_threshold
      y1 <- log(B73_length)
      y2 <- log(NX_scaffoldlength)
      NXdata <- data.frame(xdata, y1, y2)
      ggplot (NXdata, aes(x = xdata)) +
        geom_line(aes(y=y1, colour = "B73")) +
        geom_line(aes(y=y2, colour = "user_input")) +
        ylab("log scaled scaffold_length")+
        xlab("N(X)") +
        scale_x_continuous(breaks = seq(0, 100, by = 10)) +
        geom_vline(xintercept = 50)
    })
    
    
    # THIS MODULE PRODUCES METRICS TABLE WITH NG VALUES GIVEN THE ESTIMATED GENOME SIZE AND COMPARE THE VALUES WITH PRE-COMPUTED REFERENCE VALUES
    
    output$length_metrics <- renderTable({
      inFile <- input$fasta_file
      
      if (is.null(inFile))
        return(NULL)
      
      
      df <- read.fasta(inFile$datapath, as.string = TRUE)
      progress <- Progress$new(session, min=1, max=20)
      on.exit(progress$close())
      
      progress$set(message = 'Calculation in progress',
                   detail = 'This may take a while...')
      
      for (i in 1:20) {
        progress$set(value = i)
        Sys.sleep(0.5)
      }
      len <- getLength(df)
      num_scaffolds <- length(len)
      total_size_scaffolds <- sum(len)
      total_scaffold_length_as_precentage_genome_size <- (total_size_scaffolds/input$estimated_genome_size)*100
      scaffolds_greater_than_25k <- len[len>25000]
      useful_amount_scaffold_sequence <- sum(scaffolds_greater_than_25k)
      percentage_estimated_genome_useful <- (useful_amount_scaffold_sequence/input$estimated_genome_size)*100
      longest_scaffold <- max(len)
      shortest_scaffold <- min(len)
      scaffolds_greater_than_1k <- length(len[len>1000])
      scaffolds_greater_than_10k <- length(len[len>10000])
      scaffolds_greater_than_100k <- length(len[len>100000])
      scaffolds_greater_than_1M <- length(len[len>1000000])
      scaffolds_greater_than_10M <- length(len[len>10000000])
      sorted <- sort(len, decreasing = TRUE)
      N50 <- sorted[cumsum(sorted) >= sum(sorted)*50/100][1]
      L50 <- length(len[len>=N50])
      NG50 <- sorted[cumsum(sorted) >= (input$estimated_genome_size)*50/100][1]
      LG50 <- length(len[len>=NG50])
      seq_fasta <- readDNAStringSet(inFile$datapath,format="fasta")
      tab <- alphabetFrequency(seq_fasta, baseOnly=T)
      colSums(tab)
      percentage_A <- (colSums(tab)[1]/total_size_scaffolds)*100
      percentage_C <- (colSums(tab)[2]/total_size_scaffolds)*100
      percentage_G <- (colSums(tab)[3]/total_size_scaffolds)*100
      percentage_T <- (colSums(tab)[4]/total_size_scaffolds)*100
      percentage_N <- (colSums(tab)[5]/total_size_scaffolds)*100
      N_metrics <- matrix(c(27432, 1715997, 78, 15981417, 72, 209044, 64, 27251, 20481, 5713, 4, 0, 130512, 3967, 94276, 6150, 26, 23, 23, 26, 0.79, num_scaffolds, total_size_scaffolds, total_scaffold_length_as_precentage_genome_size, useful_amount_scaffold_sequence, percentage_estimated_genome_useful, longest_scaffold, shortest_scaffold, scaffolds_greater_than_1k, scaffolds_greater_than_10k, scaffolds_greater_than_100k, scaffolds_greater_than_1M, scaffolds_greater_than_10M, N50, L50, NG50, LG50, percentage_A, percentage_C, percentage_G, percentage_T, percentage_N), ncol = 2)
      colnames(N_metrics) <- c("B73", "user_input")
      row.names(N_metrics) <- c ("Number of scaffolds", "Total size of scaffolds", "Total scaffold length as percenatge of assumed genome size", "useful amount of scaffold sequences (>=25K nt)", "% of estimated genome that is useful", "Longest scaffold", "Shortest scaffold", "Number of scaffolds > 1K nt", "Number of scaffolds > 10K nt", "Number of scaffolds > 100K nt", "Number of scaffolds > 1M nt", "Number of scaffolds > 10M nt", "N50", "L50", "NG50", "LG50", "%A", "%C", "%G", "%T", "%N")
      
      N_metrics
    }, rownames = TRUE)
    
    
    
    # THIS MODULE PRODUCES BUSCO PLOT AND COMPARE THE VALUES WITH PRE-COMPUTED REFERENCE VALUES
    
    
    output$busco_plot <- renderPlot({
      inFile <- input$fasta_file
      
      if (is.null(inFile))
        return(NULL)
      example_file.fasta <- inFile$datapath
      output_name <- input$fasta_file
      busco_name <- input$busco_sets
      progress <- Progress$new(session, min=1, max=20)
      on.exit(progress$close())
      
      progress$set(message = 'Calculation in progress',
                   detail = 'This may take a while...')
      
      for (i in 1:20) {
        progress$set(value = i)
        Sys.sleep(0.5)
      }
      
      system(paste("python /home/nancy/assembly_app/busco-master/scripts/run_BUSCO.py --in", example_file.fasta, "--out", output_name, "-m genome --lineage", busco_name, "--species=E_coli_K12"), intern = TRUE)
      name_busco_directory <- paste0("run_", inFile)[1]
      file_name <- paste0("short_summary_", inFile, ".txt")[1]
      summary_file <- paste0(name_busco_directory,"/",file_name)
      m <- system(paste("tac", summary_file, "| sed '8q;d'"), intern = TRUE)
      user_input <- as.numeric(str_match_all(m,"[SDFM]:([0-9]+?\\.[0-9]+?)")[[1]][,2])
      B73 <- c(92, 2, 3, 3)
      d = tibble(B73,user_input, category = as.factor(c("C&S", "C&D", "F", "M")))
      plt <- d %>%
        gather(Group, value, -category) %>%
        ggplot(aes(x = Group, y = value, fill = category, stat = "identity")) +
        geom_col() +
        theme_bw() +
        theme(legend.position = "bottom") +
        theme(aspect.ratio = 2/1) +
        xlab("genomes") +
        ylab("percenatge_busco_genes")
      C_S <- system(paste("tac", summary_file, "| sed '5q;d'"), intern = TRUE)
      C_S_number <- as.numeric(str_match(C_S,"[0-9]+")[1])
      C_D <- system(paste("tac", summary_file, "| sed '4q;d'"), intern = TRUE)
      C_D_number <- as.numeric(str_match(C_D,"[0-9]+")[1])
      F_R <- system(paste("tac", summary_file, "| sed '3q;d'"), intern = TRUE)
      F_R_number <- as.numeric(str_match(F_R,"[0-9]+")[1])
      M_S <- system(paste("tac", summary_file, "| sed '2q;d'"), intern = TRUE)
      M_S_number <- as.numeric(str_match(M_S,"[0-9]+")[1])
      user_busco_numbers <- c(C_S_number, C_D_number, F_R_number, M_S_number)
      user_busco_numbers
      B73_busco_numbers <- c(137, 3, 4, 4)
      busco_table <- matrix(c(B73_busco_numbers, user_busco_numbers), ncol = 2, nrow = 4)
      colnames(busco_table) <- c("B73", "user_input")
      row.names(busco_table) <- c("C&S","C&D","F","M")
      tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
      tb1 <- tableGrob(busco_table, rows = rownames(busco_table), theme = tt)
      grid.arrange(plt, tb1, nrow=1, as.table=TRUE)
      name_busco_directory <- paste0("run_", inFile)[1]
      system(paste("rm -r", name_busco_directory))
    })
    
    
  }
  
  shinyApp (ui, server)
}