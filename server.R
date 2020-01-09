
library(shiny)

server <- function(input, output,session) {

  observeEvent(input$go, {
    showModal(modalDialog
      ("Thanks for submitting your job! Please click the tabs now")  
      )
      })  
  
  
  observeEvent(input$go1, {
    showModal(modalDialog(
      "Thanks for submitting your job! Please click the tabs now"   
    ))
  })  
  
  observeEvent(input$go2, {
    showModal(modalDialog(
      "Thanks for submitting your job! Please click the tabs now" 
    ))
  }) 
  
  ran_number <- eventReactive(input$go,{
    
   # uniq_num <- round(runif(1,0,1000),2)  # generates job Id for reference genome tabs
     u_num <- sample.int(1000000:999999, 1)
     uniq_num <-  as.hexmode(u_num)
     return(uniq_num)
  })
  
  ran_number1 <- eventReactive(input$go1,{
         u_num <- sample.int(1000000:999999, 1)
     uniq_num <-  as.hexmode(u_num)
    #uniq_num <- round(runif(1,0,1000),2)  # generates job Id for genome assembly tabs
    
     return(uniq_num)
  })
  
  ran_number2 <- eventReactive(input$go2,{
         u_num <- sample.int(1000000:999999, 1)
     uniq_num <-  as.hexmode(u_num)
   # uniq_num <- round(runif(1,0,1000),2) # generates job Id for genome annotation tabs
    return(uniq_num)
  })

 
  input_fasta_file1 <- eventReactive(input$go1,{
    ext <- file_ext(input$fasta_file1$datapath)
    if (ext == 'gz')
    {input_f <- gunzip(input$fasta_file1$datapath)[1]}
    else
    {input_f <- unzip(input$fasta_file1$datapath)}
    u_num <- ran_number1()
    file_f <- paste0("file_fasta1", u_num) # rename of uploaded assembly file in the genome assembly tab
    destination_directory <- paste0("/opt/shiny-server/samples/sample-apps/assembly_statistics/uploaded_files/",file_f)  
    file.copy(input_f, destination_directory)  # copy the renamed uploaded file to a directory in the shiny server folder
    return (destination_directory)
    
  })
 
    input_fasta_file3 <- eventReactive(input$go1,{
    ext <- file_ext(input$fasta_file3$datapath)
    if (ext == 'gz')
    {input_f <- gunzip(input$fasta_file3$datapath)[1]}
    else
    {input_f <- unzip(input$fasta_file3$datapath)}
    u_num <- ran_number1()
    file_f <- paste0("file_fasta3", u_num, "second") # rename of uploaded assembly file in the genome assembly tab
    destination_directory <- paste0("/opt/shiny-server/samples/sample-apps/assembly_statistics/uploaded_files/",file_f)
    file.copy(input_f, destination_directory)  # copy the renamed uploaded file to a directory in the shiny server folder
    return (destination_directory)

  })
 
  input_fasta_file2 <- eventReactive(input$go2,{
    ext <- file_ext(input$fasta_file2$datapath)
    if (ext == 'gz')
    {input_f <- gunzip(input$fasta_file2$datapath)[1]}
    else 
    {input_f <- unzip(input$fasta_file2$datapath)}
    u_num <- ran_number2()
    file_f <- paste0("file_fasta2", u_num)  # rename of uploaded assembly file in the genome annotation tab
    destination_directory <- paste0("/opt/shiny-server/samples/sample-apps/assembly_statistics/uploaded_files/",file_f)
    file.copy(input_f, destination_directory)  # copy the renamed uploaded file to a directory in the shiny server folder
    return (destination_directory)
    
  })
  
  input_gff_file <- eventReactive(input$go2,{
    
    input_g <- gunzip(input$gff_file$datapath)[1]
    u_num <- ran_number2()
    file_g <- paste0("file_gff", u_num)    # rename of uploaded assembly file in the genome annotation tab
    destination_directory_g <- paste0("/opt/shiny-server/samples/sample-apps/assembly_statistics/uploaded_files/", file_g)
    file.copy(input_g, destination_directory_g) # copy the renamed uploaded file to a directory in the shiny server folder
    return (destination_directory_g)
    
  })
  
  input_transcripts_file <- eventReactive(input$go2,{
    
    input_t <- gunzip(input$transcripts_file$datapath)[1]
    u_num <- ran_number2()
    file_t <- paste0("file_transcript", u_num)  # rename of uploaded transcripts file in the genome annotation tab
    destination_directory_t <- paste0("/opt/shiny-server/samples/sample-apps/assembly_statistics/uploaded_files/", file_t)
    file.copy(input_t, destination_directory_t) # copy the renamed uploaded file to a directory in the shiny server folder
    return (destination_directory_t)
    
  })
  
  
## THIS MODULE OUTPUTS THE REFERENCE NG(X) GRAPH 
  
  Nvalues_plot <- eventReactive(input$go,{
    
    if (is.null(input$references))
    {
      return(NULL)}
    else{
      ref_names <- input$references
      
      # progress indicator for the NG(X) sub-tab of reference genome tab
      progress <- Progress$new(session, min=1, max=15)
      on.exit(progress$close())
      
      progress$set(message = 'Calculation in progress', detail = 'This may take a while...')
      
      for (i in 1:15) {
        progress$set(value = i)
        Sys.sleep(0.5)
      }
      
      # reads pre-computed values
      list1 <- list()
      for(i in 1:length(ref_names)){
        file <- read.table(paste0("/home/nancy/assembly_app/reference_genomes/", ref_names[i]),quote="\"", comment.char="")
        #list1[i] <- log(file)
        list1[i] <- file
      }
      l <- as.vector(list1)
      
      # creates dataframe for the table
      NGX_threshold <- c(1:100)
      X1 <- as.data.frame(l)
      colnames(X1) <- ref_names
      xdata <- NGX_threshold
      NGXdata <- data.frame(xdata, X1)
      #xymelt <- melt(NGXdata, id.vars = "xdata")
      
      # plots the values
     # X <- ggplot(xymelt, aes(x = xdata, y = value, color = variable)) +
       # theme_bw() +
       # geom_line() +
        #ylab("log scaled scaffold_length")+
       # ylab("scaffold_length (bp)") +
       # xlab("NG(X)") +
       # scale_x_continuous(breaks = seq(0, 100, by = 10)) +
       # geom_vline(xintercept = 50)
      
      # adds the figure legend to the plot 
      #p2 <- add_sub(X, "Comparison of genome assembly quality based on NG(X) values.", y = 0, vjust = 0)
      #p5 <- add_sub(p2, "Each curved line charted represents contig lengths of different assemblies at different NG levels")
      return(NGXdata)
      
    }
  })
  
   # outputs the reference NG(X) plot to the tab
    output$NGX_plot <- renderPlot({
      NG <- Nvalues_plot()     
      xymelt <- melt(NG, id.vars = "xdata")

      # plots the values
      X <- ggplot(xymelt, aes(x = xdata, y = value, color = variable)) +
        theme_bw() +
        geom_line() +
        #ylab("log scaled scaffold_length")+
        ylab("scaffold_length (bp)") +
        xlab("NG(X)") +
        scale_x_continuous(breaks = seq(0, 100, by = 10)) +
        geom_vline(xintercept = 50)

      # adds the figure legend to the plot 
      p2 <- add_sub(X, "Comparison of genome assembly quality based on NG(X) values.", y = 0, vjust = 0)
      p5 <- add_sub(p2, "Each curved line charted represents contig lengths of different assemblies at different NG levels")
    ggdraw(p5)
  })
    
   # Download the plot as png file
    # output$assemblyplot <- downloadHandler(
    #          NG <- Nvalues_plot()
    #  xymelt <- melt(NG, id.vars = "xdata")

      # plots the values
    #  X <- ggplot(xymelt, aes(x = xdata, y = value, color = variable)) +
    #    theme_bw() +
    #    geom_line() +
        #ylab("log scaled scaffold_length")+
     #   ylab("scaffold_length (bp)") +
      #  xlab("NG(X)") +
      #  scale_x_continuous(breaks = seq(0, 100, by = 10)) +
      #  geom_vline(xintercept = 50)

      # adds the figure legend to the plot 
     # p2 <- add_sub(X, "Comparison of genome assembly quality based on NG(X) values.", y = 0, vjust = 0)
     # p5 <- add_sub(p2, "Each curved line charted represents contig lengths of different assemblies at different NG levels")
  
      #   filename = function() {

      # "Reference_NG(X)_plot.png"
    # },
     #content = function(file) {
      # ggsave(file, plot = p5, width = 16, height = 10.4)

    # })
 
   # Download the table as csv file
     output$assemblyplot <- downloadHandler(

     filename = function() {

       "Reference_NG(X).csv"
     },
     content = function(file) {
         write.csv(Nvalues_plot(), file)
         #ggsave(file, plot = Nvalues_plot(), width = 16, height = 10.4)

     })

 
  
## THIS MODULE OUTPUTS THE NG(X) GRAPH FOR THE USER UPLOADED GENOME ASSEMBLY AND COMPARE THOSE WITH THE PRE-COMPUTED REFRENCE VALUES
     
  Nvalues_plot1 <- eventReactive(input$go1,{
    
    # If the reference genomes are selected by the user
    if (!is.null(input$fasta_file1$datapath)) {
      if(!is.null(input$references1)) {
        #if (!is.null(input$fasta_file3$datapath)){
        estimated_genome_size <- input$estimated_genome_size
        ref_names <- input$references1
        example_file.fasta <- input_fasta_file1()
        #example_file2.fasta <- input_fasta_file3()
        u_num <- ran_number1()
        output_NG <- paste0("file_one", u_num,"_NG.txt")
        #output2_NG <- paste0("file_second", u_num,"_NG.txt")
        NGX_threshold <- c(1:100)
        user_name <- input$user_input1
        if (!is.null(input$fasta_file3$datapath)){
         example_file2.fasta <- input_fasta_file3()
         output2_NG <- paste0("file_second", u_num,"_NG.txt")
         user_name3 <- input$user_input3

        # progress indicator
        progress <- Progress$new(session, min=1, max=15)
        on.exit(progress$close())
        
        progress$set(message = 'Calculation in progress', detail = 'This may take a while...')
        
        for (i in 1:15) {
          progress$set(value = i)
          Sys.sleep(0.5)
        }
        
        # This command calls the python script which takes three input arguments: input file name, output file name and the estimated genome size and calculate NG(X) values 
        system(paste("/home/nancy/assembly_app/python3/Python-3.6.4/python /home/nancy/assembly_app/NG.py", example_file.fasta, output_NG, estimated_genome_size), intern = TRUE)
        output <- read.table(output_NG, quote="\"", comment.char="")
        NGX_test <- output$V1
        NA_vector <- rep(NA, 100-(length(NGX_test)))
        NGX_scaffoldlength <- c(NGX_test, NA_vector)  
        system(paste("/home/nancy/assembly_app/python3/Python-3.6.4/python /home/nancy/assembly_app/NG.py", example_file2.fasta, output2_NG, estimated_genome_size), intern = TRUE)
        output2 <- read.table(output2_NG, quote="\"", comment.char="")
        NGX_test2 <- output2$V1
        NA_vector2 <- rep(NA, 100-(length(NGX_test2)))
        NGX_scaffoldlength2 <- c(NGX_test2, NA_vector2)        

        # reads the reference genome values
        list1 <- list()
        for(i in 1:length(ref_names)){
          file <- read.table(paste0("/home/nancy/assembly_app/reference_genomes/", ref_names[i]),quote="\"", comment.char="")
          list1[i] <- file
           #list1[i] <- log(file)
        }
        l <- as.vector(list1)
        
        # creates dataframe for the plot
        user_n <- NGX_scaffoldlength
        user_n2 <- NGX_scaffoldlength2
        X1 <- as.data.frame(l)
        colnames(X1) <- ref_names
        comb <- cbind(X1, user_n, user_n2)
        colnames(comb) <- c(ref_names, user_name, user_name3)
        xdata <- NGX_threshold
        NGXdata <- data.frame(xdata, comb)
        return(NGXdata)
      }
       else 
       {
        system(paste("/home/nancy/assembly_app/python3/Python-3.6.4/python /home/nancy/assembly_app/NG.py", example_file.fasta, output_NG, estimated_genome_size), intern = TRUE)
        output <- read.table(output_NG, quote="\"", comment.char="")
        NGX_test <- output$V1
        NA_vector <- rep(NA, 100-(length(NGX_test)))
        NGX_scaffoldlength <- c(NGX_test, NA_vector)
        # reads the reference genome values
        list1 <- list()
        for(i in 1:length(ref_names)){
          file <- read.table(paste0("/home/nancy/assembly_app/reference_genomes/", ref_names[i]),quote="\"", comment.char="")
          list1[i] <- file
           #list1[i] <- log(file)
        }
        l <- as.vector(list1)

        # creates dataframe for the plot
        user_n <- NGX_scaffoldlength
        X1 <- as.data.frame(l)
        colnames(X1) <- ref_names
        comb <- cbind(X1, user_n)
        colnames(comb) <- c(ref_names, user_name)
        xdata <- NGX_threshold
        NGXdata <- data.frame(xdata, comb)
        return(NGXdata)
       }
      }
     # If the user does not select any reference genome 
      else{
        example_file.fasta <- input_fasta_file1()
        #example_file2.fasta <- input_fasta_file3()
        estimated_genome_size <- input$estimated_genome_size
        u_num <- ran_number1()
        output_NG <- paste0("file_first", u_num,"_NG.txt")
        #output2_NG <- paste0("file_second", u_num,"_NG.txt")
        NGX_threshold <- c(1:100)
        user_name <- input$user_input1
        #user_name3 <- input$user_input3
        
        # progress indicator
        progress <- Progress$new(session, min=1, max=15)
        on.exit(progress$close())

        progress$set(message = 'Calculation in progress', detail = 'This may take a while...')

        for (i in 1:15) {
          progress$set(value = i)
          Sys.sleep(0.5)
        }

        #if the user uploads the second assembly fasta file
        if(!is.null(input$fasta_file3$datapath))
        {
        example_file2.fasta <- input_fasta_file3()
        output2_NG <- paste0("file_second", u_num,"_NG.txt")
        user_name3 <- input$user_input3
 
  
        # progress indicator
        #progress <- Progress$new(session, min=1, max=15)
        #on.exit(progress$close())
        
        #progress$set(message = 'Calculation in progress', detail = 'This may take a while...')
        
        #for (i in 1:15) {
        #  progress$set(value = i)
        #  Sys.sleep(0.5)
        #}
        
        # This command calls the python script which takes three input arguments: input file name, output file name and the estimated genome size and calculate NG(X) values 
        system(paste("/home/nancy/assembly_app/python3/Python-3.6.4/python /home/nancy/assembly_app/NG.py", example_file.fasta, output_NG, estimated_genome_size), intern = TRUE)
        output <- read.table(output_NG, quote="\"", comment.char="")
        NGX_test <- output$V1
        NA_vector <- rep(NA, 100-(length(NGX_test)))
        NGX_scaffoldlength <- c(NGX_test, NA_vector)
        
        system(paste("/home/nancy/assembly_app/python3/Python-3.6.4/python /home/nancy/assembly_app/NG.py", example_file2.fasta, output2_NG, estimated_genome_size), intern = TRUE)
        output2 <- read.table(output2_NG, quote="\"", comment.char="")
        NGX_test2 <- output2$V1
        NA_vector2 <- rep(NA, 100-(length(NGX_test2)))
        NGX_scaffoldlength2 <- c(NGX_test2, NA_vector2)        

        
        #user_n <- log(NGX_scaffoldlength)
        user_n <- NGX_scaffoldlength
        user_n2 <- NGX_scaffoldlength2
        comb <- cbind(user_n, user_n2)
        colnames(comb) <- c(user_name, user_name3)
        xdata <- NGX_threshold
        NGXdata <- data.frame(xdata, comb)
        
        return(NGXdata)
      }
      else {
        system(paste("/home/nancy/assembly_app/python3/Python-3.6.4/python /home/nancy/assembly_app/NG.py", example_file.fasta, output_NG, estimated_genome_size), intern = TRUE)
        output <- read.table(output_NG, quote="\"", comment.char="")
        NGX_test <- output$V1
        NA_vector <- rep(NA, 100-(length(NGX_test)))
        NGX_scaffoldlength <- c(NGX_test, NA_vector)
        user_n <- NGX_scaffoldlength
        comb <- cbind(user_n)
        colnames(comb) <- c(user_name)
        xdata <- NGX_threshold
        NGXdata <- data.frame(xdata, comb)

        return(NGXdata)   
           }
      }
    }
  })
  
  
  # outputs the plot to the tab
  output$NGX_plot1 <- renderPlot({
        NG <- Nvalues_plot1()  
         xymelt <- melt(NG, id.vars = "xdata")

        # plots the values
        X <- ggplot(xymelt, aes(x = xdata, y = value, color = variable)) +
          theme_bw() +
          geom_line() +
          ylab("scaffold_length (bp)")+
          xlab("NG(X)") +
          scale_x_continuous(breaks = seq(0, 100, by = 10)) +
          geom_vline(xintercept = 50)

        # adds the figure legend to the plot
        p2 <- add_sub(X, "Comparison of genome assembly quality based on NG(X) values.", y = 0, vjust = 0)
        p5 <- add_sub(p2, "Each curved line charted represents contig lengths of different assemblies at different NG levels")

    ggdraw(p5)
  })
  

     output$assemblyplot1 <- downloadHandler(

     filename = function() {
         u_num <- ran_number1()
         file_name = paste0("NG(X)_values",u_num,".csv")
        }, 
     content = function(file) {
         write.csv(Nvalues_plot1(), file)

     })


  # Download the plot as png file
     #output$assemblyplot1 <- downloadHandler(
     
     #filename = function() {
     #  u_num <- ran_number1()
     #  file_name = paste0("NG(X)_plot",u_num,".png") 
  
    # },
    # content = function(file) {
     #  ggsave(file, plot = Nvalues_plot1(), width = 16, height = 10.4)
       
     #})
  
  
  
## THIS MODULE OUTPUTS THE REFERENCE GENOME ASSEMBLY METRICS TABLE 

  
  data_f <- eventReactive(input$go,{
    
    if (is.null(input$references))
    {
      return(NULL)}
    else{
      ref_names <- input$references
      
      # progress indicator
      progress <- Progress$new(session, min=1, max=15)
      on.exit(progress$close())
      progress$set(message = 'Calculation in progress', detail = 'This may take a while...')
      
      for (i in 1:15) {
        progress$set(value = i)
        Sys.sleep(0.5)
      }
      
      # reads the reference genome values
      list1 = list()
      for(i in 1:length(ref_names)){
        file <- read.table(paste0("/home/nancy/assembly_app/reference_genomes/", ref_names[i], "_assemblymetrics"),quote="\"", comment.char="")
        list1[i] <- file
      }
      l <- as.vector(list1)
      
      # creates dataframe for the table
      X1 <- as.data.frame(l)
      colnames(X1) <- ref_names
      N_metrics <- X1
      row.names(N_metrics) <- c ("Number of scaffolds", "Total size of scaffolds", "Total scaffold length as percentage of assumed genome size", "useful amount of scaffold sequences (>=25K nt)", "% of estimated genome that is useful", "Longest scaffold", "Shortest scaffold", "Number of scaffolds > 1K nt", "Number of scaffolds > 10K nt", "Number of scaffolds > 100K nt", "Number of scaffolds > 1M nt", "Number of scaffolds > 10M nt", "N50", "L50", "NG50", "LG50", "%A", "%C", "%G", "%T", "Total Number of Ns", "%N")
      return(N_metrics)
    }
    
  })
  
  # outputs the reference assembly metrics table to the tab 
  output$lengthmetrics = DT::renderDataTable({
    
    ref_names <- input$references
    cn <- c(ref_names)
    data_fr <- data_f()
    if (length(data_fr)==0)
      return(NULL)
    
    datatable(data_fr, selection = 'single', rownames = c ("Number of scaffolds", "Total size of scaffolds (bp)", "Total scaffold length as percentage of assumed genome size (%)", "useful amount of scaffold sequences (>=25K nt) (bp)", "% of estimated genome that is useful", "Longest scaffold (bp)", "Shortest scaffold(bp)", "Number of scaffolds > 1K nt", "Number of scaffolds > 10K nt", "Number of scaffolds > 100K nt", "Number of scaffolds > 1M nt", "Number of scaffolds > 10M nt", "N50 (bp)", "L50", "NG50 (bp)", "LG50", "%A", "%C", "%G", "%T", "Total Number of Ns", "%N"), colnames = cn)
  })
  
  
  # generates pop-up plots on selecting the rows of reference assembly metrics table
  observeEvent(input$lengthmetrics_rows_selected,
               
               {
                 
                 data_fr <- data_f()
                 ref_names <- input$references
                 len <- length(ref_names) 
                 cn <- c(ref_names)
                 new_labels <- c()
                 for (i in cn) {
                   first_word <- unlist(strsplit(i, "_"))
                   new_labels <- c(new_labels,first_word[1])
                     } 
                 if (input$lengthmetrics_rows_selected == 13 )
                 {
                   showModal(modalDialog(
                     title = "An N50 contig size of N means that 50% of the assembled bases are contained in contigs of length N or larger. N50 sizes are often used as a measure of assembly quality because they capture how much of the genome is covered by relatively large contigs.",
                     
                     
                     outputplot1 <- renderPlot({
                     plt <- barplot(as.numeric(data_fr[13,]), main="N50", ylab="N50 value (bp)", col=rainbow(len))
                     text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                     })
                     
                     
                   ))
                 }
                 if (input$lengthmetrics_rows_selected == 1 )
                 {
                   showModal(modalDialog(
                     title = "Total number of sequences in the genome assembly",                    
                     outputplot1 <- renderPlot({
                       g <- as.numeric(data_fr[1,])
                       plt <- barplot(g , ylab="Number of sequences", col=rainbow(len))
                       text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                     })
                     
                     
                   ))
                 }
                 if (input$lengthmetrics_rows_selected == 2 )
                 {
                   showModal(modalDialog(
                     title = "Total size of the scaffolds",
                     
                     
                     outputplot1 <- renderPlot({
                     plt <- barplot(as.numeric(data_fr[2,]) , ylab="scaffold length (bp)", col=rainbow(len))
                     text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                     })
                   ))
                 }
                 
                 if (input$lengthmetrics_rows_selected == 3 )
                 {
                   showModal(modalDialog(
                     title = "Total scaffold length as percentage of assumed genome size (%)",
                     
                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr[3,]), ylab="scaffold length (bp)", col=rainbow(len))
                       text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                     })
                     
                     
                   ))
                 }
                 if (input$lengthmetrics_rows_selected == 4)
                 {
                   showModal(modalDialog(
                     title = "useful amount of scaffold sequences (>=25K nt) (bp)",
                     
                     outputplot1 <- renderPlot({
                     plt <- barplot(as.numeric(data_fr[4,]), ylab="scaffold length (bp)", col=rainbow(len))
                     text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                    })
                     
                   ))
                 }
                 
                 if (input$lengthmetrics_rows_selected == 5)
                 {
                   showModal(modalDialog(
                     title =  "% of estimated genome that is useful",
                     
                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr[5,]),ylab="useful estimated genome size (%)", col=rainbow(len))
                       text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       
                       })                    
                     
                     
                     
                   ))
                 }
                 
                 if (input$lengthmetrics_rows_selected == 6)
                 {
                   showModal(modalDialog(
                     title = "Longest scaffold (bp)",
                     
                     outputplot1 <- renderPlot({
                     plt <- barplot(as.numeric(data_fr[6,]) , ylab="scaffold length (bp)", col=rainbow(len))
                     text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                     })
                   ))
                 }
                 
                 if (input$lengthmetrics_rows_selected == 7)
                 {
                   showModal(modalDialog(
                     title = "Shortest scaffold (bp)",
                     
                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr[7,]) , ylab="scaffold length (bp)", col=rainbow(len))
                       text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$lengthmetrics_rows_selected == 8)
                 {
                   showModal(modalDialog(
                     title = "Number of scaffolds > 1K nt",
                     
                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr[8,]) ,ylab="Number of scaffolds", col=rainbow(len))
                       text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$lengthmetrics_rows_selected == 9)
                 {
                   showModal(modalDialog(
                     title = "Number of scaffolds > 10K nt",
                     
                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr[9,]) ,ylab="Number of scaffolds", col=rainbow(len))
                       text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 
                 if (input$lengthmetrics_rows_selected == 10)
                 {
                   showModal(modalDialog(
                     title = "Number of scaffolds > 100K nt",
                     
                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr[10,]) , ylab="Number of scaffolds", col=rainbow(len))
                       text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$lengthmetrics_rows_selected == 11)
                 {
                   showModal(modalDialog(
                     title =  "Number of scaffolds > 1M nt",
                     
                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr[11,]) , ylab="Number of scaffolds", col=rainbow(len))
                       text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$lengthmetrics_rows_selected == 12)
                 {
                   showModal(modalDialog(
                     title = "Number of scaffolds > 10M nt",
                     
                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr[12,]) , ylab="Number of scaffolds", col=rainbow(len))
                       text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 
                 if (input$lengthmetrics_rows_selected == 14)
                 {
                   showModal(modalDialog(
                     title = "A scaffold/contig L50 is calculated as the number of sequences whose sum of lengths makes up 50% or more of the total assembly length.",
                     
                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr[14,]) ,ylab="Number of scaffolds", col=rainbow(len))
                       text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$lengthmetrics_rows_selected == 15)
                 {
                   showModal(modalDialog(
                     title = "NG50 is calculated in the same manner as N50 but uses estimated genome size instead of total assembly length, is more useful as it reduces the bias caused by assembly lengths.",
                     
                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr[15,]) , ylab="scaffold length (bp)", col=rainbow(len))
                       text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$lengthmetrics_rows_selected == 16)
                 {
                   showModal(modalDialog(
                     title = "A scaffold/contig LG50 is calculated as the number of sequences whose sum of lengths makes up 50% or more of the estimated genome size.",
                     
                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr[16,]) , ylab="Number of scaffolds", col=rainbow(len))
                       text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$lengthmetrics_rows_selected == 17)
                 {
                   showModal(modalDialog(
                     title = "Percentage A",
                     
                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr[17,]) , ylab="%A", col=rainbow(len))
                       text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$lengthmetrics_rows_selected == 18)
                 {
                   showModal(modalDialog(
                     title = "Percentage C",
                     
                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr[18,]) ,ylab="%C", col=rainbow(len))
                       text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$lengthmetrics_rows_selected == 19)
                 {
                   showModal(modalDialog(
                     title = "Percentage G",
                     
                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr[19,]) ,ylab="%G", col=rainbow(len))
                       text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$lengthmetrics_rows_selected == 20)
                 {
                   showModal(modalDialog(
                     title = "Percentage T",
                     
                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr[20,]) ,ylab="%T", col=rainbow(len))
                       text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$lengthmetrics_rows_selected == 21)
                 {
                   showModal(modalDialog(
                     title = "Number of Ns: gives an estimate of the gaps in the assembled genome sequence",
                     
                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr[21,]) ,ylab="Number of Ns", col=rainbow(len))
                       text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$lengthmetrics_rows_selected == 22)
                 {
                   showModal(modalDialog(
                     title = "Percentage N",
                     
                     outputplot1 <- renderPlot({
                     plt<- barplot(as.numeric(data_fr[22,]) , ylab="%N", col=rainbow(len))
                     text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                     })
                   ))
                 } 
               })
  
  # Download the reference assembly metrics table in csv format
     output$downloadassemblytable <- downloadHandler(
     
     filename = function() {
       
       "Reference_Assembly_metrics_table.csv"
     },
     content = function(file) {
       
       write.csv(data_f(), file)
    })  

     
## THIS MODULE OUTPUTS THE ASSEMBLY METRICS TABLE FOR THE USER UPLOADED GENOME ASSEMBLY AND COMPARE THE VALUES WITH PRE-COMPUTED REFERENCE VALUES

  data_f1 <- eventReactive(input$go1,{
    
    estimated_genome_size <- input$estimated_genome_size
    
    # progress indicator
    progress <- Progress$new(session, min=1, max=15)
    on.exit(progress$close())
    
    progress$set(message = 'Calculation in progress', detail = 'This may take a while...')
    
    for (i in 1:15) {
      progress$set(value = i)
      Sys.sleep(0.5)
    }
    
    # If the reference genomes are selected by the user
    if (!is.null(input$fasta_file1$datapath)) {
      if(!is.null(input$references1)) {
        ref_names <- input$references1
        input_name <- input$user_input1
        #input_name_second <- input$user_input3
        example_file.fasta <- input_fasta_file1()
        #example_file_second.fasta <- input_fasta_file3()
        u_num <- ran_number1()
        output_assembly_metrics <- paste0("file", u_num,"_assembly_m.txt")
        #output_assembly_metrics2 <- paste0("file_second", u_num,"_assembly_m.txt")

        #if the second genome fasta file is uploaded
         if (!is.null(input$fasta_file3$datapath)){
         example_file_second.fasta <- input_fasta_file3()
         output_assembly_metrics2 <- paste0("file_second", u_num,"_assembly_m.txt")
         input_name_second <- input$user_input3
        # This command calls the python script which takes three input arguments: input file name, output file name and the estimated genome size and calculate all the assembly metrics
        system(paste("/home/nancy/assembly_app/python3/Python-3.6.4/python /home/nancy/assembly_app/assembly_stats.py", example_file.fasta, output_assembly_metrics, estimated_genome_size), intern = TRUE)
        output <- read.table(output_assembly_metrics, quote="\"", comment.char="")
        metrics_fasta <- output$V1
        num_scaffolds <- metrics_fasta[1]
        total_size_scaffolds <- metrics_fasta[2]
        total_scaffold_length_as_precentage_genome_size <- metrics_fasta[3]
        useful_amount_scaffold_sequence <- metrics_fasta[4]
        percentage_estimated_genome_useful <- metrics_fasta[5]
        longest_scaffold <- metrics_fasta[6]
        shortest_scaffold <- metrics_fasta[7]
        scaffolds_greater_than_1k <- metrics_fasta[8]
        percentage_scaffolds_greater_than_1k <- metrics_fasta[9]
        scaffolds_greater_than_10k <- metrics_fasta[10]
        percentage_scaffolds_greater_than_10k <- metrics_fasta[11]
        scaffolds_greater_than_100k <- metrics_fasta[12]
        percentage_scaffolds_greater_than_100k <- metrics_fasta[13]
        scaffolds_greater_than_1M <- metrics_fasta[14]
        percentage_scaffolds_greater_than_1M <- metrics_fasta[15]
        scaffolds_greater_than_10M <- metrics_fasta[16]
        percentage_scaffolds_greater_than_10M <- metrics_fasta[17]
        N50 <- metrics_fasta[18]
        L50 <- metrics_fasta[19]
        NG50 <- metrics_fasta[20]
        LG50 <- metrics_fasta[21]
        percentage_A <- metrics_fasta[22]
        percentage_C <- metrics_fasta[23]
        percentage_G <- metrics_fasta[24]
        percentage_T <- metrics_fasta[25]
        Total_N <- metrics_fasta[26]
        percentage_N <- metrics_fasta[27]
        
        system(paste("/home/nancy/assembly_app/python3/Python-3.6.4/python /home/nancy/assembly_app/assembly_stats.py", example_file_second.fasta, output_assembly_metrics2, estimated_genome_size), intern = TRUE)
        output_second <- read.table(output_assembly_metrics2, quote="\"", comment.char="")
        metrics_fasta2 <- output_second$V1
        num_scaffolds2 <- metrics_fasta2[1]
        total_size_scaffolds2 <- metrics_fasta2[2]
        total_scaffold_length_as_precentage_genome_size2 <- metrics_fasta2[3]
        useful_amount_scaffold_sequence2 <- metrics_fasta2[4]
        percentage_estimated_genome_useful2 <- metrics_fasta2[5]
        longest_scaffold2 <- metrics_fasta2[6]
        shortest_scaffold2 <- metrics_fasta2[7]
        scaffolds_greater_than_1k2 <- metrics_fasta2[8]
        percentage_scaffolds_greater_than_1k2 <- metrics_fasta2[9]
        scaffolds_greater_than_10k2 <- metrics_fasta2[10]
        percentage_scaffolds_greater_than_10k2 <- metrics_fasta2[11]
        scaffolds_greater_than_100k2 <- metrics_fasta2[12]
        percentage_scaffolds_greater_than_100k2 <- metrics_fasta2[13]
        scaffolds_greater_than_1M2 <- metrics_fasta2[14]
        percentage_scaffolds_greater_than_1M2 <- metrics_fasta2[15]
        scaffolds_greater_than_10M2 <- metrics_fasta2[16]
        percentage_scaffolds_greater_than_10M2 <- metrics_fasta2[17]
        N502 <- metrics_fasta2[18]
        L502 <- metrics_fasta2[19]
        NG502 <- metrics_fasta2[20]
        LG502 <- metrics_fasta2[21]
        percentage_A2 <- metrics_fasta2[22]
        percentage_C2 <- metrics_fasta2[23]
        percentage_G2 <- metrics_fasta2[24]
        percentage_T2 <- metrics_fasta2[25]
        Total_N2 <- metrics_fasta2[26]
        percentage_N2 <- metrics_fasta2[27]

        # reads the reference genome values
        list1 = list()
        for(i in 1:length(ref_names)){
          file <- read.table(paste0("/home/nancy/assembly_app/reference_genomes/", ref_names[i], "_assemblymetrics"),quote="\"", comment.char="")
          list1[i] <- file
        }
        l <- as.vector(list1)
        
        # creates dataframe for the table
        X1 <- as.data.frame(l)
        colnames(X1) <- ref_names
       
        user_input <- c(num_scaffolds, total_size_scaffolds, total_scaffold_length_as_precentage_genome_size, useful_amount_scaffold_sequence, percentage_estimated_genome_useful, longest_scaffold, shortest_scaffold, scaffolds_greater_than_1k, scaffolds_greater_than_10k, scaffolds_greater_than_100k, scaffolds_greater_than_1M, scaffolds_greater_than_10M, N50, L50, NG50, LG50, percentage_A, percentage_C, percentage_G, percentage_T, Total_N, percentage_N)
            user_input_second <- c(num_scaffolds2, total_size_scaffolds2, total_scaffold_length_as_precentage_genome_size2, useful_amount_scaffold_sequence2, percentage_estimated_genome_useful2, longest_scaffold2, shortest_scaffold2, scaffolds_greater_than_1k2, scaffolds_greater_than_10k2, scaffolds_greater_than_100k2, scaffolds_greater_than_1M2, scaffolds_greater_than_10M2, N502, L502, NG502, LG502, percentage_A2, percentage_C2, percentage_G2, percentage_T2, Total_N2, percentage_N2)
        N_metrics <- cbind(X1, user_input, user_input_second)
        row.names(N_metrics) <- c ("Number of scaffolds", "Total size of scaffolds", "Total scaffold length as percentage of assumed genome size", "useful amount of scaffold sequences (>=25K nt)", "% of estimated genome that is useful", "Longest scaffold", "Shortest scaffold", "Number of scaffolds > 1K nt", "Number of scaffolds > 10K nt", "Number of scaffolds > 100K nt", "Number of scaffolds > 1M nt", "Number of scaffolds > 10M nt", "N50", "L50", "NG50", "LG50", "%A", "%C", "%G", "%T", "Total Number of Ns","%N")
       colnames(N_metrics) <- c(ref_names,input_name, input_name_second) 
       return(N_metrics)
    }
     else 
     {
        system(paste("/home/nancy/assembly_app/python3/Python-3.6.4/python /home/nancy/assembly_app/assembly_stats.py", example_file.fasta, output_assembly_metrics, estimated_genome_size), intern = TRUE)
        output <- read.table(output_assembly_metrics, quote="\"", comment.char="")
        metrics_fasta <- output$V1
        num_scaffolds <- metrics_fasta[1]
        total_size_scaffolds <- metrics_fasta[2]
        total_scaffold_length_as_precentage_genome_size <- metrics_fasta[3]
        useful_amount_scaffold_sequence <- metrics_fasta[4]
        percentage_estimated_genome_useful <- metrics_fasta[5]
        longest_scaffold <- metrics_fasta[6]
        shortest_scaffold <- metrics_fasta[7]
        scaffolds_greater_than_1k <- metrics_fasta[8]
        percentage_scaffolds_greater_than_1k <- metrics_fasta[9]
        scaffolds_greater_than_10k <- metrics_fasta[10]
        percentage_scaffolds_greater_than_10k <- metrics_fasta[11]
        scaffolds_greater_than_100k <- metrics_fasta[12]
        percentage_scaffolds_greater_than_100k <- metrics_fasta[13]
        scaffolds_greater_than_1M <- metrics_fasta[14]
        percentage_scaffolds_greater_than_1M <- metrics_fasta[15]
        scaffolds_greater_than_10M <- metrics_fasta[16]
        percentage_scaffolds_greater_than_10M <- metrics_fasta[17]
        N50 <- metrics_fasta[18]
        L50 <- metrics_fasta[19]
        NG50 <- metrics_fasta[20]
        LG50 <- metrics_fasta[21]
        percentage_A <- metrics_fasta[22]
        percentage_C <- metrics_fasta[23]
        percentage_G <- metrics_fasta[24]
        percentage_T <- metrics_fasta[25]
        Total_N <- metrics_fasta[26]
        percentage_N <- metrics_fasta[27]
                # reads the reference genome values
        list1 = list()
        for(i in 1:length(ref_names)){
          file <- read.table(paste0("/home/nancy/assembly_app/reference_genomes/", ref_names[i], "_assemblymetrics"),quote="\"", comment.char="")
          list1[i] <- file
        }
        l <- as.vector(list1)

        # creates dataframe for the table
        X1 <- as.data.frame(l)
        colnames(X1) <- ref_names

        user_input <- c(num_scaffolds, total_size_scaffolds, total_scaffold_length_as_precentage_genome_size, useful_amount_scaffold_sequence, percentage_estimated_genome_useful, longest_scaffold, shortest_scaffold, scaffolds_greater_than_1k, scaffolds_greater_than_10k, scaffolds_greater_than_100k, scaffolds_greater_than_1M, scaffolds_greater_than_10M, N50, L50, NG50, LG50, percentage_A, percentage_C, percentage_G, percentage_T, Total_N, percentage_N)
        N_metrics <- cbind(X1, user_input)
        row.names(N_metrics) <- c ("Number of scaffolds", "Total size of scaffolds", "Total scaffold length as percentage of assumed genome size", "useful amount of scaffold sequences (>=25K nt)", "% of estimated genome that is useful", "Longest scaffold", "Shortest scaffold", "Number of scaffolds > 1K nt", "Number of scaffolds > 10K nt", "Number of scaffolds > 100K nt", "Number of scaffolds > 1M nt", "Number of scaffolds > 10M nt", "N50", "L50", "NG50", "LG50", "%A", "%C", "%G", "%T", "Total Number of Ns","%N")
       colnames(N_metrics) <- c(ref_names,input_name)
       return(N_metrics)
     }

         
   }  
   #}
      
      # If the user does not select any reference genome 
      else{
        estimated_genome_size <- input$estimated_genome_size
        u_num <- ran_number1()
        #input_name_second <- input$user_input3
        #example_file_second.fasta <- input_fasta_file3()
        #output_assembly_metrics2 <- paste0("file_second", u_num,"_assembly_m.txt")
        example_file.fasta <- input_fasta_file1()
        input_name <- input$user_input1
        output_assembly_metrics <- paste0("file", u_num,"_assembly_m.txt")
      
        #if the user uploads the second genome fasta file
        if(!is.null(input$fasta_file3$datapath))
        {
        example_file_second.fasta <- input_fasta_file3()
        output_assembly_metrics2 <- paste0("file_second", u_num,"_assembly_m.txt")
        input_name_second <- input$user_input3

             
        # This command calls the python script which takes three input arguments: input file name, output file name and the estimated genome size and calculate all the assembly metrics
        system(paste("/home/nancy/assembly_app/python3/Python-3.6.4/python /home/nancy/assembly_app/assembly_stats.py", example_file.fasta, output_assembly_metrics, estimated_genome_size), intern = TRUE)
        output <- read.table(output_assembly_metrics, quote="\"", comment.char="")
        metrics_fasta <- output$V1
        num_scaffolds <- metrics_fasta[1]
        total_size_scaffolds <- metrics_fasta[2]
        total_scaffold_length_as_precentage_genome_size <- metrics_fasta[3]
        useful_amount_scaffold_sequence <- metrics_fasta[4]
        percentage_estimated_genome_useful <- metrics_fasta[5]
        longest_scaffold <- metrics_fasta[6]
        shortest_scaffold <- metrics_fasta[7]
        scaffolds_greater_than_1k <- metrics_fasta[8]
        percentage_scaffolds_greater_than_1k <- metrics_fasta[9]
        scaffolds_greater_than_10k <- metrics_fasta[10]
        percentage_scaffolds_greater_than_10k <- metrics_fasta[11]
        scaffolds_greater_than_100k <- metrics_fasta[12]
        percentage_scaffolds_greater_than_100k <- metrics_fasta[13]
        scaffolds_greater_than_1M <- metrics_fasta[14]
        percentage_scaffolds_greater_than_1M <- metrics_fasta[15]
        scaffolds_greater_than_10M <- metrics_fasta[16]
        percentage_scaffolds_greater_than_10M <- metrics_fasta[17]
        N50 <- metrics_fasta[18]
        L50 <- metrics_fasta[19]
        NG50 <- metrics_fasta[20]
        LG50 <- metrics_fasta[21]
        percentage_A <- metrics_fasta[22]
        percentage_C <- metrics_fasta[23]
        percentage_G <- metrics_fasta[24]
        percentage_T <- metrics_fasta[25]
        Total_N <- metrics_fasta[26]
        percentage_N <- metrics_fasta[27]
        
        system(paste("/home/nancy/assembly_app/python3/Python-3.6.4/python /home/nancy/assembly_app/assembly_stats.py", example_file_second.fasta, output_assembly_metrics2, estimated_genome_size), intern = TRUE)
        output_second <- read.table(output_assembly_metrics2, quote="\"", comment.char="")
        metrics_fasta2 <- output_second$V1
        num_scaffolds2 <- metrics_fasta2[1]
        total_size_scaffolds2 <- metrics_fasta2[2]
        total_scaffold_length_as_precentage_genome_size2 <- metrics_fasta2[3]
        useful_amount_scaffold_sequence2 <- metrics_fasta2[4]
        percentage_estimated_genome_useful2 <- metrics_fasta2[5]
        longest_scaffold2 <- metrics_fasta2[6]
        shortest_scaffold2 <- metrics_fasta2[7]
        scaffolds_greater_than_1k2 <- metrics_fasta2[8]
        percentage_scaffolds_greater_than_1k2 <- metrics_fasta2[9]
        scaffolds_greater_than_10k2 <- metrics_fasta2[10]
        percentage_scaffolds_greater_than_10k2 <- metrics_fasta2[11]
        scaffolds_greater_than_100k2 <- metrics_fasta2[12]
        percentage_scaffolds_greater_than_100k2 <- metrics_fasta2[13]
        scaffolds_greater_than_1M2 <- metrics_fasta2[14]
        percentage_scaffolds_greater_than_1M2 <- metrics_fasta2[15]
        scaffolds_greater_than_10M2 <- metrics_fasta2[16]
        percentage_scaffolds_greater_than_10M2 <- metrics_fasta2[17]
        N502 <- metrics_fasta2[18]
        L502 <- metrics_fasta2[19]
        NG502 <- metrics_fasta2[20]
        LG502 <- metrics_fasta2[21]
        percentage_A2 <- metrics_fasta2[22]
        percentage_C2 <- metrics_fasta2[23]
        percentage_G2 <- metrics_fasta2[24]
        percentage_T2 <- metrics_fasta2[25]
        Total_N2 <- metrics_fasta2[26]
        percentage_N2 <- metrics_fasta2[27]


        # creating dataframe for the table
        user_input <- c(num_scaffolds, total_size_scaffolds, total_scaffold_length_as_precentage_genome_size, useful_amount_scaffold_sequence, percentage_estimated_genome_useful, longest_scaffold, shortest_scaffold, scaffolds_greater_than_1k, scaffolds_greater_than_10k, scaffolds_greater_than_100k, scaffolds_greater_than_1M, scaffolds_greater_than_10M, N50, L50, NG50, LG50, percentage_A, percentage_C, percentage_G, percentage_T, Total_N, percentage_N)
          user_input_second <- c(num_scaffolds2, total_size_scaffolds2, total_scaffold_length_as_precentage_genome_size2, useful_amount_scaffold_sequence2, percentage_estimated_genome_useful2, longest_scaffold2, shortest_scaffold2, scaffolds_greater_than_1k2, scaffolds_greater_than_10k2, scaffolds_greater_than_100k2, scaffolds_greater_than_1M2, scaffolds_greater_than_10M2, N502, L502, NG502, LG502, percentage_A2, percentage_C2, percentage_G2, percentage_T2, Total_N2, percentage_N2)
        N_metrics <- cbind(user_input, user_input_second)
        #N_metrics <- as.matrix(user_input)
        row.names(N_metrics) <- c ("Number of scaffolds", "Total size of scaffolds", "Total scaffold length as percentage of assumed genome size", "useful amount of scaffold sequences (>=25K nt)", "% of estimated genome that is useful", "Longest scaffold", "Shortest scaffold", "Number of scaffolds > 1K nt", "Number of scaffolds > 10K nt", "Number of scaffolds > 100K nt", "Number of scaffolds > 1M nt", "Number of scaffolds > 10M nt", "N50", "L50", "NG50", "LG50", "%A", "%C", "%G", "%T", "Total Number of Ns","%N")
        colnames(N_metrics) <- c(input_name, input_name_second)
        return(N_metrics)
      }
      else
        {
         system(paste("/home/nancy/assembly_app/python3/Python-3.6.4/python /home/nancy/assembly_app/assembly_stats.py", example_file.fasta, output_assembly_metrics, estimated_genome_size), intern = TRUE)
        output <- read.table(output_assembly_metrics, quote="\"", comment.char="")
        metrics_fasta <- output$V1
        num_scaffolds <- metrics_fasta[1]
        total_size_scaffolds <- metrics_fasta[2]
        total_scaffold_length_as_precentage_genome_size <- metrics_fasta[3]
        useful_amount_scaffold_sequence <- metrics_fasta[4]
        percentage_estimated_genome_useful <- metrics_fasta[5]
        longest_scaffold <- metrics_fasta[6]
        shortest_scaffold <- metrics_fasta[7]
        scaffolds_greater_than_1k <- metrics_fasta[8]
        percentage_scaffolds_greater_than_1k <- metrics_fasta[9]
        scaffolds_greater_than_10k <- metrics_fasta[10]
        percentage_scaffolds_greater_than_10k <- metrics_fasta[11]
        scaffolds_greater_than_100k <- metrics_fasta[12]
        percentage_scaffolds_greater_than_100k <- metrics_fasta[13]
        scaffolds_greater_than_1M <- metrics_fasta[14]
        percentage_scaffolds_greater_than_1M <- metrics_fasta[15]
        scaffolds_greater_than_10M <- metrics_fasta[16]
        percentage_scaffolds_greater_than_10M <- metrics_fasta[17]
        N50 <- metrics_fasta[18]
        L50 <- metrics_fasta[19]
        NG50 <- metrics_fasta[20]
        LG50 <- metrics_fasta[21]
        percentage_A <- metrics_fasta[22]
        percentage_C <- metrics_fasta[23]
        percentage_G <- metrics_fasta[24]
        percentage_T <- metrics_fasta[25]
        Total_N <- metrics_fasta[26]
        percentage_N <- metrics_fasta[27]
        user_input <- c(num_scaffolds, total_size_scaffolds, total_scaffold_length_as_precentage_genome_size, useful_amount_scaffold_sequence, percentage_estimated_genome_useful, longest_scaffold, shortest_scaffold, scaffolds_greater_than_1k, scaffolds_greater_than_10k, scaffolds_greater_than_100k, scaffolds_greater_than_1M, scaffolds_greater_than_10M, N50, L50, NG50, LG50, percentage_A, percentage_C, percentage_G, percentage_T, Total_N, percentage_N)
        N_metrics <- cbind(user_input)
        #N_metrics <- as.matrix(user_input)
        row.names(N_metrics) <- c ("Number of scaffolds", "Total size of scaffolds", "Total scaffold length as percentage of assumed genome size", "useful amount of scaffold sequences (>=25K nt)", "% of estimated genome that is useful", "Longest scaffold", "Shortest scaffold", "Number of scaffolds > 1K nt", "Number of scaffolds > 10K nt", "Number of scaffolds > 100K nt", "Number of scaffolds > 1M nt", "Number of scaffolds > 10M nt", "N50", "L50", "NG50", "LG50", "%A", "%C", "%G", "%T", "Total Number of Ns","%N")
        colnames(N_metrics) <- c(input_name)
        return(N_metrics)
        }
      }
     }
  })
  
  # outputs the assembly metrics table to the assembly tab 
  output$lengthmetrics1 = DT::renderDataTable({
    
    ref_names <- input$references1
    
    input_name <- input$user_input1
    input_name_second <- input$user_input3

      if (!is.null(input$fasta_file1$datapath)) {
       if(!is.null(input$references1)) {
        if (!is.null(input$fasta_file3$datapath)) {
        cn1 <- c(ref_names, input_name, input_name_second)
           }
        else {
          cn1 <- c(ref_names, input_name)
             }
         }
      else {
        if (!is.null(input$fasta_file3$datapath)) {
        cn1 <- c(input_name, input_name_second)
         }
        else {
          cn1 <- c(input_name)
         }
       }
    }
#}
    
    
    data_fr <- data_f1()
    if (length(data_fr)==0)
      return(NULL)
    
    datatable(data_fr, selection = 'single', rownames = c ("Number of scaffolds", "Total size of scaffolds (bp)", "Total scaffold length as percentage of assumed genome size (%)", "useful amount of scaffold sequences (>=25K nt) (bp)", "% of estimated genome that is useful", "Longest scaffold (bp)", "Shortest scaffold(bp)", "Number of scaffolds > 1K nt", "Number of scaffolds > 10K nt", "Number of scaffolds > 100K nt", "Number of scaffolds > 1M nt", "Number of scaffolds > 10M nt", "N50 (bp)", "L50", "NG50 (bp)", "LG50", "%A", "%C", "%G", "%T", "Total Number of Ns", "%N"), colnames = cn1)
  })
  
  
  # generates pop-up plots on selecting the rows of assembly metrics table
  observeEvent(input$lengthmetrics1_rows_selected,
               
               {
                 
                 data_fr1 <- data_f1()
                 ref_names <- input$references1
                 new_labels_refnames <- c()
                 for (i in ref_names ) {
                   first_word <- unlist(strsplit(i, "_"))
                   new_labels_refnames <- c(new_labels_refnames,first_word[1])
                 } 
                 input_name <- input$user_input1
                 input_name_second <- input$user_input3
                 len <- length(ref_names) + 1
                       if (!is.null(input$fasta_file1$datapath)) {
                          if(!is.null(input$references1)) {
                            if (!is.null(input$fasta_file3$datapath)) {
                               cn2 <- c(new_labels_refnames, input_name, input_name_second)
                               }
                             else {
                               cn2 <- c(new_labels_refnames, input_name)
                               }
                              }
                            else {
                             if (!is.null(input$fasta_file3$datapath)) {
                                 cn2 <- c(input_name, input_name_second)
                                 }
                               else {
                                 cn2 <- c(input_name)
                                  }
                               }
                            }

                 
                 if (input$lengthmetrics1_rows_selected == 13 )
                 {
                   showModal(modalDialog(
                     title = "An N50 contig size of N means that 50% of the assembled bases are contained in contigs of length N or larger. N50 sizes are often used as a measure of assembly quality because they capture how much of the genome is covered by relatively large contigs.",
                     
                     
                     outputplot1 <- renderPlot({
                     plt <- barplot(as.numeric(data_fr1[13,]), main="N50", ylab="N50 value (bp)", col=rainbow(len))
                     text(plt, par("usr")[3], labels = cn2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                 })
                   ))
                 }
                 if (input$lengthmetrics1_rows_selected == 1 )
                 {
                   showModal(modalDialog(
                     title = "Total number of sequences in the genome assembly",
                     outputplot1 <- renderPlot({
                     g <- as.numeric(data_fr1[1,])
                     plt <- barplot(g ,ylab="Number of sequences", col=rainbow(len))
                     text(plt, par("usr")[3], labels = cn2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                    })
                     
                     
                   ))
                 }
                 if (input$lengthmetrics1_rows_selected == 2 )
                 {
                   showModal(modalDialog(
                     title = "Total size of the scaffolds",
                     
                    
                    outputplot1 <- renderPlot({
                    plt <- barplot(as.numeric(data_fr1[2,]) ,ylab="scaffold length (bp)", col=rainbow(len))
                    text(plt, par("usr")[3], labels = cn2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                    })
                   ))
                 }
                 
                if (input$lengthmetrics1_rows_selected == 3 )
                 {
                   showModal(modalDialog(
                     title = "Total scaffold length as percentage of assumed genome size (%)",
               
                    outputplot1 <- renderPlot({
                    plt <- barplot(as.numeric(data_fr1[3,]) , ylab="scaffold length (bp)", col=rainbow(len))
                      text(plt, par("usr")[3], labels = cn2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                      })
                    

                
                   ))
                 }
                 if (input$lengthmetrics1_rows_selected == 4)
                 {
                   showModal(modalDialog(
                     title = "useful amount of scaffold sequences (>=25K nt) (bp)",

                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr1[4,]) , ylab="scaffold length (bp)", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cn2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })


                   ))
                 }

                 if (input$lengthmetrics1_rows_selected == 5)
                 {
                   showModal(modalDialog(
                     title =  "% of estimated genome that is useful",
                     
                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr1[5,]), ylab="useful estimated genome size (%)", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cn2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })                    
                 
                     
                     
                   ))
                 }
                 
                 if (input$lengthmetrics1_rows_selected == 6)
                 {
                   showModal(modalDialog(
                     title = "Longest scaffold (bp)",
                     
                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr1[6,]) , ylab="scaffold length (bp)", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cn2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 
                 if (input$lengthmetrics1_rows_selected == 7)
                 {
                   showModal(modalDialog(
                     title = "Shortest scaffold (bp)",
                     
                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr1[7,]) , ylab="scaffold length (bp)", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cn2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$lengthmetrics1_rows_selected == 8)
                 {
                   showModal(modalDialog(
                     title = "Number of scaffolds > 1K nt",

                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr1[8,]) ,ylab="Number of scaffolds", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cn2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })


                   ))
                 }
                if (input$lengthmetrics1_rows_selected == 9)
                 {
                   showModal(modalDialog(
                     title = "Number of scaffolds > 10K nt",

                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr1[9,]) , ylab="Number of scaffolds", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cn2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })


                   ))
                 }

               if (input$lengthmetrics1_rows_selected == 10)
                 {
                   showModal(modalDialog(
                     title = "Number of scaffolds > 100K nt",

                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr1[10,]) , ylab="Number of scaffolds", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cn2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })


                   ))
                 }
               if (input$lengthmetrics1_rows_selected == 11)
                 {
                   showModal(modalDialog(
                     title =  "Number of scaffolds > 1M nt",

                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr1[11,]) , ylab="Number of scaffolds", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cn2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })


                   ))
                 }
               if (input$lengthmetrics1_rows_selected == 12)
                 {
                   showModal(modalDialog(
                     title = "Number of scaffolds > 10M nt",

                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr1[12,]) , ylab="Number of scaffolds", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cn2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })


                   ))
                 }
            
               if (input$lengthmetrics1_rows_selected == 14)
                 {
                   showModal(modalDialog(
                     title = "A scaffold/contig L50 is calculated as the number of sequences whose sum of lengths makes up 50% or more of the total assembly length.",

                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr1[14,]) , ylab="Number of scaffolds", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cn2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })


                   ))
                 }
               if (input$lengthmetrics1_rows_selected == 15)
                 {
                   showModal(modalDialog(
                     title = "NG50 is calculated in the same manner as N50 but uses estimated genome size instead of total assembly length, is more useful as it reduces the bias caused by assembly lengths.",

                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr1[15,]) , ylab="scaffold length (bp)", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cn2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })


                   ))
                 }
               if (input$lengthmetrics1_rows_selected == 16)
                 {
                   showModal(modalDialog(
                     title = "A scaffold/contig LG50 is calculated as the number of sequences whose sum of lengths makes up 50% or more of the estimated genome size.",

                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr1[16,]) , ylab="Number of scaffolds", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cn2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })


                   ))
                 }
               if (input$lengthmetrics1_rows_selected == 17)
                 {
                   showModal(modalDialog(
                     title = "Percentage A",

                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr1[17,]) , ylab="%A", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cn2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })


                   ))
                 }
               if (input$lengthmetrics1_rows_selected == 18)
                 {
                   showModal(modalDialog(
                     title = "Percentage C",

                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr1[18,]) , ylab="%C", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cn2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })


                   ))
                 }
               if (input$lengthmetrics1_rows_selected == 19)
                 {
                   showModal(modalDialog(
                     title = "Percentage G",

                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr1[19,]) , ylab="%G", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cn2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })


                   ))
                 }
                if (input$lengthmetrics1_rows_selected == 20)
                 {
                   showModal(modalDialog(
                     title = "Percentage T",

                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr1[20,]) , ylab="%T", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cn2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })


                   ))
                 }
                if (input$lengthmetrics1_rows_selected == 21)
                 {
                   showModal(modalDialog(
                     title = "Number of Ns: gives an estimate of the gaps in the assembled genome sequence",

                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr1[21,]) , ylab="Number of Ns", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cn2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })


                   ))
                 }
               if (input$lengthmetrics1_rows_selected == 22)
                 {
                   showModal(modalDialog(
                     title = "Percentage N",

                     outputplot1 <- renderPlot({
                       plt <- barplot(as.numeric(data_fr1[22,]) , ylab="%N", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cn2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })


                   ))
                 } 
               })
  
  
  # Download the assembly metrics table in csv format
   output$downloadassemblytable1 <- downloadHandler(
     
     filename = function() {
       
       u_num <- ran_number1()
       file_name = paste0("Assembly_metrics_table",u_num,".csv")
     },
     content = function(file) {
       
       write.csv(data_f1(), file)
    })
  

  
## THIS MODULE OUTPUTS THE REFERENCE GENOME ANNOTATION METRICS TABLE  
  
  data_gff <- eventReactive(input$go,{
    
    if (is.null(input$references))
    {
      return(NULL)}
    else{
      ref_names <- input$references
      
      # progress indicator
      progress <- Progress$new(session, min=1, max=15)
      on.exit(progress$close())
      progress$set(message = 'Calculation in progress', detail = 'This may take a while...')
      
      for (i in 1:15) {
        progress$set(value = i)
        Sys.sleep(0.5)
      }
      
      # reads the reference pre-computed values
      list1 = list()
      for(i in 1:length(ref_names)){
        file <- read.table(paste0("/home/nancy/assembly_app/reference_genomes/", ref_names[i], "_gmetrics"),quote="\"", comment.char="")
        list1[i] <- file
      }
      l <- as.vector(list1)
      
      # creates dataframe for the table
      X1 <- as.data.frame(l)
      colnames(X1) <- ref_names
      gff_metrics <- X1
      row.names(gff_metrics) <- c ("Number of gene models (bp)", "Minimum gene length (bp)", "Maximum gene length (bp)", "Average gene length (bp)", "Number of exons", "Average number of exons per gene model", "Average exon length (bp)", "Number of transcripts", "Average number of transcripts per gene model","Number of gene models less than 200bp length")
      data_gff_metrics <- as.data.frame(gff_metrics)
      return(data_gff_metrics)
    }
  })
  
  
  # outputs the reference annotation metrics table to the tab 
  output$annotationmetrics = DT::renderDataTable({
    
    ref_names <- input$references
    cn <- c(ref_names)
    gff_metrics_table <- data_gff()
    if (length(gff_metrics_table)==0)
      return(NULL)
    datatable(gff_metrics_table, selection = 'single', rownames = c ("Number of gene models", "Minimum gene length", "Maximum gene length", "Average gene length", "Number of exons", "Average number of exons per gene model", "Average exon length", "Number of transcripts", "Average number of transcripts per gene model","Number of gene models less than 200bp length"), colnames = cn)
  })
  
  
  # generates pop-up plots on selecting the rows of reference annotation metrics table
  observeEvent(input$annotationmetrics_rows_selected,
               
               {
                 data_ann <- data_gff()
                 ref_names <- input$references
                 len <- length(ref_names) + 1
                 cg <- c(ref_names)
                 new_labels <- c()
                 for (i in cg) {
                   first_word <- unlist(strsplit(i, "_"))
                   new_labels <- c(new_labels,first_word[1])
                 } 
                 
                 if (input$annotationmetrics_rows_selected == 1 )
                 {
                   showModal(modalDialog(
                     title = "Number of gene models",
                     
                     
                     outputplot2 <- renderPlot({
                       plt <- barplot(as.numeric(data_ann[1,]), ylab="Number of gene models", col=rainbow(len))
                     text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$annotationmetrics_rows_selected == 2 )
                 {
                   showModal(modalDialog(
                     title = "Minimum gene length",
                     
                     
                     outputplot2 <- renderPlot({
                       plt <- barplot(as.numeric(data_ann[2,]), ylab="Gene length (bp)", col=rainbow(len))
                      text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                     })
                     
                     
                   ))
                 }
                 if (input$annotationmetrics_rows_selected == 3 )
                 {
                   showModal(modalDialog(
                     title = "Maximum gene length",
                     
                     
                     outputplot2 <- renderPlot({
                       plt <- barplot(as.numeric(data_ann[3,]) , ylab="Gene length (bp)", col=rainbow(len))
                     text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                     })
                   ))
                 }
                 if (input$annotationmetrics_rows_selected == 4 )
                 {
                   showModal(modalDialog(
                     title = "Average gene length",
                     
                     
                     outputplot2 <- renderPlot({
                       plt <- barplot(as.numeric(data_ann[4,]) , ylab="Gene length (bp)", col=rainbow(len))
                       text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$annotationmetrics_rows_selected == 5 )
                 {
                   showModal(modalDialog(
                     title = "Number of exons",
                     
                     
                     outputplot2 <- renderPlot({
                       plt <- barplot(as.numeric(data_ann[5,]) , ylab="number of exons", col=rainbow(len))
                       text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$annotationmetrics_rows_selected == 6 )
                 {
                   showModal(modalDialog(
                     title = "Average number of exons per gene model",
                     
                     
                     outputplot2 <- renderPlot({
                       plt <- barplot(as.numeric(data_ann[6,]) , ylab="Average number of exons per gene model", col=rainbow(len))
                       text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$annotationmetrics_rows_selected == 7 )
                 {
                   showModal(modalDialog(
                     title = "Average exon length",
                     
                     
                     outputplot2 <- renderPlot({
                       plt <- barplot(as.numeric(data_ann[7,]) , ylab="Average exon length (bp)", col=rainbow(len))
                       text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$annotationmetrics_rows_selected == 8 )
                 {
                   showModal(modalDialog(
                     title = "Number of transcripts",
                     
                     
                     outputplot2 <- renderPlot({
                       plt <- barplot(as.numeric(data_ann[8,]) , ylab="Number of transcripts", col=rainbow(len))
                       text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$annotationmetrics_rows_selected == 9 )
                 {
                   showModal(modalDialog(
                     title = "Average number of transcripts per gene model",
                     
                     
                     outputplot2 <- renderPlot({
                       plt <- barplot(as.numeric(data_ann[9,]) , ylab="Average number of transcripts per gene model", col=rainbow(len))
                       text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$annotationmetrics_rows_selected == 10 )
                 {
                   showModal(modalDialog(
                     title = "Number of gene models less than 200bp length",
                     
                     
                     outputplot2 <- renderPlot({
                       plt <- barplot(as.numeric(data_ann[10,]) , ylab="Number of Gene models less than 200bp length", col=rainbow(len))
                       text(plt, par("usr")[3], labels = new_labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 
                 
               }
  )
  
  # Download the reference annotation metrics table in csv format
   output$downloadannotationtable <- downloadHandler(
     
     filename = function() {
       
       "Reference_annotation_metrics_table.csv"
     },
     content = function(file) {
       
       write.csv(data_gff(), file)
     })
 

  
## THIS MODULE OUTPUTS THE ANNOTATION METRICS TABLE FOR THE USER UPLOADED GENOME ANNOTATION AND COMPARE THE VALUES WITH PRE-COMPUTED REFERENCE VALUES
   
  data_gff1 <- eventReactive(input$go2,{
    
    # If the reference genomes are selected by the user
    if (!is.null(input$gff_file$datapath)) {
      inFile2 <- input$gff_file
      example_file.gff <- input_gff_file()
      if(!is.null(input$references2)) {
        
        # progress indicator
        progress <- Progress$new(session, min=1, max=15)
        on.exit(progress$close())
        progress$set(message = 'Calculation in progress', detail = 'This may take a while...')
        
        for (i in 1:15) {
          progress$set(value = i)
          Sys.sleep(0.5)
        }
        
        ref_names <- input$references2
        input_name <- input$user_input2
        u_num <- ran_number2()
        output_gff_metrics <- paste0("annotation", u_num,"_gffmetrics.txt")
        
        # This command calls the python script which takes two input arguments: input file name, output file name and calculate all the annotation metrics
        system(paste("/home/nancy/assembly_app/python3/Python-3.6.4/python /home/nancy/assembly_app/gff3_stats.py", example_file.gff, output_gff_metrics), intern = TRUE)
        output <- read.table(output_gff_metrics, quote="\"", comment.char="")
        metrics_gff <- output$V1
        number_gene_models <- metrics_gff[1]
        minimum_gene_length <- metrics_gff[2]
        maximum_gene_length <- metrics_gff[3]
        average_gene_length <- metrics_gff[4]
        number_exons <- metrics_gff[5]
        average_number_exons <- metrics_gff[6]
        average_exon_length <- metrics_gff[7]
        number_transcripts <- metrics_gff[8]
        average_number_transcripts <- metrics_gff[9]
        number_gene_models_lessthan200 <- metrics_gff[10]
        user_gff_metrics <- c(number_gene_models, minimum_gene_length, maximum_gene_length, average_gene_length, number_exons, average_number_exons, average_exon_length, number_transcripts, average_number_transcripts, number_gene_models_lessthan200)
        
        # reads the reference pre-computed values
        list1 = list()
        for(i in 1:length(ref_names)){
          file <- read.table(paste0("/home/nancy/assembly_app/reference_genomes/", ref_names[i], "_gmetrics"),quote="\"", comment.char="")
          list1[i] <- file
        }
        l <- as.vector(list1)
        
        # creates the dataframe for the table
        X1 <- as.data.frame(l)
        colnames(X1) <- ref_names
        gff_metrics <- cbind(X1, user_gff_metrics)
        colnames(gff_metrics) <- c(ref_names,"Your annotations")
        row.names(gff_metrics) <- c ("Number of gene models (bp)", "Minimum gene length (bp)", "Maximum gene length (bp)", "Average gene length (bp)", "Number of exons", "Average number of exons per gene model", "Average exon length (bp)", "Number of transcripts", "Average number of transcripts per gene model","Number of gene models less than 200bp length")
        colnames(gff_metrics) <- c(ref_names, input_name)
        data_gff_metrics <- as.data.frame(gff_metrics)
        return(data_gff_metrics)
        
      }
      
      # If the user does not select any reference genome 
      else{
        u_num <- ran_number2()
        input_name <- input$user_input2
        output_gff_metrics <- paste0("annotation", u_num,"_gffmetrics.txt")
        
        # This command calls the python script which takes two input arguments: input file name, output file name and calculate all the annotation metrics
        system(paste("/home/nancy/assembly_app/python3/Python-3.6.4/python /home/nancy/assembly_app/gff3_stats.py", example_file.gff, output_gff_metrics), intern = TRUE)
        output <- read.table(output_gff_metrics, quote="\"", comment.char="")
        metrics_gff <- output$V1
        number_gene_models <- metrics_gff[1]
        minimum_gene_length <- metrics_gff[2]
        maximum_gene_length <- metrics_gff[3]
        average_gene_length <- metrics_gff[4]
        number_exons <- metrics_gff[5]
        average_number_exons <- metrics_gff[6]
        average_exon_length <- metrics_gff[7]
        number_transcripts <- metrics_gff[8]
        average_number_transcripts <- metrics_gff[9]
        number_gene_models_lessthan200 <- metrics_gff[10]

        # creates the dataframe for the table
        user_gff_metrics <- c(number_gene_models, minimum_gene_length, maximum_gene_length, average_gene_length, number_exons, average_number_exons, average_exon_length, number_transcripts, average_number_transcripts, number_gene_models_lessthan200)
        gff_metrics <- as.matrix(user_gff_metrics)
        row.names(gff_metrics) <- c ("Number of gene models (bp)", "Minimum gene length (bp)", "Maximum gene length (bp)", "Average gene length (bp)", "Number of exons", "Average number of exons per gene model", "Average exon length (bp)", "Number of transcripts", "Average number of transcripts per gene model","Number of gene models less than 200bp length")
       colnames(gff_metrics) <- c(input_name) 
       return(gff_metrics)
      }
    }
    
  })
  
  # outputs the annotation metrics table to the tab 
  output$annotationmetrics1 = DT::renderDataTable({
    ref_names <- input$references2
    input_name <- input$user_input2
    
    if (!is.null(input$gff_file$datapath)) {
      if(!is.null(input$references2)) {
        cg1 <- c(ref_names, input_name)}
      
      else{
        cg1 <- input_name
      }
    }
    
    gff_metrics_table <- data_gff1()
    if (length(gff_metrics_table)==0)
      return(NULL)
    datatable(gff_metrics_table, selection = 'single', rownames = c ("Number of gene models", "Minimum gene length", "Maximum gene length", "Average gene length", "Number of exons", "Average number of exons per gene model", "Average exon length", "Number of transcripts", "Average number of transcripts per gene model","Number of gene models less than 200bp length"), colnames = cg1)
  })
  
  
  # generates pop-up plots on selecting the rows of the annotation metrics table
  observeEvent(input$annotationmetrics1_rows_selected,
               
               {
                 
                 data_ann <- data_gff1()
                 ref_names <- input$references2
                 new_labels_refnames <- c()
                 for (i in ref_names ) {
                   first_word <- unlist(strsplit(i, "_"))
                   new_labels_refnames <- c(new_labels_refnames,first_word[1])
                 } 
                 input_name <- input$user_input2
                 len <- length(ref_names) + 1
                 if (!is.null(input$gff_file$datapath)) {
                   if(!is.null(input$references2)) {
                     cg2 <- c(new_labels_refnames, input_name)}
                   
                   else{
                     cg2 <- input_name
                   }
                 }
               
                 
                 if (input$annotationmetrics1_rows_selected == 1 )
                 {
                   showModal(modalDialog(
                     title = "Number of gene models",
                     
                     
                     outputplot2 <- renderPlot({
                       plt <- barplot(as.numeric(data_ann[1,]), ylab="Number of gene models", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cg2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$annotationmetrics1_rows_selected == 2 )
                 {
                   showModal(modalDialog(
                     title = "Minimum gene length",
                     
                     
                     outputplot2 <- renderPlot({
                       plt <- barplot(as.numeric(data_ann[2,]), ylab="Gene length (bp)", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cg2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                     })
                     
                     
                   ))
                 }
                 if (input$annotationmetrics1_rows_selected == 3 )
                 {
                   showModal(modalDialog(
                     title = "Maximum gene length",
                     
                     
                     outputplot2 <- renderPlot({
                       plt <- barplot(as.numeric(data_ann[3,]) , ylab="Gene length (bp)", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cg2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$annotationmetrics1_rows_selected == 4 )
                 {
                   showModal(modalDialog(
                     title = "Average gene length",
                     
                     
                     outputplot2 <- renderPlot({
                       plt <- barplot(as.numeric(data_ann[4,]) , ylab="Gene length (bp)", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cg2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$annotationmetrics1_rows_selected == 5 )
                 {
                   showModal(modalDialog(
                     title = "Number of exons",
                     
                     
                     outputplot2 <- renderPlot({
                       plt <- barplot(as.numeric(data_ann[5,]) , ylab="number of exons", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cg2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$annotationmetrics1_rows_selected == 6 )
                 {
                   showModal(modalDialog(
                     title = "Average number of exons per gene model",
                     
                     
                     outputplot2 <- renderPlot({
                       plt <- barplot(as.numeric(data_ann[6,]) , ylab="Average number of exons per gene model", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cg2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$annotationmetrics1_rows_selected == 7 )
                 {
                   showModal(modalDialog(
                     title = "Average exon length",
                     
                     
                     outputplot2 <- renderPlot({
                       plt <- barplot(as.numeric(data_ann[7,]) , ylab="Average exon length (bp)", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cg2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$annotationmetrics1_rows_selected == 8 )
                 {
                   showModal(modalDialog(
                     title = "Number of transcripts",
                     
                     
                     outputplot2 <- renderPlot({
                       plt <- barplot(as.numeric(data_ann[8,]) , ylab="Number of transcripts", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cg2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$annotationmetrics1_rows_selected == 9 )
                 {
                   showModal(modalDialog(
                     title = "Average number of transcripts per gene model",
                     
                     
                     outputplot2 <- renderPlot({
                       plt <- barplot(as.numeric(data_ann[9,]) , ylab="Average number of transcripts per gene model", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cg2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 if (input$annotationmetrics1_rows_selected == 10 )
                 {
                   showModal(modalDialog(
                     title = "Number of gene models less than 200bp length",
                     
                     
                     outputplot2 <- renderPlot({
                       plt <- barplot(as.numeric(data_ann[10,]) , ylab="Number of Gene models less than 200bp length", col=rainbow(len))
                       text(plt, par("usr")[3], labels = cg2, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.9)
                       })
                     
                     
                   ))
                 }
                 
                 
               }
  ) 
  
     # Download the annotation metrics table in csv format
     output$downloadannotationtable1 <- downloadHandler(
       
       filename = function() {
       u_num <- ran_number2()
       file_name = paste0("Annotation_metrics_table",u_num,".csv")
       },
       content = function(file) {
         
         write.csv(data_gff1(), file)
       })
   
  
  
  
## THIS MODULE PLOTS AND EMAIL THE PRE-COMPUTED REFERENCE GENOME BUSCO SCORES
     reference_busco <- eventReactive(input$go,{

    if (is.null(input$references))
     {
      message2 = "ERROR: please provide both the required parameters: name of the reference genome and email address"
      return(message2)
     }else if (input$id == ""){
      message2 = "ERROR: please provide both the required parameters: name of the reference genome and email address"
      return(message2)
    }
      else{
      ref_names <- input$references
      REF_NAM <- paste((ref_names), collapse=",")
      email_id <- input$id
      u_num <- ran_number()
      output_busco_asssemblyplot <- paste0("assembly_busco", u_num)
      output_busco_gffplot <- paste0("annotation_busco", u_num)

      # This command calls the python script which takes four input arguments: user's email Id, output file name for busco assembly and annotation plots and the names of the reference genomes and plots and email the busco scores
      system(paste("/home/nancy/assembly_app/python3/Python-3.6.4/python /home/nancy/assembly_app/busco-master/scripts/png_no_fasta_gff.py", email_id, output_busco_asssemblyplot, output_busco_gffplot, REF_NAM), intern = TRUE)
      message2 = "Your job has been submitted"
      return(message2)

    }
  })


  # Prints the job submission message to the tab 
  output$reference_busco_plots <- renderPrint({

    m <- reference_busco()
    job_id <- ran_number()
    print_m <- paste0(m,". Your job ID is: ",job_id)
    #if(grepl("ERROR", print_m) == 'TRUE')
    if(grepl("ERROR", m) == 'TRUE')
    {
      HTML(paste0("<span style=\"color:red\">",m,"</span>"))
    } else {HTML(paste0("<span style=\"color:black\">",print_m,"</span>"))}

  }) 

  
## THIS MODULE CALCULATES, PLOTS, AND EMAIL THE BUSCO SCORES AND CONTAMINATION SCORE FOR THE USER UPLOADED GENOME ASSEMBLY AND COMPARE THE VALUES WITH PRE-COMPUTED REFERENCE VALUES  
  
  assembly_busco <- eventReactive(input$go1,{
    
    # If the reference genomes are selected by the user

      if (input$id1 == "") {
      message2 = "ERROR: please provide the email address"
      return(message2)
       }
      else {
      if (!is.null(input$references1)){
      ref_names <- input$references1
      email_id <- input$id1
      example_file.fasta <- input_fasta_file1()
      u_num <- ran_number1()
      output_name <- paste0("file_first", u_num,"_busco")
      db <- strsplit(input$busco_sets1, " ") [[1]]
      busco_name <- db[1]
      b_name <- paste0("/home/nancy/assembly_app/busco-master/scripts/", busco_name)
      sp_name <- input$augustus_species
      species_name <- paste0("--species=",sp_name)
      input_name <- input$user_input1
      ref_names <- input$references1
      REF_NAM <- paste((ref_names), collapse=",")
      
      #if the user uploads the second genome assembly fasta file
      if(!is.null(input$fasta_file3$datapath))
      {
        example_file2.fasta <- input_fasta_file3()
        output2_name <- paste0("file_second", u_num,"_busco")
        input_name3 <- input$user_input3 
     
      # This command calls the bash script which takes seven input arguments: input file, file name for the busco output directory, busco datasets, species name, email ID, input name to label plot, and the names of the reference genomes and calculates, plots and then email the busco scores and contamination results
      future(system(paste("bash /home/nancy/assembly_app/busco-master/scripts/combined.sh", example_file.fasta, output_name, b_name, species_name, email_id, input_name, REF_NAM, u_num, example_file2.fasta, input_name3, output2_name), intern = TRUE))
      message = "Your job has been submitted"
      return(message)
      }
      
      # if the user does not uploads the second genome assembly fasta file
      else
      {
       #This command calls the bash script which takes seven input arguments: input file, file name for the busco output directory, busco datasets, species name, email ID, input name to label plot, and the names of the reference genomes and calculates, plots and then email the busco scores and contamination results
       future(system(paste("bash /home/nancy/assembly_app/busco-master/scripts/one_assembly_contamination.sh", example_file.fasta, output_name, b_name, species_name, email_id, input_name, REF_NAM, u_num), intern = TRUE))
      message = "Your job has been submitted"
      return(message)
      }
      }

    # If the user does not select any reference genome 
    else{
      email_id <- input$id1
      example_file.fasta <- input_fasta_file1()
      
      u_num <- ran_number1()
      output_name <- paste0("file_first", u_num,"_busco")
      
      db <- strsplit(input$busco_sets1, " ") [[1]]
      busco_name <- db[1]
      b_name <- paste0("/home/nancy/assembly_app/busco-master/scripts/", busco_name)
      sp_name <- input$augustus_species
      species_name <- paste0("--species=",sp_name)
      input_name <- input$user_input1
      
      #if the user uploads the second assembly fasta file
      if(!is.null(input$fasta_file3$datapath))
      {
        example_file2.fasta <- input_fasta_file3()
        output2_name <- paste0("file_second", u_num,"_busco")
        input_name3 <- input$user_input3
      # This command calls the bash script which takes seven input arguments: input file, file name for the busco output directory, busco datasets, species name, email ID, name to label plot, and the names of the reference genomes and calculates, plots and then email the busco scores and contamination results
      future(system(paste("bash /home/nancy/assembly_app/busco-master/scripts/noref_combined.sh", example_file.fasta, output_name, b_name, species_name, email_id, input_name, u_num, example_file2.fasta, input_name3, output2_name), intern = TRUE))
      message = "Your job has been submitted"
      return(message)
      }
      # if the user does not upload the second assembly fasta file
      else 
      {
      #This command calls the bash script which takes seven input arguments: input file, file name for the busco output directory, busco datasets, species name, email ID, input name to label plot, and the names of the reference genomes and calculates, plots and then email the busco scores and contamination results
       future(system(paste("bash /home/nancy/assembly_app/busco-master/scripts/noref_one_assembly_contamination.sh", example_file.fasta, output_name, b_name, species_name, email_id, input_name, u_num), intern = TRUE))
      message = "Your job has been submitted"
      return(message)
      }
      }
     }
    #}
  })
  
  # Prints the job submission message to the tab
  output$assembly_busco_plot <- renderPrint({
    
    m <- assembly_busco()
    job_id <- ran_number1()
    print_m <- paste0(m,". Your job ID is: ",job_id)
    #if(grepl("ERROR", print_m) == 'TRUE')
    if(grepl("ERROR", m) == 'TRUE')
    {
      HTML(paste0("<span style=\"color:red\">",m,"</span>"))
    } else {HTML(paste0("<span style=\"color:black\">",print_m,"</span>"))}
    
  })
  
  
## THIS MODULE CALCULATES, PLOTS, AND EMAIL THE BUSCO SCORES FOR THE USER UPLOADED GENOME ANNOTATION AND COMPARE THE VALUES WITH PRE-COMPUTED REFERENCE VALUES  
  
  annotation_busco <- eventReactive(input$go2,{
    
    # If the reference genomes are selected by the user
    if(!is.null(input$references2)) 
    {
      ref_names <- input$references2
      email_id <- input$id2
      
      u_num <- ran_number2()
      db <- strsplit(input$busco_sets2, " ") [[1]]
      busco_name <- db[1]
      b_name <- paste0("/home/nancy/assembly_app/busco-master/scripts/", busco_name)
      input_name <- input$user_input2
      REF_NAM <- paste((ref_names), collapse=",")
      
      # If the user uploads the transcripts file, the code uses this file directly to calculate the busco scores
      if(!is.null(input$transcripts_file$datapath))
      {
        input_transcripts.fasta <- input_transcripts_file()
        output_namegff <- paste0("transcripts_file", u_num,"_buscogff")
        
        # This command calls a bash script which takes six input arguments: transcript fasta file, name of the busco output directory, busco datasets, email id, name to label plot, and the names of the reference genomes and calculates, plots and then email the busco scores 
        future(system(paste("bash /home/nancy/assembly_app/busco-master/scripts/annotation.sh", input_transcripts.fasta, output_namegff, b_name, email_id, input_name, REF_NAM, u_num), intern = TRUE))
        message = "Your job has been submitted"
        return(message)
      }
      
      # If the user does not provide the transcripts file, the code checks the formatting of the uploaded genome assembly and annotation file and extract the fasta sequences to calulate the busco scores
      else 
      {
        example_file.fasta <- input_fasta_file2()
        inFile2 <- input$gff_file
        example_file.gff <- input_gff_file()
        output_namegff <- paste0("gff", u_num,"_buscogff")
        transcripts_output <- paste0("trans", u_num,"_transcripts.fa")
        
        # This command calls the gffread program to extract fasta sequences of the gene models using the uploaded genome assembly and annotation file
        system(paste("/home/nancy/assembly_app/busco-master/scripts/gffread-0.9.12.Linux_x86_64/gffread -w", transcripts_output, "-g",  example_file.fasta, example_file.gff), intern = TRUE)
        file_t <- paste0("/opt/shiny-server/samples/sample-apps/assembly_statistics/", transcripts_output)
        
        # This command checks if the gffread generates an empty file due to incorrect formatting of the uploaded assembly and annotation file and prints error message to the tab
        if (file.info(file_t)$size == 0)
        {
          message = "ERROR: please check the format of your input files. Every contig or chromosome name found in the 1st column of the input GFF file must have a corresponding sequence entry in the header of input Fasta file and must contain transcript coordinates information."
          return(message)
        }
        # If the formatting is correct, it submits the job for busco analysis
        else 
        {
          
          # This command calls a bash script which takes six input arguments: transcripts fasta file generated by gffread, name of the busco output directory, busco datasets, email id, name to label plot, and the names of the reference genomes and calculates, plots and then email the busco scores 
          future(system(paste("bash /home/nancy/assembly_app/busco-master/scripts/annotation.sh", transcripts_output, output_namegff, b_name, email_id, input_name, REF_NAM, u_num), intern = TRUE))
          message = "Your job has been submitted"
          return(message)
        }
      }
      
    }
    
    # If the user does not select any reference genome 
    else 
    {
      email_id <- input$id2
      example_file.fasta <- input_fasta_file2()
      u_num <- ran_number2()
      db <- strsplit(input$busco_sets2, " ") [[1]]
      busco_name <- db[1]
      b_name <- paste0("/home/nancy/assembly_app/busco-master/scripts/", busco_name)
      input_name <- input$user_input2
      
      # If the user uploads the transcripts file, the code uses this file directly to calculate the busco scores
      if(!is.null(input$transcripts_file$datapath))
      {
        input_transcripts.fasta <- input_transcripts_file()
        output_namegff <- paste0("transcripts_file", u_num,"_buscogff")
        
        # This command calls a bash script which takes six input arguments: transcript fasta file, name of the busco output directory, busco datasets, email id, name to label plot, and calculates, plots and then email the busco scores 
        future(system(paste("bash /home/nancy/assembly_app/busco-master/scripts/noref_annotation.sh", input_transcripts.fasta, output_namegff, b_name, email_id, input_name, u_num), intern = TRUE))
        message = "Your job has been submitted"
        return(message)
      }
      
      # If the user does not provide the transcripts file, the code checks the formatting of the uploaded genome assembly and annotation file and extract the fasta sequences to calulate the busco scores
      else 
      {
        example_file.fasta <- input_fasta_file2()
        inFile2 <- input$gff_file
        example_file.gff <- input_gff_file()
        output_namegff <- paste0("gff", u_num,"_buscogff")
        transcripts_output <- paste0("trans", u_num,"_transcripts.fa")
        
        # This command calls the gffread program to extract fasta sequences of the gene models using the uploaded genome assembly and annotation file
        system(paste("/home/nancy/assembly_app/busco-master/scripts/gffread-0.9.12.Linux_x86_64/gffread -w", transcripts_output, "-g",  example_file.fasta, example_file.gff), intern = TRUE)
        file_t <- paste0("/opt/shiny-server/samples/sample-apps/assembly_statistics/", transcripts_output)
        
        # This command checks if the gffread generates an empty file due to incorrect formatting of the uploaded assembly and annotation file and prints error message to the tab
        if (file.info(file_t)$size == 0)
        {
          message = "ERROR: please check the format of your input files. Every contig or chromosome name found in the 1st column of the input GFF file must have a corresponding sequence entry in the header of input Fasta file and must contain transcript coordinates information."
          return(message)
        }
        # If the formatting is correct, it submits the job for busco analysis
        else 
        {
          # This command calls a bash script which takes six input arguments: transcripts fasta file generated by gffread, name of the busco output directory, busco datasets, email id, name to label plot, and calculates, plots and then email the busco scores 
          future(system(paste("bash /home/nancy/assembly_app/busco-master/scripts/noref_annotation.sh", transcripts_output, output_namegff, b_name, email_id, input_name, u_num), intern = TRUE))
          message = "Your job has been submitted"
          return(message)
        }
      }
    }
    
  })
  
  # Prints the job submission message to the tab
  output$annotation_busco_plot <- renderPrint({
    
    m <- annotation_busco()
    job_id <- ran_number2()
    print_m <- paste0(m,". Your job ID is: ",job_id)
    #if(grepl("ERROR", print_m) == 'TRUE')
    if(grepl("ERROR", m) == 'TRUE')
    {
      HTML(paste0("<span style=\"color:red\">",m,"</span>"))
    } else {HTML(paste0("<span style=\"color:black\">",print_m,"</span>"))}
    
  })
  
  
} 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
 
