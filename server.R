
### control size of input file
options(shiny.maxRequestSize=200*1024^2)

`%then%` <- shiny:::`%OR%`

server <- function(input, output, session) {
  
  #### Home Page & Resource Information
  output$resource_text1 <- renderText({
    message <- c("For a sample WHAM input file please visit the link below.")
    message
  })

  url1 <- a("Sample Input File", href="https://github.com/ruggleslab/jukebox/blob/master/wham_v1/sample_input.tsv.zip", target="_blank")
  output$samp_url <- renderUI({
    tagList("", url1)
  })
  
  
  output$resource_text2 <- renderText({
    paste("<br>", "Sample Data was derived from the HMP (Human Microbiome Project Consortium (2012) Structure, function and diversity of the healthy human microbiome. Nature, 486, 207–214.)")
  })
  
  url2 <- a("HMP HomePage", href="https://hmpdacc.org/hmp/", target="_blank")
  output$hmp_url <- renderUI({
    tagList("", url2)
  })
  
  output$humann2wham_text <- renderText({
    paste("<br>", "HUMANn2 users can readily convert HUMANn2 gene tables 
          to a wham-friendly input using our collection of python conversion scripts,
          available on the Ruggles Lab Github.")
  })
  
  url3 <- a("HUMANn2 File Converter", href="https://github.com/ruggleslab/jukebox/tree/master/wham_v1/file_conversion_scripts", target="_blank")
  output$humann2wham_url <- renderUI({
    tagList("", url3)
  })
  
   
  output$ext_text <- renderText({
    message <- c("WHAM! is an open-source project developed in the Ruggles Lab at NYU Langone Medical Center.")
  })
  
  url4 <- a("Ruggles Lab HomePage", href="http://www.ruggleslab.org/home.html", target="_blank")
  output$ruggles_url <- renderUI({
    tagList("", url4)
  })
  
  output$github <- renderText({
    paste("<br>", "Source Code can be found on the Ruggles Lab Github")
  })
  
  url5 <- a("Ruggles Lab Github", href="https://github.com/ruggleslab/jukebox", target="_blank")
  output$github_url <- renderUI({
    tagList("", url5)
  })
  
  output$citation <- renderText({
    paste("<br>", "For additional information of using WHAM! or to cite WHAM! in your work, please refer to the following paper:",
          "*Citation Info*")
  })
  
  observe({
    if (input$input_type == "Biobakery"){
      output$file_selector <- renderUI({
        fileInput('file1', 'Choose TSV File (Max=200MB)',
                  accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'))
      })
    }
    if (input$input_type == "EBI"){
      output$file_selector <- renderUI({
        fluidRow(
          fileInput('features1', 'Choose Feature File (Max=200MB)',
                    accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
          fileInput('taxa1', 'Choose Taxa File (Max=200MB)',
                    accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'))
        )
      })
    }
  })
  
  
  ## Load in the data or try a sample dataset
  full_file_feature <- reactive({
    if (input$testme) {
      start_col = 4
      #full_file <- fread("sample_input2.tsv", header=TRUE, sep=input$sep)
      full_file_feature <- fread("Relab_go_names_wham.tsv", header=TRUE, sep=input$sep) ##temp test
    }
    else {
      if (input$input_type == "Biobakery"){
        start_col = 4
        inFile <- input$file1
        if (is.null(inFile)) {
          return(NULL)}
        full_file <- try(
          {fread(inFile$datapath, header=TRUE, sep=input$sep)})
        
        correct_cols <- c("Acc", "Feature", "Taxa")
        validate(
          need(class(full_file)!="try-error", paste0("The input file provided is not an appropriate format. Please view the sample input file provided in the Home Tab.")) %then%
            need(all(colnames(full_file)[1:3]==correct_cols), paste0("The input file provided is not an appropriate format. Please view the sample input file provided in the Home Tab."))
        )
        
        nums <- data.matrix(full_file[,start_col:ncol(full_file)])
        rownames(nums) <- rownames(full_file)
        if (input$filter_level == 0){
          full_file <- full_file
        }
        else{
          keep_rows = rownames((nums[apply(nums==0,1,sum)<=input$filter_level*ncol(nums),]))
          full_file <- full_file[as.numeric(keep_rows),]
        }
        full_file <- subset(full_file, Feature != "NO_NAME")
        
        full_file$Feature <- gsub("[^[:alnum:]']", "_", full_file$Feature)
        
        full_file_feature <- full_file
      }
      if (input$input_type == "EBI"){
        start_col = 3
        inFile_taxa <- input$taxa1
        inFile_feature <- input$features1
        
        if (is.null(inFile_feature)) {
          return(NULL)}
        full_file_feature <- try(
          {fread(inFile_feature$datapath, header=TRUE, sep=input$sep)})
        
        correct_cols <- c("Acc", "Feature")
        validate(
          need(class(full_file_feature)!="try-error", paste0("The input file provided is not an appropriate format. Please view the sample input file provided in the Home Tab.")) %then%
            need(all(colnames(full_file_feature)[1:2]==correct_cols), paste0("The input file provided is not an appropriate format. Please view the sample input file provided in the Home Tab."))
        )
        
        nums <- data.matrix(full_file_feature[,start_col:ncol(full_file_feature)])
        
        nums_relab <- sweep(nums, 2, colSums(nums), '/')
        rownames(nums_relab) <- rownames(full_file_feature)
        
        full_file_feature <- cbind(full_file_feature[,1:(start_col-1)], nums_relab)
        
        if (input$filter_level == 0){
          full_file_feature <- full_file_feature
        }
        else{
          keep_rows = rownames((nums_relab[apply(nums_relab==0,1,sum)<=input$filter_level*ncol(nums_relab),]))
          full_file_feature <- full_file_feature[as.numeric(keep_rows),]
        }
        full_file_feature <- subset(full_file_feature, Feature != "NO_NAME")
        
        full_file_feature$Feature <- gsub("[^[:alnum:]']", "_", full_file_feature$Feature)
        
        full_file_feature
      }
    }
    full_file_feature
  })
  
  full_file_taxa <- reactive({
    ## now taxa
    if (input$input_type == "EBI"){
      start_col = 2
      inFile_taxa <- input$taxa1

      if (is.null(inFile_taxa)) {
        return(NULL)}
    full_file <- try(
      {fread(inFile_taxa$datapath, header=TRUE, sep=input$sep)})
    
    correct_cols <- c("Taxa")
    validate(
      need(class(full_file)!="try-error", paste0("The input file provided is not an appropriate format. Please view the sample input file provided in the Home Tab.")) %then%
        need(all(colnames(full_file)[1]==correct_cols), paste0("The input file provided is not an appropriate format. Please view the sample input file provided in the Home Tab."))
    )
    
    nums <- data.matrix(full_file[,start_col:ncol(full_file)])
    nums <- sweep(nums, 2, colSums(nums), '/')
    rownames(nums) <- rownames(full_file)
    
    full_file <- cbind(full_file[,1:(start_col-1)], nums)
    
    # if (input$filter_level == 0){
    #   full_file <- full_file
    # }
    # else{
    #   keep_rows = rownames((nums[apply(nums==0,1,sum)<=input$filter_level*ncol(nums),]))
    #   full_file <- full_file[as.numeric(keep_rows),]
    # }
    full_file <- subset(full_file, Taxa != "NO_NAME")
    
    full_file$Taxa <- gsub(";k__", "_k__", full_file$Taxa)
    #full_file$Taxa <- gsub("[^[:alnum:]']", "_", full_file$Taxa)
    
    full_file
    }
  })
  
  # Preview Full File
  output$contents <- renderDataTable({
    validate(
      need(full_file_feature(),"")
    )
    full_file <- full_file_feature()
    if (nrow(full_file)<50){
      full_file_show <- full_file[,1:ncol(full_file)]
    }
    else{
      full_file_show <- full_file[1:50, 1:ncol(full_file)]
    }
    full_file_show
  })
  
  output$contents_taxa <- renderDataTable({
    validate(
      need(full_file_taxa(),"")
    )
    full_file <- full_file_taxa()
    if (nrow(full_file)<50){
      full_file_show <- full_file[,1:ncol(full_file)]
    }
    else{
      full_file_show <- full_file[1:50, 1:ncol(full_file)]
    }
    full_file_show
  })
  
  observe({
    if (input$input_type == "Biobakery"){
      output$preview_shower <- renderUI({
        dataTableOutput('contents')
      })
    }
    if (input$input_type == "EBI"){
      output$preview_shower <- renderUI({
        fluidRow(
          dataTableOutput('contents'),
          dataTableOutput('contents_taxa')
        )
      })
    }
  })
  
  observe({
    if (input$testme){
      updateRadioButtons(session, "input_type", selected = "Biobakery")
    }
  })
  
  ### Generate prerequisites for plotting

  # Generate dataframe collapsed by Gene Family
  acc_full <- reactive({
    full_file <- full_file_feature()
    col_num <- ncol(full_file)
    if (input$testme){
      start_col = 4
      
      DT <- data.frame(full_file[,start_col:col_num], check.names = F)
      DT2 <- aggregate(DT, list(full_file$Feature), sum)
      colnames(DT2)[1] = "Feature"
      DT2 = subset(DT2, Feature != "filler")
      DT2
    } else {
      if (input$input_type == "Biobakery"){
        start_col = 4
        
        DT <- data.frame(full_file[,start_col:col_num], check.names = F)
        DT2 <- aggregate(DT, list(full_file$Feature), sum)
        colnames(DT2)[1] = "Feature"
        DT2 = subset(DT2, Feature != "filler")
        DT2
      }
      if (input$input_type == "EBI"){
        start_col = 3
        
        DT <- data.frame(full_file[,start_col:col_num], check.names = F)
        DT2 <- aggregate(DT, list(full_file$Feature), sum)
        colnames(DT2)[1] = "Feature"
        DT2 = subset(DT2, Feature != "filler")
        DT2
      }
    }
    DT2
  })
  # acc full but numeric values only
  acc_nums <- reactive({
    accs <- acc_full()
    accs2 = accs[,-1]
    rownames(accs2) = accs$Feature
    accs3 = as.matrix(accs2)
    accs3
  })

  # Generate dataframe collapsed by Species
  spec_full <- reactive({
    if (input$testme){
      start_col = 4
      full_file <- full_file_feature()
      col_num <- ncol(full_file)
      
      DT <- data.frame(full_file[,start_col:col_num], check.names=F)
      DT2 <- aggregate(DT, list(full_file$Taxa), sum)
      colnames(DT2)[1] = "Taxa"
      DT2 = subset(DT2, Taxa != "filler")
      DT2
    }
    if (input$input_type == "Biobakery"){
      start_col = 4
      full_file <- full_file_feature()
      col_num <- ncol(full_file)
      
      DT <- data.frame(full_file[,start_col:col_num], check.names=F)
      DT2 <- aggregate(DT, list(full_file$Taxa), sum)
      colnames(DT2)[1] = "Taxa"
      DT2 = subset(DT2, Taxa != "filler")
      DT2
    }
    if (input$input_type == "EBI"){
      start_col = 2
      validate(
        need(full_file_taxa(), "Taxa File Uploads are required for EBI inputs")
      )
      full_file <- full_file_taxa()
      col_num <- ncol(full_file)
      
      DT <- data.frame(full_file[,start_col:col_num], check.names=F)
      DT2 <- aggregate(DT, list(full_file$Taxa), sum)
      colnames(DT2)[1] = "Taxa"
      DT2 = subset(DT2, Taxa != "filler")
      DT2
    }
    DT2
  })

  spec_nums <- reactive({
    accs <- spec_full()
    accs2 = accs[,-1]
    rownames(accs2) = accs$Taxa
    accs3 = as.matrix(accs2)
    rownames(accs3) <- rownames(accs2)
    accs3
  })

  # Generate selection list of gene families
  acc_select <- reactive({
    accs <- acc_full()
    #accs$Feature = paste(" ", accs$Feature, " ", sep = "") #ensure unique gene families
    accs$Feature
  })
   
  # When file is uploaded update choice of gene families
  observe({
    if (input$testme) {
      updateSelectizeInput(session,'acc_list', choices = acc_select(), server = TRUE)
      updateSelectizeInput(session,'excluder', choices = acc_select(), server = TRUE)
      updateNumericInput(session, 'taxaDims', value = 2, min = 1, max = 2, step = 1)
      updateSelectizeInput(session,'sig_select', choices = sig_tab(), server = TRUE)
    }
    else {
      if (input$input_type == "Biobakery"){
        if(is.null(input$file1)){}
        else {
          updateSelectizeInput(session,'acc_list', choices = acc_select(), server = TRUE)
          updateSelectizeInput(session,'excluder', choices = acc_select(), server = TRUE)
          updateNumericInput(session, 'taxaDims', value = 2, min = 1, max = 2, step = 1)
          updateSelectizeInput(session,'sig_select', choices = sig_tab(), server = TRUE)
        }
      }
      if (input$input_type == "EBI"){
        if(is.null(input$features1)){}
        else {
          updateSelectizeInput(session,'acc_list', choices = acc_select(), server = TRUE)
          updateSelectizeInput(session,'excluder', choices = acc_select(), server = TRUE)
          updateNumericInput(session, 'taxaDims', value = 5, min = 1, max = 7, step = 1)
          updateSelectizeInput(session,'sig_select', choices = sig_tab(), server = TRUE)
        }
      }
    }
  })

  ## Group Selection Elements ##


  ###input grouped samples based on number of groups, output reordered matrix###
  results <- c()
  makeReactiveBinding('results')

  # create number of  columns based on input number of groups
  observe({
    output$inputGroup = renderUI({
      if (input$testme){
        updateNumericInput(session, 'numInputs', value = 4)
        group1 = c("SRR532024_Abundance.RPKs", "SRR532015_Abundance.RPKs", "SRR532006_Abundance.RPKs",
                   "SRR532040_Abundance.RPKs", "SRR532504_Abundance.RPKs", "SRR532507_Abundance.RPKs",
                   "SRR638753_Abundance.RPKs", "SRR640357_Abundance.RPKs", "SRR640340_Abundance.RPKs",
                   "SRR640452_Abundance.RPKs", "SRR640499_Abundance.RPKs", "SRR545546_Abundance.RPKs")
        group2 = c("SRR062353_Abundance.RPKs", "SRR062357_Abundance.RPKs", "SRR062301_Abundance.RPKs",
                   "SRR062276_Abundance.RPKs", "SRR1804686_Abundance.RPKs", "SRR1804628_Abundance.RPKs",
                   "SRR514191_Abundance.RPKs", "SRR514180_Abundance.RPKs", "SRR513168_Abundance.RPKs",
                   "SRR514231_Abundance.RPKs", "SRR513448_Abundance.RPKs")
        group3 = c("SRR062435_Abundance.RPKs", "SRR062441_Abundance.RPKs", "SRR062389_Abundance.RPKs",
                   "SRR062413_Abundance.RPKs", "SRR062402_Abundance.RPKs", "SRR062396_Abundance.RPKs",
                   "SRR346673_Abundance.RPKs", "SRR346681_Abundance.RPKs", "SRR062371_Abundance.RPKs",
                   "SRR062372_Abundance.RPKs", "SRR062462_Abundance.RPKs", "SRR062415_Abundance.RPKs")
        group4 = c("SRR528423_Abundance.RPKs", "SRR528353_Abundance.RPKs", "SRR528300_Abundance.RPKs",
                   "SRR528261_Abundance.RPKs", "SRR528183_Abundance.RPKs", "SRR528155_Abundance.RPKs",
                   "SRR532178_Abundance.RPKs", "SRR532183_Abundance.RPKs", "SRR532190_Abundance.RPKs",
                   "SRR532191_Abundance.RPKs", "SRR533152_Abundance.RPKs", "SRR533153_Abundance.RPKs")
        fluidRow(selectInput("test_groups", label = "",choices = group1,
                             multiple = TRUE, selected = group1),
                 selectInput("test_groups", label = "", choices = group2,
                             multiple = TRUE, selected = group2),
                 selectInput("test_groups", label = "", choices = group3,
                             multiple = TRUE, selected = group3),
                 selectInput("test_groups", label = "", choices = group4,
                             multiple = TRUE, selected = group4))
      } else {
        if (input$input_type == "Biobakery"){
          validate(
            need(input$file1, 'Please provide a file in the Upload Tab'))
        }
        if (input$input_type == "EBI"){
          validate(
            need(input$features1, 'Please provide a feature file in the Upload Tab'))
        }
        input_list <- lapply(1:input$numInputs, function(i) {
          inputName <- paste0("input", i)
          sampleUploadUI(inputName)
          })
      }
    })
    # Make list of resulting group allocations
    results <<- lapply(1:input$numInputs, function(i) {
      inputName <- paste0("input", i)
      callModule(sampleUpload, inputName, acc_nums)
    })
  })

  # Generate Space to Label Groups (based on input$numInputs)
  output$group_pre = renderUI({
    lapply(1:input$numInputs, function(i) {
      inputName <- paste0("group", i)
      textInput(inputName, label = "Group Prefix")
    })
  })


  # Ensure groups are named even if not input by user
  group_names <- reactive({
    if (input$testme){
      groups <- c("Arm", "Vagina", "Saliva", "Stool")
    }
    else{
      groups <- sapply(1:input$numInputs, function(i){
        cc<-input[[paste0("group", i)]][1]
        cc
      })
    }
    return(groups)
  })

  new_group_names <- reactive({
    #validate(need(group_names(), ""))
    groups <- sapply(1:input$numInputs, function(i){
      cc <- group_names()[i]
      if (cc==''){
        cc = paste0("Group", i)
      }
      cc
    })
    return(groups)
  })


  #Retain group allocations
  grouped_samps <- reactive({
    if (input$testme) {
      group1 = c("SRR532024_Abundance.RPKs", "SRR532015_Abundance.RPKs", "SRR532006_Abundance.RPKs",
                 "SRR532040_Abundance.RPKs", "SRR532504_Abundance.RPKs", "SRR532507_Abundance.RPKs",
                 "SRR638753_Abundance.RPKs", "SRR640357_Abundance.RPKs", "SRR640340_Abundance.RPKs",
                 "SRR640452_Abundance.RPKs", "SRR640499_Abundance.RPKs", "SRR545546_Abundance.RPKs")
      group2 = c("SRR062353_Abundance.RPKs", "SRR062357_Abundance.RPKs", "SRR062301_Abundance.RPKs",
                 "SRR062276_Abundance.RPKs", "SRR1804686_Abundance.RPKs", "SRR1804628_Abundance.RPKs",
                 "SRR514191_Abundance.RPKs", "SRR514180_Abundance.RPKs", "SRR513168_Abundance.RPKs",
                 "SRR514231_Abundance.RPKs", "SRR513448_Abundance.RPKs")
      group3 = c("SRR062435_Abundance.RPKs", "SRR062441_Abundance.RPKs", "SRR062389_Abundance.RPKs",
                 "SRR062413_Abundance.RPKs", "SRR062402_Abundance.RPKs", "SRR062396_Abundance.RPKs",
                 "SRR346673_Abundance.RPKs", "SRR346681_Abundance.RPKs", "SRR062371_Abundance.RPKs",
                 "SRR062372_Abundance.RPKs", "SRR062462_Abundance.RPKs", "SRR062415_Abundance.RPKs")
      group4 = c("SRR528423_Abundance.RPKs", "SRR528353_Abundance.RPKs", "SRR528300_Abundance.RPKs",
                 "SRR528261_Abundance.RPKs", "SRR528183_Abundance.RPKs", "SRR528155_Abundance.RPKs",
                 "SRR532178_Abundance.RPKs", "SRR532183_Abundance.RPKs", "SRR532190_Abundance.RPKs",
                 "SRR532191_Abundance.RPKs", "SRR533152_Abundance.RPKs", "SRR533153_Abundance.RPKs")
      g1 = list(group1)
      g2 = list(group2)
      g3 = list(group3)
      g4 = list(group4)
      glist = list(g1,g2,g3,g4)
      glist
    } else {
      #print(acc_nums()[1:2,1:2])
      result_call <- try(lapply(1:input$numInputs, function(i) {
        results[[i]]()}))
      result_call
    }
  })

  ## Store dimensions of each group for plottng
  group_dims <- reactive({
    req(grouped_samps())
    tl <- sapply(1:input$numInputs, function(i){
      sapply(grouped_samps()[[i]], length)
    })
    sample_num <- c(tl)
    return(sample_num)
  })

  # Reorder input matrix based on group allocations
  reorder_mat <- reactive({
    # validate(
    #   need(input$file1, "") %then%
    #     need(!is.null(grouped_samps()), "")
    #   )
    req(grouped_samps())
    exprs_reorder = acc_nums()[,c(unlist(grouped_samps()))]
    exp2 = data.frame(exprs_reorder)
    rownames(exp2) <- acc_full()$Feature
    return(exp2)
  })

  # Generate warning if sample is assigned to multiple groups
  #observeEvent(input$input1, {
  output$group_warning <- renderText({
    validate(
      need(length(unlist(grouped_samps()))>0, 'Select groups!')
    )
    group_list = c(unlist(grouped_samps()))
    duplicates = c(duplicated(group_list))
    dups <- group_list[which(duplicates==TRUE)]
    if ('TRUE' %in% duplicates) {
      message = paste0("Warning: ", dups, " assigned to more than one group!")
    } else {
      message = c(" ")
    }
    message
  })
  #})

  output$tutorialGroup <- renderText({
    if (input$testme){
      message = c("Groups have been named and assigned!\n Please continue to the next Tab")
      message
    }
  })

  #### Begin Plotting ! ####

  #
  #
  #
  #
  #
  #
  ## Ultimate ColorBar
  color_select <- reactive({
    groupings = new_group_names()
    grouping_nums = group_dims()
    
    bin <- c()
    for (i in 1:length(groupings)){
      use <- rep(groupings[i], grouping_nums[i])
      bin <- c(bin, use)
    }
    
    ##### make fake color bar to add to the plotly object!
    names <- bin
    uniq_names <- unique(names)
    
    cols_keep2 <- c("#f70c1c", "#6d89e8", "#91dd68", "#482e91", "#fc9207",
                    "#fcdb23", "#e87ac3", "#5beabd", "#01871e", "#a0080f")
    
    cols_cols <- cols_keep2[1:length(uniq_names)]
    cols_cols
    
  })
  
  colorer <- reactive({
    groupings = new_group_names()
    grouping_nums = group_dims()
    
    cols_cols <- color_select()
    
    bin <- c()
    for (i in 1:length(groupings)){
      use <- rep(groupings[i], grouping_nums[i])
      bin <- c(bin, use)
    }
    
    ##### make fake color bar to add to the plotly object!
    names <- bin
    uniq_names <- unique(names)
    
    meta <- data.frame(names)
    
    for (i in 1:length(uniq_names)){
      new_bin <- as.integer(meta$names == uniq_names[i])
      meta <- cbind(meta, new_bin)
    }
    
    series_mat <- meta[,-1]
    
    if (input$numInputs > 1){
      
      colnames(series_mat) <- uniq_names
      #series_mat <- t(series_mat)
      series_mat[series_mat==0] <- NA
      
      #cols_cols <- randomColor(length(uniq_names))
      
      g1 = t(as.matrix(series_mat[,1]))
      g1col <- data.frame(x = c(0,1), y = c(cols_cols[1], cols_cols[1]))
      colnames(g1col) <- NULL
      
      ax <- list(
        title = "",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        tickcolor = 'white',
        showgrid = FALSE
      )
      
      colorer <- plot_ly(
        type = "heatmap"
      ) %>% add_trace(
        x=1:length(bin),
        z = g1,
        colorscale = g1col, showscale = F,
        hoverinfo = 'all'
      ) %>%
        layout(yaxis = ax, xaxis = ax)
      
      for (i in 2:length(uniq_names)){
        g = t(as.matrix(series_mat[,i]))
        gcol <- data.frame(x = c(0,1), y = c(cols_cols[i], cols_cols[i]))
        colnames(gcol) <- NULL
        
        colorer <- colorer %>% add_trace(
          x=1:length(bin),
          z = g,
          colorscale = gcol, showscale = F,
          hoverinfo = 'all'
        )
      }
    }
    else {
      names(series_mat) <- uniq_names
      #series_mat <- t(series_mat)
      series_mat[series_mat==0] <- NA
      
      #cols_cols <- randomColor(length(uniq_names))
      
      g1 = t(as.matrix(series_mat))
      g1col <- data.frame(x = c(0,1), y = c(cols_cols[1], cols_cols[1]))
      colnames(g1col) <- NULL
      
      ax <- list(
        title = "",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        tickcolor = 'white',
        showgrid = FALSE
      )
      
      colorer <- plot_ly(
        type = "heatmap"
      ) %>% add_trace(
        x=1:length(bin),
        z = g1,
        colorscale = g1col, showscale = F,
        hoverinfo = 'all'
      ) %>%
        layout(yaxis = ax, xaxis = ax)
    }
    colorer
  })
  
  
  ##### color bar legend element!
  
  colorer_key <- reactive({
    validate(
      need(!is.null(new_group_names()[[1]]), "")
    )
    #print(is.null(new_group_names()[[1]]))
    groupings = new_group_names()
    cols_cols <- color_select()
    
    n <- length(groupings)
    Group = groupings
    
    plotplot <- data.frame(Group,cols_cols)
    plotplot$Group <- factor(plotplot$Group, levels = groupings)
    
    gg = ggplot(plotplot, aes(1:n, 1:n, color = Group)) +
      geom_point() +
      guides(color = guide_legend(direction = "horizontal",
                                  override.aes = list(shape = 15, size=9),
                                  ncol  = 4))+
      scale_color_manual(values = cols_cols) + 
      theme(legend.key = element_rect(fill = "transparent"),
            legend.text = element_text(size = 40),
            legend.title = element_text(size = 42))
    
    legend <- g_legend(gg) 
    grid.arrange(legend) 
  })
  
  
  output$key1_plot <- renderPlot({
    validate(
      need(!is.null(new_group_names()[[1]]), "")
    )
    #print(is.null(new_group_names()[[1]]))
    groupings = new_group_names()
    cols_cols <- color_select()
    
    n <- length(groupings)
    Group = groupings
    
    plotplot <- data.frame(Group,cols_cols)
    plotplot$Group <- factor(plotplot$Group, levels = groupings)
    
    gg = ggplot(plotplot, aes(1:n, 1:n, color = Group)) +
      geom_point() +
      guides(color = guide_legend(direction = "horizontal",
                                  override.aes = list(shape = 15, size=9),
                                  ncol  = 4))+
      scale_color_manual(values = cols_cols) + 
      theme(legend.key = element_rect(fill = "transparent"),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 18))
    
    legend <- g_legend(gg) 
    grid.draw(legend) 
  }, bg="transparent")
  
  output$key1 <- renderUI({
    pheight = paste0('50px')
    
    if (input$numInputs > 4){
      pheight = paste0('100px')
    }
    if (input$numInputs > 8){
      pheight = paste0('150px')
    }
    
    plotOutput("key1_plot", height = pheight)
  })
  
  output$key2_plot <- renderPlot({
    validate(
      need(!is.null(new_group_names()[[1]]), "")
    )
    #print(is.null(new_group_names()[[1]]))
    groupings = new_group_names()
    cols_cols <- color_select()
    
    n <- length(groupings)
    Group = groupings
    
    plotplot <- data.frame(Group,cols_cols)
    plotplot$Group <- factor(plotplot$Group, levels = groupings)
    
    gg = ggplot(plotplot, aes(1:n, 1:n, color = Group)) +
      geom_point() +
      guides(color = guide_legend(direction = "horizontal",
                                  override.aes = list(shape = 15, size=9),
                                  ncol  = 4))+
      scale_color_manual(values = cols_cols) + 
      theme(legend.key = element_rect(fill = "transparent"),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 18))
    
    legend <- g_legend(gg) 
    grid.draw(legend) 
  }, bg="transparent")
  
  output$key2 <- renderUI({
    pheight = paste0('50px')
    
    if (input$numInputs > 4){
      pheight = paste0('100px')
    }
    if (input$numInputs > 8){
      pheight = paste0('150px')
    }
    
    plotOutput("key2_plot", height = pheight)
  })
  
  output$key3_plot <- renderPlot({
    validate(
      need(!is.null(new_group_names()[[1]]), "")
    )
    #print(is.null(new_group_names()[[1]]))
    groupings = new_group_names()
    cols_cols <- color_select()
    
    n <- length(groupings)
    Group = groupings
    
    plotplot <- data.frame(Group,cols_cols)
    plotplot$Group <- factor(plotplot$Group, levels = groupings)
    
    gg = ggplot(plotplot, aes(1:n, 1:n, color = Group)) +
      geom_point() +
      guides(color = guide_legend(direction = "horizontal",
                                  override.aes = list(shape = 15, size=9),
                                  ncol  = 4))+
      scale_color_manual(values = cols_cols) + 
      theme(legend.key = element_rect(fill = "transparent"),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 18))
    
    legend <- g_legend(gg) 
    grid.draw(legend)
  }, bg="transparent")
  
  output$key3 <- renderUI({
    pheight = paste0('50px')
    
    if (input$numInputs > 4){
      pheight = paste0('100px')
    }
    if (input$numInputs > 8){
      pheight = paste0('150px')
    }
    
    plotOutput("key3_plot", height = pheight)
  })
  
  output$key4_plot <- renderPlot({
    validate(
      need(!is.null(new_group_names()[[1]]), "")
    )
    #print(is.null(new_group_names()[[1]]))
    groupings = new_group_names()
    cols_cols <- color_select()
    
    n <- length(groupings)
    Group = groupings
    
    plotplot <- data.frame(Group,cols_cols)
    plotplot$Group <- factor(plotplot$Group, levels = groupings)
    
    gg = ggplot(plotplot, aes(1:n, 1:n, color = Group)) +
      geom_point() +
      guides(color = guide_legend(direction = "horizontal",
                                  override.aes = list(shape = 15, size=9),
                                  ncol  = 4))+
      scale_color_manual(values = cols_cols) + 
      theme(legend.key = element_rect(fill = "transparent"),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 18))
    
    legend <- g_legend(gg) 
    grid.draw(legend)
  }, bg="transparent")
  
  output$key4 <- renderUI({
    pheight = paste0('50px')
    
    if (input$numInputs > 4){
      pheight = paste0('100px')
    }
    if (input$numInputs > 8){
      pheight = paste0('150px')
    }
    
    plotOutput("key4_plot", height = pheight)
  })
  
  output$key5_plot <- renderPlot({
    validate(
      need(!is.null(new_group_names()[[1]]), "")
    )
    #print(is.null(new_group_names()[[1]]))
    groupings = new_group_names()
    cols_cols <- color_select()
    
    n <- length(groupings)
    Group = groupings
    
    plotplot <- data.frame(Group,cols_cols)
    plotplot$Group <- factor(plotplot$Group, levels = groupings)
    
    gg = ggplot(plotplot, aes(1:n, 1:n, color = Group)) +
      geom_point() +
      guides(color = guide_legend(direction = "horizontal",
                                  override.aes = list(shape = 15, size=9),
                                  ncol  = 4))+
      scale_color_manual(values = cols_cols) + 
      theme(legend.key = element_rect(fill = "transparent"),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 18))
    
    legend <- g_legend(gg) 
    grid.draw(legend)
  }, bg="transparent")
  
  output$key5 <- renderUI({
    pheight = paste0('50px')
    
    if (input$numInputs > 4){
      pheight = paste0('100px')
    }
    if (input$numInputs > 8){
      pheight = paste0('150px')
    }
    
    plotOutput("key5_plot", height = pheight)
  })
  
  
  output$key6_plot <- renderPlot({
    validate(
      need(!is.null(new_group_names()[[1]]), "")
    )
    #print(is.null(new_group_names()[[1]]))
    groupings = new_group_names()
    cols_cols <- color_select()
    
    n <- length(groupings)
    Group = groupings
    
    plotplot <- data.frame(Group,cols_cols)
    plotplot$Group <- factor(plotplot$Group, levels = groupings)
    
    gg = ggplot(plotplot, aes(1:n, 1:n, color = Group)) +
      geom_point() +
      guides(color = guide_legend(direction = "horizontal",
                                  override.aes = list(shape = 15, size=9),
                                  ncol  = 4))+
      scale_color_manual(values = cols_cols) + 
      theme(legend.key = element_rect(fill = "transparent"),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 28))
    
    legend <- g_legend(gg) 
    grid.draw(legend)
  }, bg="transparent")
  
  output$key6 <- renderUI({
    pheight = paste0('50px')
    
    if (input$numInputs > 4){
      pheight = paste0('100px')
    }
    if (input$numInputs > 8){
      pheight = paste0('150px')
    }
    
    plotOutput("key6_plot", height = pheight)
  })
  
  
  output$key7_plot <- renderPlot({
    validate(
      need(!is.null(new_group_names()[[1]]), "")
    )
    #print(is.null(new_group_names()[[1]]))
    groupings = new_group_names()
    cols_cols <- color_select()
    
    n <- length(groupings)
    Group = groupings
    
    plotplot <- data.frame(Group,cols_cols)
    plotplot$Group <- factor(plotplot$Group, levels = groupings)
    
    gg = ggplot(plotplot, aes(1:n, 1:n, color = Group)) +
      geom_point() +
      guides(color = guide_legend(direction = "horizontal",
                                  override.aes = list(shape = 15, size=9),
                                  ncol  = 4))+
      scale_color_manual(values = cols_cols) + 
      theme(legend.key = element_rect(fill = "transparent"),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 28))
    
    legend <- g_legend(gg) 
    grid.draw(legend)
  }, bg="transparent")
  
  output$key7 <- renderUI({
    pheight = paste0('50px')
    
    if (input$numInputs > 4){
      pheight = paste0('100px')
    }
    if (input$numInputs > 8){
      pheight = paste0('150px')
    }
    
    plotOutput("key7_plot", height = pheight)
  })
  
  #
  #
  #
  #
  #
  ##
  ## Gene Plots
  ##

  # Gene Explore Heatmap!
  feat_mat <- reactive({
    
    if (input$testme) {
      validate(
        need(length(group_dims())==4, "Please visit the group tab to verify group assignment")
      )
    } else {
      if (input$input_type == "Biobakery"){
        validate(
          need(input$file1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
      if (input$input_type == "EBI"){
        validate(
          need(input$features1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
    }
    
    feat_plot_df = reorder_mat()
    feat_var <- rowVars(as.matrix(feat_plot_df))
    #print(summary(feat_var))
    #summ <- summary(feat_var)
    keep_quant <- quantile(feat_var, input$var_select)
    keep_filt <- which(feat_var > keep_quant)
    print(length(keep_filt))
    feat_var2 <- feat_plot_df[keep_filt,]
    feat_var2
  })
  
  # output$gene_explore <- renderPlotly({
  #   nums <- feat_mat()
  #   
  #   groupings = new_group_names()
  #   grouping_nums = group_dims()
  #   
  #   bin <- c()
  #   for (i in 1:length(groupings)){
  #     use <- rep(groupings[i], grouping_nums[i])
  #     bin <- c(bin, use)
  #   }
  #   colss = bin
  #   cc = length(unique(bin))
  #   cols <- randomColor(length(groupings))
  #   
  #   for (i in 1:cc){
  #     colss = gsub(groupings[i], cols[i], colss)
  #   }
  #   
  #   col_col <- t(matrix(colss))
  #   col_col <- rbind(col_col, col_col)
  #   
  #   heatty = heatmap.3(t(t(nums)), scale = 'row', margins = c(5,5), col = viridis(100),
  #             labCol = "", labRow = "",
  #             Colv = F, ColSideColors = t(col_col))
  #   
  #   #legend(0.22,0.99,      
  #   #       legend = groupings,
  #   #       col = cols, 
  #   #       lty= 1,             
  #   #       lwd = 5,           
  #   #       cex=1,
  #   #       horiz = T
  #   #)
  #   
  #   gg_nums <- nums
  #   
  #   gg_nums_feat <- gg_nums[heatty$rowInd,]
  #   gg_nums_samp <- gg_nums_feat[,heatty$colInd]
  #   
  #   gg_names <- factor(colnames(gg_nums_samp), levels = colnames(gg_nums_samp))
  #   gg_feat <- factor(rownames(gg_nums_samp), levels = rownames(gg_nums_samp))
  #   
  #   x <- sweep(gg_nums_samp, 1L, rowMeans(gg_nums_samp, na.rm = T), check.margin = FALSE)
  #   sx <- apply(x, 1L, sd, na.rm = T)
  #   x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
  #   
  #   ## for taxa plot later
  #   #gene_order= hclust(dist(x))
  #   #plot(gene_order)
  #   ###
  #   
  #   Z_scored_CPM <- unlist(x)
  #   head(Z_scored_CPM)
  #   
  #   
  #   names <- rep(gg_names, each = nrow(gg_nums))
  #   feat <- rep(gg_feat, ncol(gg_nums))
  #   groups <- rep(bin, each = nrow(gg_nums))
  #   
  #   
  #   plotter <- data.frame(Z_scored_CPM, names, groups, feat)
  #   ###
  #   ##
  #   ### make gmain in plotly!!!!
  #   library(plotly)
  #   
  #   ax <- list(
  #     title = "",
  #     zeroline = FALSE,
  #     showline = FALSE,
  #     showticklabels = FALSE,
  #     tickcolor = 'white',
  #     showgrid = FALSE
  #   )
  #   p <- plot_ly(source = "g_exp",
  #                x = plotter$names, y = plotter$feat,
  #                text = plotter$groups,
  #                z = plotter$Z_scored_CPM, type = "heatmap",
  #                hoverinfo = 'y+text',
  #                colorbar = list(x = -0.2, 
  #                                xanchor = 'left',
  #                                tickmode='array',
  #                                tickvals = c(min(plotter$Z_scored_CPM),max(plotter$Z_scored_CPM)),
  #                                ticktext = c("low", "high"))) %>%
  #     layout(yaxis = ax, xaxis = ax)
  #   
  #   ####
  #   p3 <- subplot(colorer(),
  #                 p, nrows = 2, margin = c(0,0,-0.01,0),
  #                 heights = c(0.1, 0.9), shareX = T)
  #   
  #   p3
  #   
  # })
  # 
  
  da_feat_stat <- reactive({
    withProgress({
      #if (input$input_type == "Biobakery"){
        #   if (input$testme){
      #     validate(
      #       need(length(group_dims())==4, "Please visit the group tab to verify group assignment")
      #     )
      #   }
      #   else {validate(
      #     need(input$file1, "Please provide a file in the Upload Tab") %then%
      #       need(unlist(grouped_samps()), "Please select samples in the Group Tab")
      #   )
      #   }
      feat_order <- feat_mat()
      
      groupings = new_group_names()
      grouping_nums = group_dims()
      
      #feat_trans <- clr(feat_order)
      feat_trans <- feat_order
      cc <- data.frame(feat_trans)
      #boxplot.matrix(spec_trans)
      
      RA <- unlist(feat_order)
      RA_clr <- unlist(cc)
      sample_id <- rep(colnames(feat_order), each = nrow(feat_order))
      
      group_titles = c()
      for (i in 1:length(groupings)) {
        title = groupings[i]
        reps = grouping_nums[i]
        title_rep = rep(title, reps)
        group_titles = c(group_titles, title_rep)
      }
      
      group_titles2 = rep(group_titles, each = nrow(feat_order))
      
      Feature <- rep(rownames(cc), ncol(cc))
      
      stat_df <- data.frame(Sample = sample_id, Group = group_titles2, RA = RA,
                            RA_clr = RA_clr, Feature = Feature)
      
      
      #a_test <- anova(lm(RA_clr ~ Group + Taxa, stat_df))
      uniq_feat <- as.character(unlist(unique(stat_df$Feature)))
      stat_results <- data.frame(Feature = NA, F_val = NA, P_val = NA)
      hoc_results <- data.frame(Feat = NA, Int = NA, p_adj = NA)
      for (i in 1:length(uniq_feat)){
        feat <- uniq_feat[i]
        test <- subset(stat_df, Feature == feat)
        a_test <- try(anova(lm(RA ~ Group, test)))
        result <- data.frame(Feature = uniq_feat[i],
                             F_val = a_test$`F value`[1],
                             P_val = a_test$`Pr(>F)`[1])
        stat_results <- rbind(stat_results, result)
        if (is.na(result$P_val)){}
        else{
          if (result$P_val < 0.05){
            h_test <- aov(RA ~ Group, test)
            tk <- TukeyHSD(h_test, "Group")
            hresult <- data.frame(Feat = rep(uniq_feat[i], ncol(combn(length(unique(test$Group)), 2))),
                                  Int = rownames(tk$Group),
                                  p_adj = tk$Group[,4])
            hoc_results <- rbind(hoc_results, hresult)
          }
        }
      }
      stat_results$p_adj <- p.adjust(stat_results$P_val, method = "BH")
      stat2 <- subset(stat_results, Feature %in% unique(hoc_results$Feat))
      
      #stat_results
      stat_results <- hoc_results
      #}
      
      # if (input$input_type == "EBI"){
      #   # validate(
      #   #   need(input$taxa1, "Please provide a Taxa File in the Upload Tab")%then%
      #   #     need(unlist(grouped_samps()), "Please select samples in the Group Tab")
      #   # )
      #   
      #   feat_order <- reorder_mat()
      #   
      #   groupings = new_group_names()
      #   grouping_nums = group_dims()
      #   
      #   meta <- data.frame(sample_id = colnames(feat_order),
      #                      sample = 1:ncol(feat_order))
      #   
      #   bin <- c()
      #   for (i in 1:length(groupings)){
      #     use <- rep(groupings[i], grouping_nums[i])
      #     bin <- c(bin, use)
      #   }
      #   meta$bins <- bin
      #   
      #   rownames(meta) <- meta$sample_id
      #   
      #   ## deseq2
      #   
      #   countData <- feat_order
      #   condition <- factor(meta$bins)
      #   
      #   dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~condition)
      #   gm_mean = function(x, na.rm=TRUE){
      #     exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
      #   }
      #   geoMeans = apply(counts(dds), 1, gm_mean)
      #   diagdds = estimateSizeFactors(dds, geoMeans = geoMeans)
      #   diagdds = DESeq(diagdds, fitType="local")
      #   res<-results(diagdds)
      #   res<-res[order(res$padj),]
      #   cc <- data.frame(res)
      #   
      #   stat_results <- cc
      #   #keepers <- unique(rownames(go_de_sig))
      #   #print(length(keepers))
      #   #keepers
      # }
    }, message = "Calculating Differentially Abundant Features")
    stat_results
  })
  
  da_feat <- reactive({
    stat_results <- da_feat_stat()
    
    # if (input$input_type == "EBI") {
    #   go_de <- subset(stat_results, log2FoldChange > input$feat_FC_up | log2FoldChange < -input$feat_FC_up)
    #   go_de_sig <- subset(go_de, padj < input$feat_pval_up)
    #   keepers <- unique(rownames(go_de_sig))
    #   print(length(keepers))
    #   keepers
    # } else {
    keepers <- unique(subset(stat_results, p_adj < input$feat_pval_up)$Feat)
    print(length(keepers))
    keepers
    # }
    # keepers
  })
  
  gene_da_plotly <- reactive({
    if(input$numInputs > 1){
      da_gene <- da_feat()
      feat_order <- feat_mat()
      go_show <- feat_order[da_gene,]
    } else {
      feat_order <- feat_mat()
      go_show <- feat_order
    }
    
    
    
    groupings = new_group_names()
    grouping_nums = group_dims()
    
    
    bin <- c()
    for (i in 1:length(groupings)){
      use <- rep(groupings[i], grouping_nums[i])
      bin <- c(bin, use)
    }
    
    go_show_nums2 <- log10(go_show + 1)
    
    gg_nums <- go_show_nums2
    
    feat_clust <- hclust(dist(gg_nums))
    
    gg_nums_feat <- gg_nums[feat_clust$order,]
    gg_nums_samp <- gg_nums_feat
    
    gg_names <- factor(colnames(gg_nums_samp), levels = colnames(gg_nums_samp))
    gg_feat <- factor(rownames(gg_nums_samp), levels = rownames(gg_nums_samp))
    
    x <- sweep(gg_nums_samp, 1L, rowMeans(gg_nums_samp, na.rm = T), check.margin = FALSE)
    sx <- apply(x, 1L, sd, na.rm = T)
    x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
    
    ## for taxa plot later
    gene_order= hclust(dist(x))
    #plot(gene_order)
    ###
    
    #Z_scored_CPM = unlist(data.frame(clr(gg_nums_samp)))
    Z_scored_CPM <- unlist(x)
    #print(head(Z_scored_CPM))
    
    
    names <- rep(gg_names, each = nrow(gg_nums))
    feat <- rep(gg_feat, ncol(gg_nums))
    groups <- rep(bin, each = nrow(gg_nums))
    
    plotter <- data.frame(Z_scored_CPM, names, groups, feat)
    
    ###
    ##
    ### make gmain in plotly!!!!
    library(plotly)
    
    ax <- list(
      title = "",
      zeroline = FALSE,
      showline = FALSE,
      showticklabels = FALSE,
      tickcolor = 'white',
      showgrid = FALSE
    )
    p <- plot_ly(source = "g_exp",
                 x = plotter$names, y = plotter$feat,
                 text = plotter$groups,
                 z = plotter$Z_scored_CPM, type = "heatmap",
                 hoverinfo = 'y+text',
                 colorbar = list(x = -0.2, 
                                 xanchor = 'left',
                                 tickmode='array',
                                 tickvals = c(min(plotter$Z_scored_CPM),max(plotter$Z_scored_CPM)),
                                 ticktext = c("low", "high"))) %>%
      layout(yaxis = ax, xaxis = ax)
    
    ####
    p3 <- subplot(colorer(),
                  p, nrows = 2, margin = c(0,0,-0.01,0),
                  heights = c(0.1, 0.9), shareX = T)
    
    p3
    
  })
    
  observeEvent({
    input$feat_pval_up
  },{
    
    output$gene_da <- renderPlotly({
      plott <- gene_da_plotly()
      plott
    })
  })
  
  plot3_df <- reactive({
    event.data <- event_data("plotly_click", source = "g_exp")
    
    # If NULL dont do anything
    #if(is.null(event.data) == T) return(NULL)
    validate(
      need(is.null(event.data) == F, "Click on the Heatmap") %then%
        need(event.data$y %in% full_file_feature()$Feature, "Click on the Heatmap")
    )
    
    select_gfam <- event.data$y
    
    full_full <- full_file_feature()
    full_select <- data.frame(subset(full_full, Feature %in% select_gfam))
    
    go_nums <- full_select[,4:ncol(full_select)]
    go_reorder <- go_nums[,unlist(grouped_samps())]
    
    RA <- unlist(go_reorder)
    
    samp_order <- factor(colnames(go_reorder), 
                         levels = c(colnames(go_reorder)))
    Sample_id <- rep(samp_order, each = nrow(go_reorder))
    
    Taxa <- rep(full_select$Taxa, ncol(go_reorder))
  

    ## groups
    groupings = new_group_names()
    grouping_nums = group_dims()
    bin <- c()
    for (i in 1:length(groupings)){
      use <- rep(groupings[i], grouping_nums[i])
      bin <- c(bin, use)
    }
    Group = rep(bin, each = nrow(go_reorder))
    ###
    
    plotter <- data.frame(Sample_id, Taxa, Group, RA)
    plotter 
  })
  
  plot3_ggplot <- reactive({
    plotter <- plot3_df()
    
    #### filter to speed up plotting!!!!
    plotter2 <- aggregate(plotter$RA, list(plotter$Taxa), sum)
    
    summer <- sum(plotter2[,2])
    plotter2$prop <- plotter2[,2]/summer
    plotter3 <- subset(plotter2, 
                           prop <= quantile(plotter2$prop, input$taxa_limit[2]))
    plotter4 <- subset(plotter3, 
                           prop > quantile(plotter3$prop, input$taxa_limit[1]))
    
    keep_taxa <- as.character(unlist(plotter4[,1]))
    plotter2 <- subset(plotter, Taxa %in% keep_taxa)
    
    ####
    plotter <- plotter2

    bug_sorter <- plotter[order(plotter$RA, decreasing = T),]
    bug_sort <- as.character(unique(bug_sorter$Taxa))
    plotter$Taxa <- factor(plotter$Taxa, levels = bug_sort)
    
    tax_sel <- ggplot(plotter, aes(x=Sample_id, y = RA, fill = Taxa, text = Group)) + 
      stat_summary(fun.y = "mean", geom = "bar", position = "fill") +
      scale_fill_manual(values = randomColor(length(unique(plotter$Taxa))),
                        na.value = "grey") +
      #ggtitle(select_gfam) + 
      xlab("Sample") + ylab("Relative Abundance") +
      theme(panel.background = element_blank(),
            axis.text.x = element_blank())
    tax_sel
    # taxly <- ggplotly(tax_sel, tooltip = c('fill', 'text'))
    # 
    # 
    # p3 <- subplot(colorer(), taxly,
    #               nrows = 2, margin = c(0,0,-0.01,0),
    #               heights = c(0.1, 0.9), shareX = T)
    # p3
  })
  
  output$Plot3 <- renderPlotly({
    event.data <- event_data("plotly_click", source = "g_exp")
    
    # If NULL dont do anything
    #if(is.null(event.data) == T) return(NULL)
    validate(
      need(is.null(event.data) == F, "Click on the Heatmap") %then%
        need(event.data$y %in% full_file_feature()$Feature, "Click on the Heatmap")
    )
    tax_sel <- plot3_ggplot()+theme(legend.position='none')
    
    taxly <- ggplotly(tax_sel, tooltip = c('fill', 'text'))
    
    
    p3 <- subplot(colorer(), taxly,
                  nrows = 2, margin = c(0,0,-0.01,0),
                  heights = c(0.1, 0.9), shareX = T)
    p3
  })

  
  
  ##### DA gene post hoc test!
  
  da_feat_stat_mat <- reactive({
    event.data <- event_data("plotly_click", source = "g_exp")
    
    validate(
      need(is.null(event.data) == F, "Click on the Heatmap")
    )
    
    #if (input$input_type == "Biobakery"){
    validate(
      need(event.data$y %in% da_feat_stat()$Feat, "Click on the Heatmap")
    )
    
    groupings = new_group_names()
    grouping_nums = group_dims()
    
    #### now select bug
    
    select_feat <- event.data$y
    
    statter <- da_feat_stat()
    statter2 <- subset(statter, Feat %in% select_feat)
    
    group_num <- length(groupings)
    group1 <- gsub('.*-', "", statter2$Int[1])
    group_rest <- gsub("-.*", "", statter2$Int[1:(group_num-1)])
    group_names <- c(group1, group_rest)
    
    resm <- matrix(NA, group_num, group_num)
    resm[lower.tri(resm) ] <-round(statter2$p_adj, 4)
    resm <- t(resm)
    resm[lower.tri(resm) ] <-round(statter2$p_adj, 4)
    rownames(resm) <- group_names
    colnames(resm) <- group_names
    print(resm)
    resm
      
    # } else {
    #   dd_obj <- da_feat_stat()
    #   groupings = new_group_names()
    #   
    #   #contrasts
    #   select_feat <- event.data$y
    #   
    #   contrast_results <- data.frame(Int = NULL, P_adj = NULL)
    #   cond_mat <- combn(groupings, 2)
    #   print(cond_mat)
    #   for (i in 1:ncol(cond_mat)){
    #     use <- cond_mat[,i]
    #     contrast <- c("condition", use[1], use[2])
    #     res<-results(dd_obj, contrast = contrast)
    #     cc <- data.frame(res)
    #     cc2 <- subset(cc, rownames(cc) %in% select_feat)
    #     
    #     interaction <- paste0(use[1], "_", use[2])
    #     P_adj <- cc2$padj
    #     results <- data.frame(Int = interaction, P_adj = P_adj)
    #     contrast_results <- rbind(contrast_results, results)
    #   }
    #   
    #   group_num <- length(groupings)
    #   
    #   resm <- matrix(NA, group_num, group_num)
    #   resm[lower.tri(resm) ] <-round(contrast_results$P_adj, 4)
    #   resm <- t(resm)
    #   resm[lower.tri(resm) ] <-round(contrast_results$P_adj, 4)
    #   rownames(resm) <- groupings
    #   colnames(resm) <- groupings
    #   resm
    # }
    # resm
  })
  
  
  observe({
    event.data <- event_data("plotly_click", source = "g_exp")
    
    # If NULL dont do anything
    if (is.null(event.data) == T){
    } else {
      output$da_feat_stat_heat <- renderPlot({
        resm <- da_feat_stat_mat()
        
        resm_lab <- resm
        resm_lab[resm_lab < 0.0001] <- "<0.0001"
        
        
        if (any(resm<0.05, na.rm = T)){
          breaker = seq(0, 0.05, by = 0.0005)
          coler = c(colorRampPalette(c("red", "white"))(n=100))
        } else {
          breaker <- seq(0, 1, by = 0.01)
          coler <- c(colorRampPalette(c("white"))(n=100))
        }
        
        select_feat <- event.data$y
        #new_name <- unlist(strsplit(select_feat, ";"))
        #new_name <- new_name[length(new_name)]
        
        if (input$numInputs == 2){
          resm[1,2] = resm[1,2] + 0.0000001
        }
        
        par(cex.main=0.8)
        heatmap.3(data.matrix(resm),
                  cellnote = resm_lab,
                  notecol="black",
                  #main=select_feat,
                  #key = TRUE,
                  #keysize = 1.0,
                  key.title = NULL,
                  sepcolor="black",
                  breaks = breaker,
                  col = coler,
                  #breaks = seq(0, 0.05, by = 0.0005),
                  #col=c(colorRampPalette(c("red", "white"))(n=100)),
                  dendrogram = 'none',
                  Rowv=F,
                  Colv=F,
                  margins=c(5,5),
                  cexRow=1.2,
                  cexCol=1.2#,
                  #na.color="gray60"
        )
        title(select_feat, line= -2.5)
      })
    }
  })
  
  
  observe({
    event.data <- event_data("plotly_click", source = "g_exp")
    
    # If NULL dont do anything
    if (is.null(event.data) == T){
    } else {
      output$da_feat_stat_tab <- renderTable({
        resm <- da_feat_stat_mat()
        resm_lab <- resm
        resm_lab[resm_lab < 0.0001] <- "<0.0001"
        resm_lab <- cbind(Names = rownames(resm_lab), resm_lab)
        resm_lab
      }, caption = {
        select_taxa <- event.data$y
        return(select_taxa)
      })
    }
  })
  
  observe({
    validate(
      need(input$numInputs > 1, "Need at least 2 groups for ANOVA Test")
    )
      output$da_feat_stat_ui <- renderUI({
        fluidPage(
          plotOutput("da_feat_stat_heat"),
          fluidRow(fluidPage(downloadButton("feature_stat_download", "Download Heatmap"),
                    downloadButton("feature_stat_results", "Download ANOVA Results")))
        )
      })
  })
  
  
  ### download uis
  
  output$gene_explore_download <- downloadHandler(
    filename = function() { paste("feature_explore", '.png', sep='') },
    content = function(file) {
      p3 <- gene_da_plotly()
      
      export(p3, file = file)
      
    })
  
  
  
  
  output$gene_explore_taxa_download <- downloadHandler(
    filename = function() { paste("taxa_explore", '.png', sep='') },
    content = function(file) {
      
      tax_sel <- plot3_ggplot()+theme(legend.position='none')
      
      taxly <- ggplotly(tax_sel, tooltip = c('fill', 'text'))
      
      
      p3 <- subplot(colorer(), taxly,
                    nrows = 2, margin = c(0,0,-0.01,0),
                    heights = c(0.1, 0.9), shareX = T)
      p3
      
      export(p3, file = file)
      
    }
  )
  
  gene_explore_taxa_legend <- reactive({
    tax_sel <- plot3_ggplot()
    legend1 <- g_legend(tax_sel+guides(fill=guide_legend(ncol=3)))
    legend1
  })
  
  output$gene_explore_taxa_legend_download <- downloadHandler(
    filename = function() { paste("gene_explore_taxa_legend", '.png', sep='') },
    content = function(file) {
      species = length(unique(plot3_df()$Taxa))
      png(file, width = 35, height = 0.3*(species/3), units ='cm', res = 300)
      grid.draw(gene_explore_taxa_legend())
      dev.off()
    }
  )
  
  
  
  output$feature_stat_download <- downloadHandler(
    filename = function() { paste("differential_features_heat", '.png', sep='') },
    content = function(file) {
      event.data <- event_data("plotly_click", source = "g_select")
      resm <- da_feat_stat_mat()
      
      resm_lab <- resm
      resm_lab[resm_lab < 0.0001] <- "<0.0001"
      
      
      if (any(resm<0.05, na.rm = T)){
        breaker = seq(0, 0.05, by = 0.0005)
        coler = c(colorRampPalette(c("red", "white"))(n=100))
      } else {
        breaker <- seq(0, 1, by = 0.01)
        coler <- c(colorRampPalette(c("white"))(n=100))
      }
      
      select_gfam <- event.data$y
      
      if (input$numInputs == 2){
        resm[1,2] = resm[1,2] + 0.0000001
      }
      
      png(file, width = 25, height = 20, units = 'cm', res = 300)
      par(cex.main=0.8)
      v = heatmap.3(data.matrix(resm),
                    cellnote = resm_lab,
                    notecol="black",
                    #main=new_name,
                    #key = TRUE,
                    #keysize = 1.0,
                    key.title = NULL,
                    breaks = breaker,
                    col = coler,
                    #breaks = seq(0, 0.05, by = 0.0005),
                    #col=c(colorRampPalette(c("red", "white"))(n=100)),
                    dendrogram = 'none',
                    Rowv=F,
                    Colv=F,
                    margins=c(10,10),
                    cexRow=1.2,
                    cexCol=1.2#,
                    #na.color="gray60"
      )
      v
      title(select_gfam, line= -2.5)
      dev.off()
    }
  )
  
  output$feature_stat_results <- downloadHandler(
    filename = function() { paste("differential_features_anova", '.txt', sep='') },
    content = function(file) {
      stat_results <- da_feat_stat()
      keepers <- subset(stat_results, p_adj < input$feat_pval_up)
      write.table(keepers, file, row.names = FALSE, sep = '\t', quote = FALSE)
    }
  )
  
  
  ####### Query Your Data ######
  
  
  #### new plot1 for selectable features
  
  ####
  
  plot1_df <- reactive({
    validate(
      need(length(input$acc_list) > 0, "")
    )
    
    da_gene <- input$acc_list
    feat_order <- reorder_mat()
    
    groupings = new_group_names()
    grouping_nums = group_dims()
    
    
    bin <- c()
    for (i in 1:length(groupings)){
      use <- rep(groupings[i], grouping_nums[i])
      bin <- c(bin, use)
    }
    
    go_show <- subset(feat_order, rownames(feat_order) %in% da_gene)
    
    go_show_nums2 <- log10(go_show + 1)
    
    x <- sweep(go_show_nums2, 1L, rowMeans(go_show_nums2, na.rm = T), check.margin = FALSE)
    sx <- apply(x, 1L, sd, na.rm = T)
    x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
    
    Z_scored_CPM <- unlist(x)
    

    names <- rep(colnames(go_show_nums2), each = nrow(go_show_nums2))
    feat <- rep(rownames(go_show_nums2), ncol(go_show_nums2))
    groups <- rep(bin, each = nrow(go_show_nums2))
    
    plotter <- data.frame(Z_scored_CPM, names, groups, feat)
    
    plotter$names <- factor(plotter$names, levels = colnames(go_show_nums2))
    
    plotter
    
  })
  
  plot1_plotly <- reactive({
      plotter <- plot1_df()

    ax <- list(
      title = "",
      zeroline = FALSE,
      showline = FALSE,
      showticklabels = FALSE,
      tickcolor = 'white',
      showgrid = FALSE
    )
    p <- plot_ly(source = "g_select",
                 x = plotter$names, y = plotter$feat,
                 text = plotter$groups,
                 z = plotter$Z_scored_CPM, type = "heatmap",
                 hoverinfo = 'y+text',
                 colorbar = list(x = -0.2, 
                                 xanchor = 'left',
                                 tickmode='array',
                                 tickvals = c(min(plotter$Z_scored_CPM),max(plotter$Z_scored_CPM)),
                                 ticktext = c("low", "high"))) %>%
      layout(yaxis = ax, xaxis = ax)
    
    ####
    p3 <- subplot(colorer(),
                  p, nrows = 2, margin = c(0,0,-0.01,0),
                  heights = c(0.1, 0.9), shareX = T)
    
    p3
    
  })
  
  output$plot1 <- renderPlotly({
    plot1_plotly <- plot1_plotly()
    plot1_plotly
  })
  
  
  ##### new spec select plot
  
  spec_select_df <- reactive({
    event.data <- event_data("plotly_click", source = "g_select")
    
    # If NULL dont do anything
    #if(is.null(event.data) == T) return(NULL)
    validate(
      need(is.null(event.data) == F, "Click on the Heatmap") %then%
        need(event.data$y %in% input$acc_list, "Click on the Heatmap")
    )
   
    select_gfam <- event.data$y
    
    full_full <- full_file_feature()
    full_select <- data.frame(subset(full_full, Feature %in% select_gfam))
    
    go_nums <- full_select[,4:ncol(full_select)]
    go_reorder <- go_nums[,unlist(grouped_samps())]
    
    RA <- unlist(go_reorder)
    
    samp_order <- factor(colnames(go_reorder), 
                         levels = c(colnames(go_reorder)))
    Sample_id <- rep(samp_order, each = nrow(go_reorder))
    
    Taxa <- rep(full_select$Taxa, ncol(go_reorder))
    
    
    ## groups
    groupings = new_group_names()
    grouping_nums = group_dims()
    bin <- c()
    for (i in 1:length(groupings)){
      use <- rep(groupings[i], grouping_nums[i])
      bin <- c(bin, use)
    }
    Group = rep(bin, each = nrow(go_reorder))
    ###
    #print(length(Group))
    #print(length(Sample_id))
    #print(length(Taxa))
    #print(length(RA))
    plotter <- data.frame(Sample_id, Taxa, Group, RA)
    
  })
  
  spec_select_ggplot <- reactive({
    
    plotter <- spec_select_df()
    
    #### filter to speed up plotting!!!!
    plotter2 <- aggregate(plotter$RA, list(plotter$Taxa), sum)
    
    summer <- sum(plotter2[,2])
    plotter2$prop <- plotter2[,2]/summer
    plotter3 <- subset(plotter2, 
                       prop <= quantile(plotter2$prop, input$taxa_limit[2]))
    plotter4 <- subset(plotter3, 
                       prop > quantile(plotter3$prop, input$taxa_limit[1]))
    
    keep_taxa <- as.character(unlist(plotter4[,1]))
    plotter2 <- subset(plotter, Taxa %in% keep_taxa)
    
    ####
    plotter <- plotter2
    
    bug_sorter <- plotter[order(plotter$RA, decreasing = T),]
    bug_sort <- as.character(unique(bug_sorter$Taxa))
    plotter$Taxa <- factor(plotter$Taxa, levels = bug_sort)
    
    tax_sel <- ggplot(plotter, aes(x=Sample_id, y = RA, fill = Taxa, text = Group)) + 
      stat_summary(fun.y = "mean", geom = "bar", position = "fill") +
      scale_fill_manual(values = randomColor(length(unique(plotter$Taxa))),
                        na.value = "grey") +
      #ggtitle(select_gfam) + 
      xlab("Sample") + ylab("Relative Abundance") +
      theme(panel.background = element_blank(),
            axis.text.x = element_blank())
    tax_sel
  })
  
  
  output$spec_select <- renderPlotly({
    tax_sel <- spec_select_ggplot() + theme(legend.position='none')
    
    taxly <- ggplotly(tax_sel, tooltip = c('fill', 'text'))
    
    
    p3 <- subplot(colorer(), taxly,
                  nrows = 2, margin = c(0,0,-0.01,0),
                  heights = c(0.1, 0.9), shareX = T)
    p3
  })
  
  output$curr_select_search <- renderText({
    event.data <- event_data("plotly_click", source = "g_select")
    select_gfam <- event.data$y
    message <- paste0("Current Selection from Plot: ", select_gfam)
    message
  })
  
  output$curr_select_exp <- renderText({
    event.data <- event_data("plotly_click", source = "g_exp")
    select_gfam <- event.data$y
    message <- paste0("Current Selection from Plot: ", select_gfam)
    message
  })
  
  output$curr_taxa_exp <- renderText({
    event.data <- event_data("plotly_click", source = "t_exp")
    select_taxa <- event.data$y
    message <- paste0("Current Selection from Plot: ", select_taxa)
    message
  })
  
  
  
  #### select anova and post hoc results
  
  select_stat_results <- reactive({
    event.data <- event_data("plotly_click", source = "g_select")
    
    # If NULL dont do anything
    #if(is.null(event.data) == T) return(NULL)
    validate(
      need(is.null(event.data) == F, "Click on the Heatmap") %then%
        need(event.data$y %in% input$acc_list, "Click on the Heatmap")
    )
    
    select_gfam <- event.data$y
    
    feat_order <- reorder_mat()
    feat_order2 <- subset(feat_order, rownames(feat_order) %in% select_gfam)
    
    groupings = new_group_names()
    grouping_nums = group_dims()
    
    
    
    RA <- unlist(feat_order2)
    sample_id <- rep(colnames(feat_order2), each = nrow(feat_order2))
    
    group_titles = c()
    for (i in 1:length(groupings)) {
      title = groupings[i]
      reps = grouping_nums[i]
      title_rep = rep(title, reps)
      group_titles = c(group_titles, title_rep)
    }
    
    group_titles2 = rep(group_titles, each = nrow(feat_order2))
    
    Feature <- rep(rownames(feat_order2), ncol(feat_order2))
    
    stat_df <- data.frame(Sample = sample_id, Group = group_titles2, 
                          RA = RA, Feature = Feature)
    
    #a_test <- anova(lm(RA_clr ~ Group + Taxa, stat_df))
    hoc_results <- data.frame(Feat = NA, Int = NA, p_adj = NA)
    a_test <- anova(lm(RA ~ Group, stat_df))
    result <- data.frame(Feature = select_gfam,
                         F_val = a_test$`F value`[1],
                         P_val = a_test$`Pr(>F)`[1])
    validate(
      need(result$P_val < 0.05, "This feature was not significant at the set threshold")
    )
    
    h_test <- aov(RA ~ Group, stat_df)
    tk <- TukeyHSD(h_test, "Group")
    hresult <- data.frame(Feat = rep(select_gfam, ncol(combn(length(unique(stat_df$Group)), 2))),
                          Int = rownames(tk$Group),
                          p_adj = tk$Group[,4])
    hoc_results <- hresult
    
    hoc_results
    
  })
  
  select_stat_plot <- reactive({
    hoc_results <- select_stat_results()
    
    event.data <- event_data("plotly_click", source = "g_select")
    select_gfam <- event.data$y
    
    groupings = new_group_names()
    grouping_nums = group_dims()
    
    group_num <- length(groupings)
    group1 <- gsub('.*-', "", hoc_results$Int[1])
    group_rest <- gsub("-.*", "", hoc_results$Int[1:(group_num-1)])
    group_names <- c(group1, group_rest)
    
    resm <- matrix(NA, group_num, group_num)
    resm[lower.tri(resm) ] <-round(hoc_results$p_adj, 4)
    resm <- t(resm)
    resm[lower.tri(resm) ] <-round(hoc_results$p_adj, 4)
    rownames(resm) <- group_names
    colnames(resm) <- group_names
    print(resm)
    resm
    
    resm_lab <- resm
    resm_lab[resm_lab < 0.0001] <- "<0.0001"
    
    
    if (any(resm<0.05, na.rm = T)){
      breaker = seq(0, 0.05, by = 0.0005)
      coler = c(colorRampPalette(c("red", "white"))(n=100))
    } else {
      breaker <- seq(0, 1, by = 0.01)
      coler <- c(colorRampPalette(c("white"))(n=100))
    }
    
    if (input$numInputs == 2){
      resm[1,2] = resm[1,2] + 0.0000001
    }
    
    par(cex.main=0.8)
    
    v = heatmap.3(data.matrix(resm),
              cellnote = resm_lab,
              notecol="black",
              #main=new_name,
              #key = TRUE,
              #keysize = 1.0,
              key.title = NULL,
              breaks = breaker,
              col = coler,
              #breaks = seq(0, 0.05, by = 0.0005),
              #col=c(colorRampPalette(c("red", "white"))(n=100)),
              dendrogram = 'none',
              Rowv=F,
              Colv=F,
              margins=c(10,10),
              cexRow=1.2,
              cexCol=1.2#,
              #na.color="gray60"
    )
    v
    title(select_gfam, line= -2.5)
  })
  
  output$select_stat_heat <- renderPlot({
    select_mat <- select_stat_plot()
    select_mat
  })
  
  
  observe({
    validate(
      need(input$numInputs > 1, "Need at least 2 groups for ANOVA Test")
    )
    output$select_feat_stat_ui <- renderUI({
      fluidPage(
        plotOutput("select_stat_heat"),
        fluidRow(fluidPage(downloadButton("sel_stat_download", "Download Plot"),
                 downloadButton("sel_stat_results", "Download ANOVA Results")))
      )
    })
  })
  
  
  #####
  ### download uis
  
  output$sel_explore_download <- downloadHandler(
    filename = function() { paste("select_explore", '.png', sep='') },
    content = function(file) {
      p3 <- plot1_plotly()
      
      export(p3, file = file)
      
    })
  
  
  
  output$sel_explore_taxa_download <- downloadHandler(
    filename = function() { paste("select_taxa", '.png', sep='') },
    content = function(file) {
      
      tax_sel <- spec_select_ggplot()+theme(legend.position='none')
      
      taxly <- ggplotly(tax_sel, tooltip = c('fill', 'text'))
      
      
      p3 <- subplot(colorer(), taxly,
                    nrows = 2, margin = c(0,0,-0.01,0),
                    heights = c(0.1, 0.9), shareX = T)
      p3
      
      export(p3, file = file)
      
    }
  )
  
  sel_explore_taxa_legend <- reactive({
    tax_sel <- spec_select_ggplot()
    legend1 <- g_legend(tax_sel+guides(fill=guide_legend(ncol=3)))
    legend1
  })
  
  output$sel_explore_taxa_legend_download <- downloadHandler(
    filename = function() { paste("select_taxa_legend", '.png', sep='') },
    content = function(file) {
      species = length(unique(spec_select_df()$Taxa))
      png(file, width = 35, height = 0.3*(species/3), units ='cm', res = 300)
      grid.draw(sel_explore_taxa_legend())
      dev.off()
    }
  )
  
  
  output$sel_stat_download <- downloadHandler(
    filename = function() { paste("select_differential_features_heat", '.png', sep='') },
    content = function(file) {
      sel_ggplot <- select_stat_plot()
      png(file, height = 20, width = 25, units = 'cm', res = 300)
      print(sel_ggplot)
      dev.off()
    }
  )
  
  output$sel_stat_results <- downloadHandler(
    filename = function() { paste("select_differential_features_anova", '.txt', sep='') },
    content = function(file) {
      stat_results <- select_stat_results()
      write.table(stat_results, file, row.names = FALSE, sep = '\t', quote = FALSE)
    }
  )
  
  
  
  # gene_explore_plot = reactive({
  #   if (input$testme) {
  #     validate(
  #       need(length(group_dims())==4, "Please visit the group tab to verify group assignment")
  #     )
  #   } else {
  #     validate(
  #       need(input$file1, "Please provide a file in the Upload Tab") %then%
  #       need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
  #   }
  #   
  #   orig_df = acc_full()
  #   #"Acc", "Gene.Family", "SpeciesNumber", "Species"
  #   
  #   plot_df = reorder_mat()
  #   plot_df$Gene_Family = paste0(" ", orig_df$Gene_Family, " ")
  #   
  #   keep_genes = paste0(" ", top_pca_list(), " ")
  #   
  #   g_data1 = plot_df[grep(paste(keep_genes, collapse="|"), plot_df$Gene_Family, invert=FALSE), ]
  #   ### old code shows top 500 gene families
  #   #plot_df$Sum = rowSums(plot_df[,1:ncol(plot_df)-1])
  #   #plot_df = plot_df[order(-plot_df$Sum),]
  #   
  #   #if (nrow(plot_df) < 500){
  #   #  out_bound <- nrow(plot_df)
  #   #}
  #   #else {
  #   #  out_bound = 500
  #   #}
  #   #plot_df2 <- plot_df[1:out_bound,]
  #   
  #   heyo = g_data1[,1:ncol(g_data1)]
  # 
  #   col = (ncol(heyo))-1
  #   row = nrow(heyo)
  # 
  #   RelExp = data.frame(heyo[1:(col)])
  #   RelExp2 = t(RelExp)
  #   df_RelExp = data.frame(RelExp2)
  #   
  #   datas = c()
  #   for (i in 1:row){
  #     data=as.vector(t(df_RelExp[i]))
  #     datas = c(datas,data)
  #   }
  # 
  #   groupings = new_group_names()
  #   grouping_nums = group_dims()
  #   
  #   group_titles = c()
  #   for (i in 1:length(groupings)) {
  #     title = groupings[i]
  #     reps = grouping_nums[i]
  #     title_rep = rep(title, reps)
  #     group_titles = c(group_titles, title_rep)
  #   }
  #   
  #   group_titles2 = rep(group_titles, row)
  #   
  #   gfams = as.character(unlist(heyo$Gene_Family))
  #   Gfam = rep(gfams, each = col)
  #   Gfam = paste0(" ",Gfam," ")
  #   sample_num = paste(1:(col))
  #   samp=strtoi(sample_num)
  #   isamp = rep(samp,row)
  #   
  #   g_data = data.frame(datas, Gfam, isamp, group_titles2)
  #   colnames(g_data) = c("Relative_Abundance", "Gene_Family", "Sample_num", "Group")
  #   #g_data without zeros
  #   
  #   g_data2 = subset(g_data, Relative_Abundance > 0, select=c(Relative_Abundance,Gene_Family,Sample_num, Group))
  #   
  #   if (length(input$excluder)>0){
  #     exclude_list = input$excluder
  #     g_data2 = g_data2[grep(paste(exclude_list, collapse="|"), g_data2$Gene_Family, invert=TRUE), ]
  #   }
  #   
  #   uniq_gfam_num = length(unique(g_data2$Gene_Family))
  #   
  #   library(RColorBrewer)
  #   qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  #   col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  #   col_vector[10] = "blue2" 
  #   
  #   if (uniq_gfam_num > 74) {
  #     col_vector_edit = sample(col_vector, uniq_gfam_num, TRUE)
  #   } else {
  #     col_vector_edit = col_vector
  #   }
  #   
  #   exp_plot = ggplot(g_data2, aes(x = Group, y = Relative_Abundance, fill = Gene_Family)) +
  #     #geom_bar(position = 'fill', stat = 'identity') +
  #     stat_summary(fun.y = "mean", geom = "bar", position = "fill") +
  #     scale_fill_manual(values = col_vector_edit) +
  #     ylab("Relative Abundance") +
  #     theme(legend.position = "none")
  # })
  # 
  # output$gene_explore <- renderPlotly({
  #   gene_explore_plot()
  # })
  # 
  # output$gene_explore_download <- downloadHandler(
  #   filename = function() { paste("gene_explore", '.png', sep='') },
  #   content = function(file) {
  #     ggsave(file, plot = gene_explore_plot() +
  #              theme(axis.title.y = element_text(size = 22),
  #                    axis.title.x = element_text(size = 22),
  #                    axis.text.y = element_text(size = 22),
  #                    axis.text.x = element_text(size = 22)),
  #            device = 'png', 
  #            width = 30, height = 24, units = "cm")
  #   }
  # )
  # #####
  
  
  
  # 
  # output$instructions <- renderText({
  #   instruction = c("Hover over bars to view Gene Family Name")
  #   instruction
  # })
  # 
  # # Gene Search Plot
  # Gene search dataframe for plot
  # exp_plot <- reactive({
  #   if (input$testme) {
  #     validate(
  #       need(length(group_dims())==4, "Please visit the group tab to verify group assignment")
  #     )
  #   }
  #   else {
  #     validate(
  #       need(input$file1, "Please provide a file in the Upload Tab") %then%
  #         need(unlist(reorder_mat()), "Please select samples in the Group Tab")
  #     )
  #   }
  #   validate(
  #     need(input$acc_list, 'Select at least one Gene Family!')
  #   )
  # 
  #   orig_df = acc_full()
  #   #"Acc", "Gene.Family", "SpeciesNumber", "Species"
  #   plot_df = reorder_mat()
  #   plot_df$Feature = paste(" ", orig_df$Feature, " ", sep="")
  # 
  #   gfam_list = input$acc_list
  #   #heyo = plot_df[grep(paste(gfam_list, collapse="|"), acc_select()), ]
  # 
  #   heyo = subset(plot_df, Feature %in% gfam_list)
  #   
  #   #print(gfam_list)
  #   #print(acc_select()[1:5])
  #   
  #   col = (ncol(heyo))-1
  #   row = nrow(heyo)
  # 
  #   RelExp = data.frame(heyo[1:(col)])
  #   RelExp2 = t(RelExp)
  #   df_RelExp = data.frame(RelExp2)
  # 
  #   datas = c()
  #   for (i in 1:row){
  #     data=as.vector(t(df_RelExp[i]))
  #     datas = c(datas,data)
  #   }
  # 
  #   groupings = new_group_names()
  #   grouping_nums = group_dims()
  # 
  #   group_titles = c()
  #   for (i in 1:length(groupings)) {
  #     title = groupings[i]
  #     reps = grouping_nums[i]
  #     title_rep = rep(title, reps)
  #     group_titles = c(group_titles, title_rep)
  #   }
  # 
  #   group_titles2 = rep(group_titles, row)
  # 
  #   gfams = as.vector(t(heyo[1 + col]))
  #   Gfam = rep(gfams, each = col)
  #   sample_num = paste(1:(col))
  #   samp=strtoi(sample_num)
  #   isamp = rep(samp,row)
  #   g_data = data.frame(datas, Gfam, isamp, group_titles2)
  #   colnames(g_data) = c("Exp", "Gfam", "Sample_num", "groups")
  #   #g_data without zeros
  # 
  # 
  #   g_data2 = subset(g_data, select=c(Exp,Gfam,Sample_num, groups))
  #   g_data2
  # })
  # 
  # # Gene search plot object
  # expression_plot <- reactive({
  #   g_data2 = data.frame(exp_plot())
  #   class(g_data2$Exp) = "numeric"
  # 
  #   g_data2$Gfam2 <- gsub(" ", "", g_data2$Gfam)
  #   #g_data2$Gfam2 <- str_trunc(g_data2$Gfam2, 30, "center")
  # 
  #   g_data2$group_gfam <- paste(g_data2$group, g_data2$Gfam2, sep=":")
  # 
  #   # necessary to keep boxplot boxes at correct width, psuedo data at 1,000,000
  #   combos <- unique(g_data2$group_gfam)
  #   g_data_new <- data.frame()
  #   for (i in 1:length(combos)){
  #     checker <- subset(g_data2, group_gfam == combos[i])
  #     if (sum(checker$Exp) > 0){
  #       new_df <- checker
  #     } else {
  #       new_df <- checker[1,]
  #       new_df$Exp = 1000000
  #     }
  #     g_data_new <- rbind(g_data_new, new_df)
  #   }
  #   colnames(g_data_new) = colnames(g_data2)
  # 
  #   g_data_new <- subset(g_data_new, Exp > 0)
  #   glogs <- g_data_new$Exp[g_data_new$Exp < 1000000]
  #   gmax <- max(glogs)
  #   gmin <- min(glogs)
  # 
  #   if (input$xy_switch){
  #     plot_exp = ggplot(g_data_new, aes(x = Gfam2, y = Exp, fill = groups)) +
  #       geom_boxplot(outlier.shape=3) +
  #       ggtitle("Logarthimic Gene Abundance") +
  #       coord_cartesian(ylim = c(gmin, gmax)) +
  #       scale_y_log10() +
  #       xlab("Gene Family") +
  #       ylab("Log Relative Abundance") +
  #       guides(color=FALSE, fill = guide_legend(title = "Group")) +
  #       theme(legend.position = 'none',
  #         #axis.text.x = element_blank(),
  #         #plot.margin = unit(c(.1, .1, .1, .1), "cm"),
  #         plot.title = element_text(hjust = 0, size = 22),
  #         axis.title.y = element_text(size = 22),
  #         axis.title.x = element_text(size = 18),
  #         axis.text.y = element_text(size = 22),
  #         axis.text.x = element_text(size = 18,
  #                                    angle = 20+5*(length(unique(g_data_new$Gfam2))),
  #                                    hjust=1),
  #         legend.title = element_text(size = 22),
  #         legend.text = element_text(size = 18))
  #     plot_exp
  #   } else {
  #     plot_exp = ggplot(g_data_new, aes(x = groups, y = Exp)) +
  #       geom_boxplot(outlier.shape=3, aes(fill = Gfam2)) +
  #       geom_point(position=position_dodge(width=0.75), 
  #                  aes(group=Gfam2),
  #                  color = "grey20") +
  #       ggtitle("Logarthimic Gene Abundance") +
  #       coord_cartesian(ylim = c(gmin, gmax)) +
  #       scale_y_log10() +
  #       xlab("Groups") +
  #       ylab("Log Relative Abundance") +
  #       guides(color=TRUE, fill = guide_legend(title = "Feature",
  #                                               ncol = 1)) +
  #       scale_fill_manual(values = randomColor(length(unique(g_data_new$Gfam2)))) +
  #       theme(legend.position = 'bottom',
  #             #panel.background = element_blank(),
  #             #panel.grid.minor = theme_blank(), 
  #             #panel.grid.major = theme_blank(),
  #             #plot.background = element_rect(fill = "transparent",colour = NA),
  #             #axis.text.x = element_blank(),
  #             plot.margin = unit(c(1, 1, 1, 1), "cm"),
  #             plot.title = element_text(hjust = 0, size = 22),
  #             axis.title.y = element_text(size = 22),
  #             axis.title.x = element_text(size = 22),
  #             axis.text = element_text(size =22, color = 'black'),
  #             legend.title = element_text(size = 22),
  #             legend.text = element_text(size = 12))
  #     plot_exp
  #   }
  # })
  # 
  # output$plot1 <- renderPlot({
  #   expression_plot()
  # })

  output$expression_download <- downloadHandler(
    filename = function() { paste("gene_family_abundance", '.png', sep='') },
    content = function(file) {
      ggsave(file, plot = expression_plot(), device = 'png',
             width = 35, height = 20, units = "cm")
    }
  )
  


  # Significance Table for Gene search t-tests
  #
  #change text annoations to "<0.0001"
  # expression_table <- reactive({
  #   g_data2 = data.frame(exp_plot())
  #   g_data2 = g_data2[order(g_data2$groups),]
  #   
  #   Gfam_uniq = unique(g_data2$Gfam)
  #   group_uniq = as.character(unique(g_data2$groups))
  #   if (length(group_uniq) > 1) {
  #     hio=c()
  #     for (i in Gfam_uniq){
  #       #acc_split = g_data2[grep(i, g_data2$Gfam),]
  #       acc_split <- subset(g_data2, Gfam == i)
  #       acc_split[c("Exp")][is.na(acc_split[c("Exp")])] <- 0
  #       hi = pairwise.t.test(acc_split$Exp, acc_split$groups, p.adjust = 'BH')
  #       hiP <- hi$p.value
  #       nas <- rep(NA, length(group_uniq))
  #       hiP <- rbind(nas[-1], hiP)
  #       hiP <- cbind(hiP, nas)
  #       hiP[upper.tri(hiP)] = t(hiP)[upper.tri(hiP)]
  #       hiP <- round(hiP, 4)
  #       hiP[hiP < 0.0001] <- "<0.0001"
  #       hiP <- data.frame(hiP)
  #       hio = rbind.fill(hio, hiP)
  #     }
  #     hi3 = data.frame(rep(Gfam_uniq, each = length(group_uniq)))
  #     colnames(hi3) = "Gfam"
  #     hi3$group_comparisons = c(rep(group_uniq, length(Gfam_uniq)))
  #     hi3 = cbind.fill(hi3, hio)
  #     colnames(hi3) <- c("Gfam", "Group Comparison", group_uniq)
  #   }
  #   else if (length(group_uniq)==1){
  #     hi3 <- data.frame(Seletion=1:length(Gfam_uniq), Gene_Family=Gfam_uniq)
  #     hi3
  #   }
  #   hi3
  # })
  # 
  # #actual numbers for correct heat map colors
  # expression_table_orig <- reactive({
  #   g_data2 = data.frame(exp_plot())
  #   g_data2 = g_data2[order(g_data2$groups),]
  # 
  #   Gfam_uniq = unique(g_data2$Gfam)
  #   group_uniq = as.character(unique(g_data2$groups))
  #   if (length(group_uniq) > 1) {
  #     hio=c()
  #     for (i in Gfam_uniq){
  #       #acc_split = g_data2[grep(i, g_data2$Gfam),]
  #       acc_split <- subset(g_data2, Gfam == i)
  #       acc_split[c("Exp")][is.na(acc_split[c("Exp")])] <- 0
  #       hi = pairwise.t.test(acc_split$Exp, acc_split$groups, p.adjust = 'BH')
  #       hiP <- hi$p.value
  #       nas <- rep(NA, length(group_uniq))
  #       hiP <- rbind(nas[-1], hiP)
  #       hiP <- cbind(hiP, nas)
  #       hiP[upper.tri(hiP)] = t(hiP)[upper.tri(hiP)]
  #       hiP <- round(hiP, 5)
  #       hiP <- data.frame(hiP)
  #       hio = rbind.fill(hio, hiP)
  #     }
  #     hi3 = data.frame(rep(Gfam_uniq, each = length(group_uniq)))
  #     colnames(hi3) = "Gfam"
  #     hi3$group_comparisons = c(rep(group_uniq, length(Gfam_uniq)))
  #     hi3 = cbind.fill(hi3, hio)
  #     colnames(hi3) <- c("Gfam", "Group Comparison", group_uniq)
  #   }
  #   else if (length(group_uniq)==1){
  #     hi3 <- data.frame(Seletion=1:length(Gfam_uniq), Gene_Family=Gfam_uniq)
  #     hi3
  #   }
  #   hi3
  # })
  # 
  # output$exp_table = renderTable({
  #   expression_table()
  # })
  # 
  # output$expression_table_download <- downloadHandler(
  #   filename = function() { paste("gene_family_abundance_table", '.txt', sep='') },
  #   content = function(file) {
  #     exp = expression_table()
  #     exp2 = data.frame(exp)
  #     write.table(exp2, file, row.names = FALSE, sep = '\t', quote = FALSE)
  #   }
  # )
  # 
  # output$exp_heat <- renderUI({
  #   g_data2 = data.frame(exp_plot())
  #   g_data2 = g_data2[order(g_data2$groups),]
  # 
  #   Gfam_uniq = unique(g_data2$Gfam)
  #   group_uniq = as.character(unique(g_data2$groups))
  #   exp_table <- expression_table_orig()
  #   exp_labs <- expression_table()
  #   if (length(group_uniq) > 1) {
  #     get_exp_heat_list(max_plots, length(Gfam_uniq), exp_table, exp_labs)
  #   }
  # })
  # 
  # output$exp_heat_download <- renderUI({
  #   g_data2 = data.frame(exp_plot())
  #   g_data2 = g_data2[order(g_data2$groups),]
  # 
  #   Gfam_uniq = unique(g_data2$Gfam)
  #   group_uniq = as.character(unique(g_data2$groups))
  #   if (length(group_uniq) > 1) {
  #     lapply(1:length(Gfam_uniq), function(i) {
  #       display_name = Gfam_uniq
  #       downloadButton(paste0("downloadExp", i), paste("Download", display_name[i], sep=" "))
  #     })
  #   }
  # })
  # 
  # #download gene search t-test heat maps with variable width according to Gene name
  # observe({
  #   g_data2 = data.frame(exp_plot())
  #   g_data2 = g_data2[order(g_data2$groups),]
  # 
  #   Gfam_uniq = unique(g_data2$Gfam)
  #   group_uniq = as.character(unique(g_data2$groups))
  # 
  # 
  #   lapply(1:length(Gfam_uniq), function(i) {
  #     if (nchar(as.character(Gfam_uniq[i])) > 30){
  #       fixed_width <- 0.4*nchar(as.character(Gfam_uniq[i]))
  #     } else {fixed_width <- 15}
  # 
  #     if (nchar(as.character(Gfam_uniq[i])) > 30){
  #       fixed_height <- 0.35*nchar(as.character(Gfam_uniq[i]))
  #     } else {fixed_height <- 15}
  # 
  #     output[[paste0("downloadExp", i)]] <- downloadHandler(
  #       filename = function() { paste(Gfam_uniq[i], "_p_value_heat", '.png', sep='') },
  #       content = function(file) {
  #         png(file, width = fixed_width,
  #             height = fixed_height, units ='cm', res = 300)
  #         exp_table <- expression_table_orig()
  #         exp_labs <- expression_table()
  #         download_exp_heat_list(max_plots, i, exp_table, exp_labs)
  #         dev.off()
  #       }
  #     )
  #   })
  # })
   
  ## Plots for Taxa! ##
  
  ## Taxa instructions
  
  output$ex_delimiter <- renderUI({
    str1 = c("Select Taxa Level of Interest")
    str2 = c("Biobakery: 1 = Genus, 2 = Species")
    str3 = c("EBI: 1= Kingdom, 2=Phylum, 3=Class...7=Species")
    HTML(paste(str1, str2, str3, '<br/>', sep = '<br/>'))
  })
  
  
  ## UI based on group nums
  
  output$da_taxa_ui <- renderUI({
    if (input$numInputs > 1){
      fluidPage(
        fluidRow(uiOutput("taxa_selectors")),
        fluidRow(textOutput("curr_taxa_exp")),
        fluidRow(column(12, uiOutput("key2"))),
        fluidRow(column(12, uiOutput("da_taxa_heat_UI"))),
        fluidRow(downloadButton("species_heat_download", "Download Heatmap")),
        fluidRow(column(12, uiOutput("da_taxa_stat_ui"))),
        fluidRow(downloadButton("species_stat_download", "Download Heatmap"),
                 downloadButton("species_stat_results", "Download ANOVA Results"))
      )
    }
  })
  
  
  # Taxa explore #

  # matrix reordered by group and collapsed by species
  reorder_spec_mat <- reactive({
    req(grouped_samps())
    exprs_reorder = spec_nums()[,c(unlist(grouped_samps()))]
    exp2 = data.frame(exprs_reorder)
    rownames(exp2) <- spec_full()$Taxa
    return(exp2)
  })

  # Explore Taxa dataframe
  spec_taxa_data <- reactive({
    # validate(
    #   need(input$taxaSep, "Please provide a Character Delimiter")
    # )
    if (input$testme) {
      validate(
        need(length(group_dims())==4, "Please visit the group tab to verify group assignment")
      )
    } else {
      if (input$input_type == "Biobakery"){
        validate(
          need(input$file1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
      if (input$input_type == "EBI"){
        validate(
          need(input$features1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
    }
    
    withProgress({

      spec_plot_df = reorder_spec_mat()
      specs = spec_plot_df
      
      specs1 <- cbind(rownames(specs), specs)
      colnames(specs1)[1] <- "Taxa"
      
      if (input$input_type == "EBI"){
        bugg_fix <- gsub(";k__", "_k__", specs1$Taxa)
        bug_levels <- strsplit(bugg_fix, ";")
        
        tax_lev <- c()
        for (i in 1:length(bug_levels)){
          bug_len <- length(bug_levels[[i]])
          if (bug_len == input$taxaDims){
            tax_lev <- c(tax_lev, i)
          }
        }
        specs <- specs1[tax_lev,]
        
        Taxa2 <- as.character(unlist(rownames(specs)))
        Taxa2 <- strsplit(Taxa2, ";")
        Taxa3 <- c()
        for (i in 1:length(Taxa2)){
          bugg <- Taxa2[i]
          keep <- length(bugg[[1]])
          Taxa3 <- c(Taxa3, bugg[[1]][keep])
        }
        specs$Taxa <- Taxa3
        
      }
      if (input$input_type == "Biobakery"){
        if (input$taxaDims == 1){
          specs1$Taxa <- gsub("\\..*", "", specs1$Taxa)
        }
        specs <- specs1
      }
      
      
      spec_col = (ncol(specs))-1
      spec_row = nrow(specs)
      
      spec_RelExp = data.frame(specs[,-1])
      spec_RelExp2 = t(spec_RelExp)
      df_spec_RelExp = data.frame(spec_RelExp2)
      
      spec_datas = c()
      for (i in 1:spec_row) {
        spec_data=as.vector(t(df_spec_RelExp[i]))
        spec_datas = c(spec_datas,spec_data)
      }
      
      groupings = new_group_names()
      grouping_nums = group_dims()
      
      group_titles = c()
      for (i in 1:length(groupings)) {
        title = groupings[i]
        reps = grouping_nums[i]
        title_rep = rep(title, reps)
        group_titles = c(group_titles, title_rep)
      }
      
      group_titles2 = rep(group_titles, spec_row)
      
      samp_order <- factor(colnames(spec_RelExp), 
                           levels = c(colnames(spec_RelExp)))
      sample_id <- rep(samp_order, nrow(spec_RelExp))
      
      #spec_sample_num = paste(1:(spec_col))
      #spec_samp=strtoi(spec_sample_num)
      #spec_isamp = rep(spec_samp,spec_row)
      
      spec_spec <- rep(specs$Taxa, each = spec_col)
      
      spec_g_data = data.frame(spec_datas,
                               sample_id, group_titles2, spec_spec)
      colnames(spec_g_data) = c("Relative_Abundance", "Sample_num", "Group", "Taxa")
      
      spec_g_data
    }, message = "Collecting Taxanomic Data")
  })

  # Explore Taxa Species Plot Object
  species_explore_plot <- reactive({
    withProgress({
      spec_g_data <- spec_taxa_data()
      uniq_spec_num = length(unique(spec_g_data$Taxa))
      
      spec_g_data$Taxa <- reorder(spec_g_data$Taxa, -spec_g_data$Relative_Abundance)
      
      spec_g_data2 <- aggregate(spec_g_data$Relative_Abundance, list(spec_g_data$Taxa), sum)
      
      summer <- sum(spec_g_data2[,2])
      spec_g_data2$prop <- spec_g_data2[,2]/summer
      spec_g_data3 <- subset(spec_g_data2, 
                             prop <= quantile(spec_g_data2$prop, input$taxa_limit[2]))
      spec_g_data4 <- subset(spec_g_data3, 
                             prop > quantile(spec_g_data2$prop, input$taxa_limit[1]))
      
      keep_taxa <- as.character(unlist(spec_g_data4[,1]))
      spec_g_data_filt <- subset(spec_g_data, Taxa %in% keep_taxa)
      
      
      exp_plot = ggplot(spec_g_data_filt, aes(x = Sample_num,
                                              y = Relative_Abundance,
                                              fill = Taxa, text = Group)) +
        stat_summary(fun.y = "mean", geom = "bar", position = "fill")+
        scale_fill_manual(values = randomColor(uniq_spec_num)) +
        ylab("Relative Abundance") + xlab("Sample")
    }, message = "Organizing Taxa")
    exp_plot
  })

  output$species_explore <- renderPlotly({
    withProgress({
      taxa_all <- species_explore_plot() + 
        theme(legend.position = 'none',
              panel.background = element_blank(),
              axis.text.x = element_blank())#legend.text = element_text(size = 8))
      
      taxly <- ggplotly(taxa_all, tooltip = c("fill", "text"))
      
      p3 <- subplot(colorer(), taxly,
                    nrows = 2, margin = c(0,0,-0.01,0),
                    heights = c(0.1, 0.9), shareX = T)
      p3
    }, message = "Rendering")
  })

  species_explore_legend <- reactive({
    legend1 <- g_legend(species_explore_plot()+guides(fill=guide_legend(ncol=3)))
    legend1
  })

  output$species_download <- downloadHandler(
    filename = function() { paste("taxa_explore", '.png', sep='') },
    content = function(file) {
      species = length(unique(spec_taxa_data()$Taxa))
      # ggsave(file, plot = species_explore_plot() +
      #          theme(legend.position = "none",
      #                axis.title.y = element_text(size = 22),
      #                axis.title.x = element_text(size = 22),
      #                axis.text.y = element_text(size = 22),
      #                axis.text.x = element_text(size = 22)),
      #        device = 'png',
      #        width = 30, height = 24, units = "cm")
      taxa_all <- species_explore_plot() + 
        theme(legend.position = 'none',
              panel.background = element_blank(),
              axis.text.x = element_blank())#legend.text = element_text(size = 8))
      
      taxly <- ggplotly(taxa_all, tooltip = c("fill", "text"))
      
      p3 <- subplot(colorer(), taxly,
                    nrows = 2, margin = c(0,0,-0.01,0),
                    heights = c(0.1, 0.9), shareX = T)
      
      export(p3, file = file)
      
    }
  )

  output$species_legend_download <- downloadHandler(
    filename = function() { paste("taxa_explore_legend", '.png', sep='') },
    content = function(file) {
      species = length(unique(spec_taxa_data()$Taxa))
      png(file, width = 35, height = 0.3*(species/3), units ='cm', res = 300)
      grid.draw(species_explore_legend())
      dev.off()
    }
  )

  output$species_raw_data <- downloadHandler(
    filename = function() { paste("taxa_explore_raw_data", '.txt', sep='') },
    content = function(file) {
      spec_g_data <- spec_taxa_data()
      uniq_spec_num = length(unique(spec_g_data$Taxa))

      spec_g_data$Taxa <- reorder(spec_g_data$Taxa, -spec_g_data$Relative_Abundance)
      write.table(spec_g_data, file, quote=F, sep='\t', row.names = F)
    }
  )

  observe({
    #if (input$input_type == "Biobakery"){
    output$taxa_selectors <- renderUI({
      sliderInput("taxa_pval_up", "Adjusted P Value", min = 0, max = 0.1,
                  step = 0.01, value = 0.05, width = '35%')
    })
    output$feat_selectors <- renderUI({
      fluidRow(
        column(4, sliderInput("feat_pval_up", "Adjusted P Value", 
                              min = 0, max = 0.1,
                              step = 0.01, value = 0.05)),
        column(4, sliderInput("var_select", "Variance Filter",
                              min = 0, max = 1,
                              value = 0.75))
      )
    })
    #}
    # if (input$input_type == "EBI") {
    #   output$taxa_selectors <- renderUI({
    #     fluidRow(column(3, 
    #                     numericInput("taxa_FC_up", "Fold Change Level", min = 0, max = 10,
    #                                  step = 0.5, value = 1.5)),
    #              column(3, 
    #                     numericInput("taxa_pval_up", "Adjusted P Value", min = 0, max = 0.1,
    #                                  step = 0.01, value = 0.05)))
    #   })
    #   output$feat_selectors <- renderUI({
    #     fluidRow(column(3, 
    #                     numericInput("feat_FC_up", "Fold Change Level", min = 0, max = 10,
    #                                  step = 0.5, value = 1.5)),
    #              column(3, 
    #                     numericInput("feat_pval_up", "Adjusted P Value", min = 0, max = 0.1,
    #                                  step = 0.01, value = 0.05)))
    #   })
    # }
  })
  
 ##### Test Differential Abundance for Taxa
  
  da_taxa_stat <- reactive({
    withProgress({
      if (input$input_type == "Biobakery"){
        if (input$testme){
          validate(
            need(length(group_dims())==4, "Please visit the group tab to verify group assignment")
          )
        }
        else {validate(
          need(input$file1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab")
        )
        }
      } else {
        validate(
          need(input$taxa1, "Please provide a Taxa File in the Upload Tab")%then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab")
        )
      }
      spec_order <- reorder_spec_mat()
      
      groupings = new_group_names()
      grouping_nums = group_dims()
      
      #spec_trans <- clr(spec_order)
      spec_trans <- spec_order
      cc <- data.frame(spec_trans)
      #boxplot.matrix(spec_trans)
      
      x <- sweep(spec_order, 1L, rowMeans(spec_order, na.rm = T), check.margin = FALSE)
      sx <- apply(x, 1L, sd, na.rm = T)
      x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
      
      
      RA <- unlist(x)
      RA_clr <- unlist(cc)
      sample_id <- rep(colnames(spec_order), each = nrow(spec_order))
      
      group_titles = c()
      for (i in 1:length(groupings)) {
        title = groupings[i]
        reps = grouping_nums[i]
        title_rep = rep(title, reps)
        group_titles = c(group_titles, title_rep)
      }
      
      group_titles2 = rep(group_titles, each = nrow(spec_order))
      
      Taxa <- rep(rownames(spec_order), ncol(spec_order))
      
      stat_df <- data.frame(Sample = sample_id, Group = group_titles2, RA = RA,
                            RA_clr = RA_clr, Taxa = Taxa)
      
      
      #a_test <- anova(lm(RA_clr ~ Group + Taxa, stat_df))
      uniq_bugs <- as.character(unlist(unique(stat_df$Taxa)))
      stat_results <- data.frame(Taxa = NA, F_val = NA, P_val = NA)
      hoc_results <- data.frame(Taxa = NA, Int = NA, p_adj = NA)
      for (i in 1:length(uniq_bugs)){
        taxa <- uniq_bugs[i]
        test <- subset(stat_df, Taxa == taxa)
        a_test <- try(anova(lm(RA ~ Group, test)))
        if (class(a_test)[1] != "try-error"){
          result <- data.frame(Taxa = uniq_bugs[i],
                               F_val = a_test$`F value`[1],
                               P_val = a_test$`Pr(>F)`[1])
          stat_results <- rbind(stat_results, result)
          if (result$P_val < 0.05){
            h_test <- aov(RA ~ Group, test)
            tk <- TukeyHSD(h_test, "Group")
            hresult <- data.frame(Taxa = rep(uniq_bugs[i], ncol(combn(length(unique(test$Group)), 2))),
                                  Int = rownames(tk$Group),
                                  p_adj = tk$Group[,4])
            hoc_results <- rbind(hoc_results, hresult)
          }
        }
      }
      stat_results <- hoc_results
      #}
      
      # if (input$input_type == "EBI"){
      #   validate(
      #     need(input$taxa1, "Please provide a Taxa File in the Upload Tab")%then%
      #       need(unlist(grouped_samps()), "Please select samples in the Group Tab")
      #   )
      #   
      #   spec_order <- reorder_spec_mat()
      #   
      #   groupings = new_group_names()
      #   grouping_nums = group_dims()
      #   
      #   meta <- data.frame(sample_id = colnames(spec_order),
      #                      sample = 1:ncol(spec_order))
      #   
      #   bin <- c()
      #   for (i in 1:length(groupings)){
      #     use <- rep(groupings[i], grouping_nums[i])
      #     bin <- c(bin, use)
      #   }
      #   meta$bins <- bin
      #   
      #   rownames(meta) <- meta$sample_id
      #   
      #   ## deseq2
      #   
      #   countData <- spec_order
      #   condition <- factor(meta$bins)
      #   
      #   dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~condition)
      #   gm_mean = function(x, na.rm=TRUE){
      #     exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
      #   }
      #   geoMeans = apply(counts(dds), 1, gm_mean)
      #   diagdds = estimateSizeFactors(dds, geoMeans = geoMeans)
      #   diagdds = DESeq(diagdds, fitType="local")
      #   #res<-results(diagdds)
      #   #res<-res[order(res$padj),]
      #   #cc <- data.frame(res)
      #   
      #   #stat_results <- cc
      #   stat_results <- diagdds
      # }
      stat_results
    }, message = "Calculating Differentially Abundant Taxa")
  })
  
  da_taxa <- reactive({
    # if (input$input_type == "EBI"){
    #   diagdds <- da_taxa_stat()
    #   res<-results(diagdds)
    #   res<-res[order(res$padj),]
    #   cc <- data.frame(res)
    #   
    #   stat_results <- cc
    #   go_de <- subset(stat_results, log2FoldChange > input$taxa_FC_up | log2FoldChange < -input$taxa_FC_up)
    #   go_de_sig <- subset(go_de, padj < input$taxa_pval_up)
    #   keepers <- unique(rownames(go_de_sig))
    #   print(length(keepers))
    #   keepers
    # } else {
    stat_results <- da_taxa_stat()
    keepers <- unique(subset(stat_results, p_adj < input$taxa_pval_up)$Taxa)
    print(length(keepers))
    keepers
    # }
    # keepers
  })
  
  da_taxa_heat_plotly <- reactive({
    da_taxa <- da_taxa()
    spec_order <- reorder_spec_mat()
    
    groupings = new_group_names()
    grouping_nums = group_dims()
    
    ### colors
    bin <- c()
    for (i in 1:length(groupings)){
      use <- rep(groupings[i], grouping_nums[i])
      bin <- c(bin, use)
    }
    
    go_show <- spec_order[da_taxa,]

    go_show_nums2 <- log10(go_show + 1)
    
    gg_nums <- go_show_nums2
    
    feat_clust <- hclust(dist(gg_nums))
    
    gg_nums_feat <- gg_nums[feat_clust$order,]
    gg_nums_samp <- gg_nums_feat
    
    gg_names <- factor(colnames(gg_nums_samp), levels = colnames(gg_nums_samp))
    gg_feat <- factor(rownames(gg_nums_samp), levels = rownames(gg_nums_samp))
    
    x <- sweep(gg_nums_samp, 1L, rowMeans(gg_nums_samp, na.rm = T), check.margin = FALSE)
    sx <- apply(x, 1L, sd, na.rm = T)
    x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
    
    ## for taxa plot later
    #gene_order= hclust(dist(x))
    #plot(gene_order)
    ###
    
    Z_scored_CPM <- unlist(x)
    head(Z_scored_CPM)
    
    
    names <- rep(gg_names, each = nrow(gg_nums))
    feat <- rep(gg_feat, ncol(gg_nums))
    groups <- rep(bin, each = nrow(gg_nums))
    
    
    plotter <- data.frame(Z_scored_CPM, names, groups, feat)
    ###
    ##
    ### make gmain in plotly!!!!
    library(plotly)
    
    ax <- list(
      title = "",
      zeroline = FALSE,
      showline = FALSE,
      showticklabels = FALSE,
      tickcolor = 'white',
      showgrid = FALSE
    )
    p <- plot_ly(source = "t_exp",
                 x = plotter$names, y = plotter$feat,
                 text = plotter$groups,
                 z = plotter$Z_scored_CPM, type = "heatmap",
                 hoverinfo = 'y+text',
                 colorbar = list(x = -0.2, 
                                 xanchor = 'left',
                                 tickmode='array',
                                 tickvals = c(min(plotter$Z_scored_CPM),max(plotter$Z_scored_CPM)),
                                 ticktext = c("low", "high"))) %>%
      layout(yaxis = ax, xaxis = ax)
    
    
    ####
    p3 <- subplot(colorer(),
                  p, nrows = 2, margin = c(0,0,-0.01,0),
                  heights = c(0.1, 0.9), shareX = T)
    
    p3
    
  })
  
  output$da_taxa_heat <- renderPlotly({
    p3 <- da_taxa_heat_plotly()
    p3
  })
  
  output$da_taxa_heat_UI <- renderUI({
    plotlyOutput("da_taxa_heat", height = 450)
  })
  
  da_taxa_stat_mat <- reactive({
    event.data <- event_data("plotly_click", source = "t_exp")
    
    validate(
      need(is.null(event.data) == F, "Click on the Heatmap")
    )
    
    #if (input$input_type == "Biobakery"){
    validate(
      need(event.data$y %in% da_taxa_stat()$Taxa, "Click on the Heatmap")
    )
    
    groupings = new_group_names()
    grouping_nums = group_dims()
    
    #### now select bug
    
    select_taxa <- event.data$y
    
    statter <- da_taxa_stat()
    statter2 <- subset(statter, Taxa %in% select_taxa)
    
    group_num <- length(groupings)
    group1 <- gsub('.*-', "", statter2$Int[1])
    group_rest <- gsub("-.*", "", statter2$Int[1:(group_num-1)])
    group_names <- c(group1, group_rest)
    
    resm <- matrix(NA, group_num, group_num)
    resm[lower.tri(resm) ] <-round(statter2$p_adj, 4)
    resm <- t(resm)
    resm[lower.tri(resm) ] <-round(statter2$p_adj, 4)
    rownames(resm) <- group_names
    colnames(resm) <- group_names
    print(resm)
    resm
    
    # } else {
    #   dd_obj <- da_taxa_stat()
    #   groupings = new_group_names()
    #   
    #   #contrasts
    #   select_taxa <- event.data$y
    #   
    #   contrast_results <- data.frame(Int = NULL, P_adj = NULL)
    #   cond_mat <- combn(groupings, 2)
    #   print(cond_mat)
    #   for (i in 1:ncol(cond_mat)){
    #     use <- cond_mat[,i]
    #     contrast <- c("condition", use[1], use[2])
    #     res<-results(dd_obj, contrast = contrast)
    #     cc <- data.frame(res)
    #     cc2 <- subset(cc, rownames(cc) %in% select_taxa)
    #     
    #     interaction <- paste0(use[1], "_", use[2])
    #     P_adj <- cc2$padj
    #     results <- data.frame(Int = interaction, P_adj = P_adj)
    #     contrast_results <- rbind(contrast_results, results)
    #   }
      
      # group_num <- length(groupings)
      # 
      # resm <- matrix(NA, group_num, group_num)
      # resm[lower.tri(resm) ] <-round(contrast_results$P_adj, 4)
      # resm <- t(resm)
      # resm[lower.tri(resm) ] <-round(contrast_results$P_adj, 4)
      # rownames(resm) <- groupings
      # colnames(resm) <- groupings
      # resm
    #}
    #resm
  })
  
  
  # observe({
  #   event.data <- event_data("plotly_click", source = "t_exp")
  da_taxa_resm <- reactive({
    #validate(
    #  need(is.null(event.data) == F, "Click on the Heatmap")
    #)
    event.data <- event_data("plotly_click", source = "t_exp")
    resm <- da_taxa_stat_mat()
    
    resm_lab <- resm
    resm_lab[resm_lab < 0.0001] <- "<0.0001"
    
    
    if (any(resm<0.05, na.rm = T)){
      breaker = seq(0, 0.05, by = 0.0005)
      coler = c(colorRampPalette(c("red", "white"))(n=100))
    } else {
      breaker <- seq(0, 1, by = 0.01)
      coler <- c(colorRampPalette(c("white"))(n=100))
    }
    
    select_taxa <- event.data$y
    new_name <- unlist(strsplit(select_taxa, ";"))
    new_name <- new_name[length(new_name)]
    
    if (input$numInputs == 2){
      resm[1,2] = resm[1,2] + 0.0000001
    }
    
    par(cex.main=0.8)
    v = heatmap.3(data.matrix(resm),
              cellnote = resm_lab,
              notecol="black",
              #main=new_name,
              #key = TRUE,
              #keysize = 1.0,
              key.title = NULL,
              breaks = breaker,
              col = coler,
              #breaks = seq(0, 0.05, by = 0.0005),
              #col=c(colorRampPalette(c("red", "white"))(n=100)),
              dendrogram = 'none',
              Rowv=F,
              Colv=F,
              margins=c(10,10),
              cexRow=1.2,
              cexCol=1.2#,
              #na.color="gray60"
    )
    v
    title(new_name, line= -2.5)
  })
 # })
  
    
    # If NULL dont do anything
  output$da_taxa_stat_heat <- renderPlot({
    event.data <- event_data("plotly_click", source = "t_exp")
    validate(
      need(is.null(event.data) == F, "Click on the Heatmap")
    )
    plott <- da_taxa_resm()
    plott
  })

  
  output$da_taxa_stat_ui <- renderUI({
    plotOutput("da_taxa_stat_heat")
  })
  
  
  output$species_heat_download <- downloadHandler(
    filename = function() { paste("differential_taxa", '.png', sep='') },
    content = function(file) {
      species = length(unique(spec_taxa_data()$Taxa))
      
      p3 <- da_taxa_heat_plotly()

      export(p3, file = file)
    }
  )
  
  output$species_stat_download <- downloadHandler(
    filename = function() { paste("differential_taxa_heat", '.png', sep='') },
    content = function(file) {
      event.data <- event_data("plotly_click", source = "t_exp")
      resm <- da_taxa_stat_mat()
      
      resm_lab <- resm
      resm_lab[resm_lab < 0.0001] <- "<0.0001"
      
      
      if (any(resm<0.05, na.rm = T)){
        breaker = seq(0, 0.05, by = 0.0005)
        coler = c(colorRampPalette(c("red", "white"))(n=100))
      } else {
        breaker <- seq(0, 1, by = 0.01)
        coler <- c(colorRampPalette(c("white"))(n=100))
      }
      
      select_taxa <- event.data$y
      new_name <- unlist(strsplit(select_taxa, ";"))
      new_name <- new_name[length(new_name)]
      
      if (input$numInputs == 2){
        resm[1,2] = resm[1,2] + 0.0000001
      }
      
      png(file, width = 25, height = 20, units = 'cm', res = 300)
      par(cex.main=0.8)
      v = heatmap.3(data.matrix(resm),
                    cellnote = resm_lab,
                    notecol="black",
                    #main=new_name,
                    #key = TRUE,
                    #keysize = 1.0,
                    key.title = NULL,
                    breaks = breaker,
                    col = coler,
                    #breaks = seq(0, 0.05, by = 0.0005),
                    #col=c(colorRampPalette(c("red", "white"))(n=100)),
                    dendrogram = 'none',
                    Rowv=F,
                    Colv=F,
                    margins=c(10,10),
                    cexRow=1.2,
                    cexCol=1.2#,
                    #na.color="gray60"
      )
      v
      title(new_name, line= -2.5)
      dev.off()
    }
  )
  
  output$species_stat_results <- downloadHandler(
    filename = function() { paste("differential_taxa_anova", '.txt', sep='') },
    content = function(file) {
      stat_results <- da_taxa_stat()
      keepers <- subset(stat_results, p_adj < input$taxa_pval_up)
      write.table(keepers, file, row.names = FALSE, sep = '\t', quote = FALSE)
    }
  )
  
 
  # #### Taxa Search UI elements! ####
  # 
  # # Taxa Search Plot
  # 
  # #reorder full file with group allocations
  # reorder_full_mat <- reactive({
  #   req(grouped_samps())
  #   full_nums <- full_file()[,4:ncol(full_file())]
  #   exprs_reorder = full_nums[,c(unlist(grouped_samps()))]
  #   exp2 = data.frame(exprs_reorder)
  #   return(exp2)
  # })
  # 
  # # Taxa Search Dataframe
  # spec <- reactive({
  #   validate(
  #     need(input$acc_list, 'Please Select at least one Gene Family in the Gene Search Tab')
  #   )
  # 
  #   spec_list = input$acc_list
  #   acc_column = as.vector(acc_select())
  #   gfam_df = data.frame(full_file())
  #   
  #   #print(dim(full_file()))
  #   #print(dim(reorder_full_mat()))
  #   
  #   spec_plot_df = reorder_full_mat()
  #   spec_plot_df$Acc = gfam_df$Acc
  #   spec_plot_df$Gene_Family = paste(" ", gfam_df$Gene_Family, " ", sep='')
  #   spec_plot_df$Species = gfam_df$Species
  # 
  #   specs = spec_plot_df[which(spec_plot_df$Gene_Family %in% spec_list), ]
  #   
  #   spec_col = (ncol(specs))-3
  #   spec_row = nrow(specs)
  # 
  #   spec_RelExp = data.frame(specs[1:(spec_col)])
  #   spec_RelExp2 = t(spec_RelExp)
  #   df_spec_RelExp = data.frame(spec_RelExp2)
  # 
  #   spec_datas = c()
  #   for (i in 1:spec_row) {
  #     spec_data=as.vector(t(df_spec_RelExp[i]))
  #     spec_datas = c(spec_datas,spec_data)
  #   }
  #   
  #   groupings = new_group_names()
  #   grouping_nums = group_dims()
  # 
  #   group_titles = c()
  #   for (i in 1:length(groupings)) {
  #     title = groupings[i]
  #     reps = grouping_nums[i]
  #     title_rep = rep(title, reps)
  #     group_titles = c(group_titles, title_rep)
  #   }
  #   
  #   group_titles2 = rep(group_titles, spec_row)
  # 
  #   spec_Accs = as.vector(t(specs[1+spec_col]))
  #   spec_Acc = rep(spec_Accs, each = spec_col)
  #   spec_sample_num = paste(1:(spec_col))
  #   spec_samp=strtoi(spec_sample_num)
  #   spec_isamp = rep(spec_samp,spec_row)
  # 
  #   spec_specs = as.vector(t(specs[3+spec_col]))
  #   spec_spec = rep(spec_specs, each = spec_col)
  # 
  #   spec_gfams0 = as.vector(t(specs[2+spec_col]))
  #   spec_gfam = rep(spec_gfams0, each = spec_col)
  #   
  # 
  #   spec_g_data = data.frame(spec_datas, spec_Acc, spec_gfam, 
  #                            spec_isamp, group_titles2, spec_spec)
  #   colnames(spec_g_data) = c("Exp","Acc", "Gfam", "Sample_num", "Groups", "Species")
  #   
  #   spec_g_data
  # })
  # 
  # # Generate necessary color for Taxa Search Plot
  # spec_colors <- reactive({
  #   spec_df <- spec()
  #   most_specs <- max(table(spec_df$Gfam))
  #   cols <- randomColor(most_specs)
  #   cols
  # })
  # 
  # # Generate necessary UI elements for correct number of Taxa Search Plots
  # observeEvent({
  #   input$acc_list
  #   new_group_names()
  #   }, {
  #   if (input$testme) {}
  #   else{
  #     validate(
  #       need(input$file1, "Please provide a file in the Upload Tab")
  #       )
  #   }
  #   output$spec_plot <- renderUI({
  #     if (input$input_type == "Biobakery"){
  #       validate(
  #         need(length(input$acc_list) > 0, 'Please select at least one gene family on the previous tab')
  #       )
  #       get_spec_plot_list(max_plots, 
  #                          length(input$acc_list), 
  #                          spec(),
  #                          spec_colors())
  #     }
  #     else{
  #       textOutput("ebi_message")
  #     }
  #     })
  # })
  # 
  # output$ebi_message <- renderText({
  #   paste("Features by Taxa not Available from EBI Metagenomic Inputs", sep='\n')
  # })
  # 
  # output$search_taxa_download <- renderUI({
  #   if (input$input_type == "Biobakery"){
  #     lapply(1:length(input$acc_list), function(i) {
  #       display_name = input$acc_list
  #       downloadButton(paste0("downloadTaxa", i), paste("Download", display_name[i], sep=" "))
  #     })
  #   }
  #   else{
  #     textOutput("ebi_message")
  #   }
  # })
  # 
  # observe({
  #   if (input$input_type == "Biobakery"){
  #     lapply(1:length(input$acc_list), function(i) {
  #       output[[paste0("downloadTaxa", i)]] <- downloadHandler(
  #         filename = function() { paste(input$acc_list[i], "_taxa", '.png', sep='') },
  #         content = function(file) {
  #           spec <- spec()
  #           spec_curr <- subset(spec, spec$Gfam==input$acc_list[i])
  #           specs_nums <- length(unique(spec_curr$Species))
  #           png(file, width = 25+0.3*specs_nums, height = 30, units ='cm', res = 300)
  #           plotter<-download_spec_plot_list(max_plots, 
  #                                            1,
  #                                            #subset(spec, spec$Gfam==input$acc_list[i]),
  #                                            spec_curr,
  #                                            spec_colors())
  #           print(plotter)
  #           dev.off()
  #         }
  #       )
  #     })
  #   }
  #   else{
  #     textOutput("ebi_message")
  #   }
  # })
  # 
  # 
  #### Correlation Plot #####

  ## Correlation selections
  sig_tab <- reactive({
    gfam_DF = acc_full()
    samp_paths = gfam_DF
    #samp_paths = gfam_DF[apply(gfam_DF==0,1,sum)<=(0.90*length(group_names())),]
    #samp_paths$Gene_Family <- paste(" ", samp_paths$Gene_Family, " ", sep="")
    samp_paths$Feature
  })

  ## Update selections ##
  # observe({
  #   if (input$testme) {
  #     updateSelectizeInput(session,'sig_select', choices = sig_tab(), server = TRUE)
  #   }
  #   else {
  #     if (is.null(input$file1)) {}
  #     else{
  #       updateSelectizeInput(session,'sig_select', choices = sig_tab(), server = TRUE)
  #     }
  #   }
  #   
  # })

  # All sample Correlation Plot #

  actual_corr_plot <- reactive({
    if (input$testme) {
      validate(
        need(length(group_dims())==4, "Please visit the group tab to verify group assignment")
      )
    }
    else{
      if (input$input_type == "Biobakery"){
      validate(
        need(input$file1, "Please provide a file in the Upload Tab") %then%
        need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
      if (input$input_type == "EBI"){
        validate(
          need(input$features1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
    }

    validate(
      need(length(input$sig_select) > 1, 'Select at least two Features!')
    )
    corr_list = input$sig_select

    gfam_DF1 = acc_full()
    col_num = ncol(gfam_DF1)
    row_num = nrow(gfam_DF1)
    reorder_samps = gfam_DF1[,-1]
    reorder_samps = reorder_samps[,c(unlist(grouped_samps()))]
    header_DF = gfam_DF1[,1]

    gfam_DF = cbind(header_DF, reorder_samps)
    samp_paths = gfam_DF
    #samp_paths$Gene_Family = paste(" ", samp_paths$Gene_Family, " ", sep="")
    
    #heyo = samp_paths[grep(paste(corr_list, collapse = '|'), samp_paths$Feature), ]
    heyo <- subset(samp_paths, header_DF %in% corr_list)
    colnames(heyo)[1] = "Feature"
    
    col = (ncol(samp_paths))-1
    row = nrow(heyo)

    heyo_small <- heyo
    rownames(heyo_small) = heyo$Feature
    heyo_small = heyo_small[,-1]

    heyo_side = data.frame(t(heyo_small))
    colnames(heyo_side) = rownames(heyo_small)

    library(psych)
    hah2 = corr.test(heyo_side, method = "spearman")

    corr_mat = as.matrix(hah2$r)
    rownames(corr_mat) = paste(1:nrow(corr_mat))
    colnames(corr_mat) = paste(1:ncol(corr_mat))

    orig_p = hah2$p
    orig_p[col(orig_p) == row(orig_p) | upper.tri(orig_p)] <- NA
    new_ps = p.adjust(orig_p, method = "BH", n = length(orig_p))
    new_p_mat = matrix(new_ps, ncol = ncol(orig_p), nrow = nrow(orig_p))

    new_p_mat[upper.tri(new_p_mat)] = t(new_p_mat)[upper.tri(new_p_mat)]
    diag(new_p_mat) = rep(1,nrow(new_p_mat))

    p_list = c(new_p_mat)
    p_dims = dim(new_p_mat)

    sym_list = c()
    for (i in p_list){
      if (is.na(i)){
        sym = ""
      }
      else {
        if (i > 0){
          sym = ""
          if (i < 0.05){
            sym = "*"
            if (i<0.01){
              sym = "**"
              if(i<0.001){
                sym = "***"
              }
            }
          }
        }
      }
      sym_list = c(sym_list, sym)
    }

    sym_mat = matrix(sym_list, ncol = p_dims[1], nrow = p_dims[2])

    library(gplots)
    my_palette <-colorRampPalette(c("blue", "white", "red"))(n=100)
    corr_mat[is.na(corr_mat)] <- 0
    v <- heatmap.2(corr_mat,
                   cellnote = sym_mat, notecol="black",
                   main="Correlation Values across \nAll Samples", #heatmap title
                   density.info="none",
                   #key = TRUE,
                   #keysize = 1.0,
                   breaks = seq(-1, 1, by = 0.02),
                   col=my_palette,
                   dendrogram="both",
                   margins=c(5,5),
                   cexRow=1,
                   cexCol=1.2,
                   trace=c("none"),
                   na.color="gray60"
    )
    v
  })

  output$corr_plot <- renderPlot({
    actual_corr_plot()
  })

  actual_corr_names <- reactive({
    if (input$testme) {
      validate(
        need(length(group_dims())==4, "Please visit the group tab to verify group assignment")
      )
    }
    else{
      if (input$input_type == "Biobakery"){
        validate(
          need(input$file1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
      if (input$input_type == "EBI"){
        validate(
          need(input$features1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
    }

    validate(
      need(length(input$sig_select) > 1, 'Select at least two Gene Families!')
    )
    corr_list = input$sig_select

    gfam_DF1 = acc_full()
    col_num = ncol(gfam_DF1)
    row_num = nrow(gfam_DF1)
    reorder_samps = gfam_DF1[,-1]
    reorder_samps = reorder_samps[,c(unlist(grouped_samps()))]
    header_DF = gfam_DF1[,1]

    gfam_DF = cbind(header_DF, reorder_samps)
    samp_paths = gfam_DF
    #samp_paths$Gene_Family = paste(" ", samp_paths$Gene_Family, " ", sep="")


    #heyo = samp_paths[grep(paste(corr_list, collapse = '|'), samp_paths$Gene_Family), ]
    heyo <- subset(samp_paths, header_DF %in% corr_list)
    colnames(heyo)[1] = "Feature"

    actual_corr_names <- heyo$Feature
    actual_corr_names
  })


  # Generate list of correlation matrices for each Group
  group_corr_plist <- reactive({
    if (input$testme) {}
    else{
      if (input$input_type == "Biobakery"){
        validate(
          need(input$file1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
      if (input$input_type == "EBI"){
        validate(
          need(input$features1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
    }

    validate(
      need(length(input$sig_select) > 1, 'Select at least two Gene Families!')
    )
    corr_list2 = input$sig_select

    groupings = new_group_names()
    grouping_nums = group_dims()

    group_titles = c()
    for (i in 1:length(groupings)) {
      title = groupings[i]
      reps = grouping_nums[i]
      title_rep = rep(title, reps)
      group_titles = c(group_titles, title_rep)
    }

    group_titles_uniq = unique(group_titles)

    corr_mat_list=list()
    for (j in 1:input$numInputs){

      gfam_DF1 = acc_full()
      col_num = ncol(gfam_DF1)
      row_num = nrow(gfam_DF1)
      reorder_samps = gfam_DF1[,2:col_num]
      reorder_samps = reorder_samps[,c(unlist(grouped_samps()))]
      header_DF = gfam_DF1[,1]

      gfam_DF = cbind(header_DF, reorder_samps)
      samp_paths = gfam_DF
      #samp_paths$Gene_Family = paste(" ", samp_paths$Gene_Family, " ", sep="")

      #heyo = samp_paths[grep(paste(corr_list2, collapse = '|'), samp_paths$Gene_Family), ]
      heyo <- subset(samp_paths, header_DF %in% corr_list2)
      colnames(heyo)[1] = "Feature"

      col = (ncol(samp_paths))-1
      row = nrow(heyo)

      heyo_small<-heyo
      rownames(heyo_small) = heyo$Feature
      heyo_small = heyo_small[,-1]

      heyo_side = data.frame(t(heyo_small))
      colnames(heyo_side) = rownames(heyo_small)

      heyo_side$Group = group_titles
      heyo_side_real = subset(heyo_side, Group == group_titles_uniq[j], select = -c(Group))


      library(psych)
      hah2 = corr.test(heyo_side_real, method = "spearman")

      corr_mat = as.matrix(hah2$r)
      rownames(corr_mat) = paste(1:nrow(corr_mat))
      colnames(corr_mat) = paste(1:ncol(corr_mat))

      corr_mat_list[[j]] <- corr_mat
    }
    corr_mat_list
  })

  # Generate list of correlation significance symbols for each group
  group_sym_plist <- reactive({
    if (input$testme) {}
    else{
      if (input$input_type == "Biobakery"){
        validate(
          need(input$file1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
      if (input$input_type == "EBI"){
        validate(
          need(input$features1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
    }

    validate(
      need(length(input$sig_select) > 1, 'Select at least two Gene Families!')
    )
    corr_list2 = input$sig_select

    groupings = new_group_names()
    grouping_nums = group_dims()

    group_titles = c()
    for (i in 1:length(groupings)) {
      title = groupings[i]
      reps = grouping_nums[i]
      title_rep = rep(title, reps)
      group_titles = c(group_titles, title_rep)
    }

    group_titles_uniq = unique(group_titles)

    sym_mat_list=list()
    for (j in 1:input$numInputs){

      gfam_DF1 = acc_full()
      col_num = ncol(gfam_DF1)
      row_num = nrow(gfam_DF1)
      reorder_samps = gfam_DF1[,2:col_num]
      reorder_samps = reorder_samps[,c(unlist(grouped_samps()))]
      header_DF = gfam_DF1[,1]

      gfam_DF = cbind(header_DF, reorder_samps)
      samp_paths = gfam_DF
      #samp_paths = gfam_DF[apply(gfam_DF==0,1,sum)<=(0.90*length(group_names())),]
      #samp_paths$Gene_Family = paste(" ", samp_paths$Gene_Family, " ", sep="")

      #heyo = samp_paths[grep(paste(corr_list2, collapse = '|'), samp_paths$Gene_Family), ]
      
      heyo <- subset(samp_paths, header_DF %in% corr_list2)
      colnames(heyo)[1] = "Feature"

      col = (ncol(samp_paths))-1
      row = nrow(heyo)

      heyo_small <- heyo
      rownames(heyo_small) = heyo$Feature
      heyo_small = heyo_small[,-1]

      heyo_side = data.frame(t(heyo_small))
      colnames(heyo_side) = rownames(heyo_small)

      heyo_side$Group = group_titles
      heyo_side_real = subset(heyo_side, Group == group_titles_uniq[j], select = -c(Group))


      library(psych)
      hah2 = corr.test(heyo_side_real, method = "spearman")


      corr_mat = as.matrix(hah2$r)
      rownames(corr_mat) = paste(1:nrow(corr_mat))
      colnames(corr_mat) = paste(1:ncol(corr_mat))

      orig_p = hah2$p
      orig_p[col(orig_p) == row(orig_p) | upper.tri(orig_p)] <- NA
      new_ps = p.adjust(orig_p, method = "BH", n = length(orig_p))
      new_p_mat = matrix(new_ps, ncol = ncol(orig_p), nrow = nrow(orig_p))

      new_p_mat[upper.tri(new_p_mat)] = t(new_p_mat)[upper.tri(new_p_mat)]
      diag(new_p_mat) = rep(1,nrow(new_p_mat))

      p_list = c(new_p_mat)
      p_dims = dim(new_p_mat)

      sym_list = c()
      for (i in p_list){
        if (is.na(i)){
          sym = ""
        }
        else {
          if (i > 0){
            sym = ""
            if (i < 0.05){
              sym = "*"
              if (i<0.01){
                sym = "**"
                if(i<0.001){
                  sym = "***"
                }
              }
            }
          }
        }
        sym_list = c(sym_list, sym)
      }

      sym_mat = matrix(sym_list, ncol = p_dims[1], nrow = p_dims[2])
      sym_mat_list[[j]] <- sym_mat
    }
    sym_mat_list
  })


  output$sig_message <- renderText({
    paste("* = p < 0.05", "** = p < 0.01", "*** = p < 0.001", sep='\n')
  })

  # Generate Table of Selected Gene Families
  corr_label_table <- reactive({
    if (input$testme) {}
    else{
      if (input$input_type == "Biobakery"){
        validate(
          need(input$file1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
      if (input$input_type == "EBI"){
        validate(
          need(input$features1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
    }

    validate(
      need(length(input$sig_select) > 1, 'Select at least two Gene Families!')
    )

    Gene_Families = actual_corr_names()
    nums = length(Gene_Families)
    Label = paste(1:nums)
    Label_Matrix = data.frame(Gene_Families, Label)
    #print(Label_Matrix)
    Label_Matrix
  })

  output$corr_labels <- renderTable({
    corr_label_table()
  })

  output$corr_table_download <- downloadHandler(
    filename = function() { paste("correlation_sample_key", '.txt', sep='') },
    content = function(file) {
      corr_tab = corr_label_table()
      corr_tab2 = data.frame(corr_tab)
      write.table(corr_tab2, file, row.names = FALSE, sep = '\t', quote = FALSE)
    }
  )

  observeEvent({
    #new_group_names()
    input$sig_select
    }, {
    if (input$testme) {
      validate(
        need(length(group_dims())==4, "")
      )
    }
      else{
        if (input$input_type == "Biobakery"){
          validate(
            need(input$file1, "Please provide a file in the Upload Tab") %then%
              need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
        }
        if (input$input_type == "EBI"){
          validate(
            need(input$features1, "Please provide a file in the Upload Tab") %then%
              need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
        }
      }

    if (input$numInputs > 1) {
      #checker = input$sig_select
      #display_name = group_names()
      output$group_corrs <- renderUI({ get_plot_output_list(max_plots,
                                                            input$numInputs,
                                                            group_corr_plist(),
                                                            group_sym_plist(),
                                                            new_group_names())
      })
    }
  })


  # Download full correlation heatmap
  output$corr_download <- downloadHandler(
    filename = function() { paste("all_gene_family_correlation", '.png', sep='') },
    content = function(file) {
      png(file, width = 25, height = 20, units ='cm', res = 300)
      corr_list = input$sig_select

      gfam_DF1 = acc_full()
      col_num = ncol(gfam_DF1)
      row_num = nrow(gfam_DF1)
      reorder_samps = gfam_DF1[,2:col_num]
      reorder_samps = reorder_samps[,c(unlist(grouped_samps()))]
      header_DF = gfam_DF1[,1]

      gfam_DF = cbind(header_DF, reorder_samps)
      samp_paths = gfam_DF
      #samp_paths = gfam_DF[apply(gfam_DF==0,1,sum)<=(0.90*length(group_names())),]
      #samp_paths$Gene_Family = paste(" ", samp_paths$Gene_Family, " ", sep="")


      #heyo = samp_paths[grep(paste(corr_list, collapse = '|'), samp_paths$Gene_Family), ]
      
      heyo <- subset(samp_paths, header_DF %in% corr_list)
      colnames(heyo)[1] = "Feature"

      col = (ncol(samp_paths))-1
      row = nrow(heyo)

      heyo_small<-heyo
      rownames(heyo_small) = heyo$Feature
      heyo_small = heyo_small[,-1]

      heyo_side = data.frame(t(heyo_small))
      colnames(heyo_side) = rownames(heyo_small)

      library(psych)
      hah2 = corr.test(heyo_side, method = "spearman")


      corr_mat = as.matrix(hah2$r)
      rownames(corr_mat) = paste(1:nrow(corr_mat))
      colnames(corr_mat) = paste(1:ncol(corr_mat))

      orig_p = hah2$p
      orig_p[col(orig_p) == row(orig_p) | upper.tri(orig_p)] <- NA
      new_ps = p.adjust(orig_p, method = "BH", n = length(orig_p))
      new_p_mat = matrix(new_ps, ncol = ncol(orig_p), nrow = nrow(orig_p))

      new_p_mat[upper.tri(new_p_mat)] = t(new_p_mat)[upper.tri(new_p_mat)]
      diag(new_p_mat) = rep(1,nrow(new_p_mat))
      #new_p_dat = data.frame(new_p_mat)

      p_list = c(new_p_mat)
      p_dims = dim(new_p_mat)

      sym_list = c()
      for (i in p_list){
        if (is.na(i)){
          sym = ""
        }
        else {
          if (i > 0){
            sym = ""
            if (i < 0.05){
              sym = "*"
              if (i<0.01){
                sym = "**"
                if(i<0.001){
                  sym = "***"
                }
              }
            }
          }
        }
        sym_list = c(sym_list, sym)
      }

      sym_mat = matrix(sym_list, ncol = p_dims[1], nrow = p_dims[2])

      library(gplots)
      my_palette <-colorRampPalette(c("blue", "white", "red"))(n=100)
      corr_mat[is.na(corr_mat)] <- 0
      v <- heatmap.2(corr_mat,
                     cellnote = sym_mat, notecol="black",
                     main="Correlation Values across \nAll Samples", #heatmap title
                     density.info="none",
                     #key = TRUE,
                     #keysize = 1.0,
                     breaks = seq(-1, 1, by = 0.02),
                     col=my_palette,
                     dendrogram="both",
                     margins=c(5,5),
                     cexRow=1,
                     cexCol=1.2,
                     trace=c("none"),
                     na.color="gray60"
      )
      v
      dev.off()
    }
  )

  ## UI Elements for downloading group correlation ##

  output$group_download <- renderUI({
    lapply(1:input$numInputs, function(i) {
      display_name = new_group_names()
      downloadButton(paste0("downloadData", i), paste("Download", display_name[i], sep=" "))
    })
  })

  observe({
    lapply(1:input$numInputs, function(i) {
      output[[paste0("downloadData", i)]] <- downloadHandler(
        filename = function() { paste(new_group_names()[i], "_correlation", '.png', sep='') },
        content = function(file) {
          png(file, width = 25, height = 20, units ='cm', res = 300)
          download_plot_output_list(max_plots,
                                    1,
                                    group_corr_plist()[i],
                                    group_sym_plist()[i],
                                    new_group_names()[i])
          dev.off()
        }
      )
    })
  })
  
}