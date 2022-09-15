###-------------AneufinderCopyNumberCounter------###

#----Thomas van Ravesteyn
#----Kops Group
#----Hubrecht Institute


AneufinderCopyNumberCounter <- function(sampleIDs){
  
  ###-----Required packages----------------------###
  require("tidyr")
  require("ggplot2")
  require("AneuFinder")
  require("openxlsx")

  ###--------------------------------------------###
  message("=== Run Aneufinder Copy Number Counter ===")

  #check if output folder exists, if not it wil create a new output folder
  if (dir.exists(output_dir) == FALSE){
    message("=== Create output directory ")
    message("- ", output_dir )
    dir.create(output_dir)
  }
  
  if (analysis.type == "absolute"){
    file_end <- "abs"
  }
  if (analysis.type == "relative"){
    file_end <- "rel"
  }
  
  #verify presence of model files for each sample, if not remove sample
  message("=== Verify presence of model files for each sample ===")
  for (sample in sampleIDs){
    if (selected.model == "dnaCopy"){
      if (length(dnaFiles[[sample]]) == 0){
        sampleIDs <- sampleIDs[! sampleIDs %in% sample]
        message(paste0("- ", sample, " contains no dnaCopy model files"))
      } 
    } else if (selected.model == "hmm") {
      if (length(hmmFiles[[sample]]) == 0){
        sampleIDs <- sampleIDs[! sampleIDs %in% sample]
        message(paste0("- ",sample, " contains no HMM model files"))
      } 
    } else if (selected.model == "edivisive"){
      if (length(edivisiveFiles[[sample]]) == 0){
        sampleIDs <- sampleIDs[! sampleIDs %in% sample]
        message(paste0("- ",sample, " contains no edivisive model files"))
      } 
    } else {
      stop("Incorrect model selection")
    }
    
  }
  
  #create lists to combine all data
  my_alterations <- list()
  my_alterations_sum <- list()
  my_alterations_sum_table <- data.frame()
  my_sample_length <- data.frame()
  
  #Loop over all samples and analyze per sample
  message("=== Count CNAs for each sample and write to .xlxs file ===")
  for(sample in sampleIDs) {
    
    #select files based on model selection
    message("- Load files for ", sample)
    if (selected.model == "dnaCopy"){
      my_files <- unlist(dnaFiles[[sample]])
    } else if (selected.model == "hmm") {
      my_files <- unlist(hmmFiles[[sample]])
    } else if (selected.model == "edivisive"){
      my_files <- unlist(edivisiveFiles[[sample]])
    } else {
      stop("Incorrect model selection")
    }
    
    #Determine number of model files and append to sample length dataframe
    sample_name <- sample
    sample_length <- length(my_files)
    temp_sample_length <- data.frame(sample_name ,sample_length)
    my_sample_length <- rbind(my_sample_length, temp_sample_length )
    
    #Determine copy number alterations for each model file
    for (my_model in my_files){
      
      #load data from aneufinder model
      load(my_model)
      
      #extract all segments from model file
      my_segments <- as.data.frame(model[["segments"]])
      
      #determine number of chromsomes
      my_chromosomes <- as.character(unique(my_segments$seqnames))
      
      ##loop over all chromosomes to identify which copy numbers are present
      chr_list <- list()
      
      for (my_chr in my_chromosomes){
        
        #set n-state for current chromosome
        if (exists(paste0("chr",my_chr,"_n"))){
          n_base <- get0(paste0("chr",my_chr,"_n"))
        } else {
          n_base <- chr_n
        }
        
        #identify the number of segments and corresponding copy number
        initial_segments <- my_segments[my_segments$seqnames == my_chr, ]
        
        #select segments that meet minimal segment size criteria
        temp_segments <- initial_segments[initial_segments$width > segment_min_size, ]
        
        #if there are multiple segments after size selection, 
        #check if the 2nd and later elements have identical copy numbers to the previous one,
        #if Yes, collect in segment number in remove_segments vector
        remove_segments <- c()
        if (nrow(temp_segments)>1){
          for (i in 2:nrow(temp_segments))
            if ((temp_segments[i,]$copy.number) == (temp_segments[c(i-1),]$copy.number)){
              remove_segments <- c(remove_segments, c(i-1))
            }
        }
        #if there are segments to be removed, remove segments from temp_segments to get final segment selection
        if (length(remove_segments > 0)){
          final_segments <- temp_segments[-c(remove_segments),]
        } else (final_segments <- temp_segments)
        
        #determine number of segments in final selection
        number_of_segments <- nrow(final_segments)
        
        #identify copynumber for each selected segment
        segment_cn <- final_segments[final_segments$seqnames == my_chr, "copy.number"]
        
        #add statistics to chromosome dataframe, one chr per row
        chr_list[[my_chr]]$number_of_segments <-  number_of_segments
        chr_list[[my_chr]]$segment_cn <- segment_cn
        
        #determine for each chr if there are (partial) gains and (partial) losses
        #first - determine if there are one or more segments per chromosome
        #second - assign chromosome gains and losses based on copy number for each segment
        #start with zero gains and losses
        chr_list[[my_chr]]$chr_whole_loss <- c(0)
        chr_list[[my_chr]]$chr_whole_gain <- c(0)
        chr_list[[my_chr]]$chr_partial_loss <- c(0)
        chr_list[[my_chr]]$chr_partial_gain <- c(0)
        
        #RUN analysis type "absolute"
        if (analysis.type == "absolute"){
        
          if (number_of_segments == 1){
            if (segment_cn < n_base){
              chr_list[[my_chr]]$chr_whole_loss <- chr_list[[my_chr]]$chr_whole_loss + abs(segment_cn - n_base)
            } 
            if (segment_cn > n_base){
              chr_list[[my_chr]]$chr_whole_gain <- chr_list[[my_chr]]$chr_whole_gain + abs(segment_cn - n_base)
            }
          } else {
            #if there are multiple segments
            
            #check if ALL segments have copy number below baseline,
            #if TRUE use highest copy number to calculate distance from baseline for whole chr loss
            if (!(FALSE %in% (segment_cn < n_base))){
              chr_list[[my_chr]]$chr_whole_loss <-  chr_list[[my_chr]]$chr_whole_loss + abs(max(segment_cn) - n_base)
            
              #remove segment from segment_cn that was already counted as whole chr loss
              segment_cn <- segment_cn[segment_cn != max(segment_cn)]
              
            }
            #check if ALL segments have copy number above baseline
            #if TRUE use lowest copy number to calculate distance from baseline for whole chr gain
            if (!(FALSE %in% (segment_cn > n_base))){
              chr_list[[my_chr]]$chr_whole_gain <- chr_list[[my_chr]]$chr_whole_gain + abs(min(segment_cn) - n_base)
              
              #remove segment from segment_cn that was already counted as whole chr gain
              segment_cn <- segment_cn[segment_cn != min(segment_cn)]
              
            }
          
            #loop over all segments to count partial losses and gains
            for (segment in segment_cn){
              if (segment < n_base){
                chr_list[[my_chr]]$chr_partial_loss <- chr_list[[my_chr]]$chr_partial_loss + abs(segment - n_base)
              } 
              if (segment > n_base){
                chr_list[[my_chr]]$chr_partial_gain <- chr_list[[my_chr]]$chr_partial_gain + abs(segment - n_base)
              }
                
             
              
            }
            
          }
        
        }
        #RUN analysis type "relative"
        if (analysis.type == "relative"){
          #if there is only 1 segment
          if (number_of_segments == 1){
            if (segment_cn < n_base){
              chr_list[[my_chr]]$chr_whole_loss <- 1
            } 
            if (segment_cn > n_base){
              chr_list[[my_chr]]$chr_whole_gain <- 1
            }
          } else {
            #if there are multiple segments
            
            #check if ALL segments have copy number below baseline
            if (!(FALSE %in% (segment_cn < n_base))){
              chr_list[[my_chr]]$chr_whole_loss <- 1
            }
            #check if ALL segments have copy number above baseline
            if (!(FALSE %in% (segment_cn > n_base))){
              chr_list[[my_chr]]$chr_whole_gain <- 1
            }
            
            #loop over all segments to count partial losses and gains
            for (segment in segment_cn){
              if (segment < n_base){
                chr_list[[my_chr]]$chr_partial_loss <- 1
              } 
              if (segment > n_base){
                chr_list[[my_chr]]$chr_partial_gain <- 1
              }
              
            }
            
          }
          
        }
        
      }
      
      #store identified alterations for current model
      my_alterations[[sample]][[my_model]]  <- chr_list
    }
    
    #save my alterations for current sample as R data object
    my_set <- my_alterations[[sample]]
    save(my_set, list = "my_alterations", file = paste0(output_dir,"/",Sys.Date(),"_" , sample ,"_alterations per cell_", file_end, ".Rdata"))
    rm(my_set)
    
    #Calculate total number of CNA within current sample
    #Loop over all chromosomes, collect data for each model 
    CNA_overview <- data.frame()
    sample_counts <- c()
    for (my_chr in names(my_alterations[[sample]][[1]])){
     
      sum_chr_whole_loss <- c(0)
      sum_chr_whole_gain <- c(0)
      sum_chr_partial_loss <- c(0)
      sum_chr_partial_gain <- c(0)
      
      #Loop over all cells, sum all CNA alterations for current chromosome
      for (my_cell in names(my_alterations[[sample]])){
        
        sum_chr_whole_loss <- sum_chr_whole_loss + my_alterations[[sample]][[my_cell]][[my_chr]]$chr_whole_loss
        sum_chr_whole_gain <- sum_chr_whole_gain + my_alterations[[sample]][[my_cell]][[my_chr]]$chr_whole_gain
        sum_chr_partial_loss <- sum_chr_partial_loss + my_alterations[[sample]][[my_cell]][[my_chr]]$chr_partial_loss
        sum_chr_partial_gain <- sum_chr_partial_gain + my_alterations[[sample]][[my_cell]][[my_chr]]$chr_partial_gain
      }
      
      #create temporary dataframe with total for CNAs and combine with CNA overview dataframe -- CNA summary list
      temp_CNA_overview <- data.frame()
      temp_CNA_overview <- data.frame("chr" = my_chr,
                                      "sum_chr_whole_loss" = sum_chr_whole_loss,
                                      "sum_chr_whole_gain" = sum_chr_whole_gain,
                                      "sum_chr_partial_loss" = sum_chr_partial_loss,
                                      "sum_chr_partial_gain" = sum_chr_partial_gain) 
      rownames(temp_CNA_overview) <- c(my_chr)
      CNA_overview <- rbind(CNA_overview, temp_CNA_overview)
      
      #Collect CNA counts for current chr in a single vector and append to previous chromosomes -- CNA summary table
      chr_counts <- c(sum_chr_whole_loss, sum_chr_whole_gain, sum_chr_partial_loss, sum_chr_partial_gain)
      sample_counts <- append(sample_counts, chr_counts)    
      }
    
    #store CNA overview -- CNA summary list
    my_alterations_sum[[sample]]$CNA_overview <- CNA_overview
    
    #save CNA overview as sample xlsx file -- CNA summary list
    write.xlsx(CNA_overview, file = paste0(output_dir,"/",Sys.Date(),"_" , sample ,"_CNA-overview_", file_end, ".xlsx") , sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE, showNA = FALSE)
   
    #Append CNA counts for current sample -- CNA summary table
    my_alterations_sum_table <- rbind(my_alterations_sum_table, sample_counts)
    
  }
  
  #create vector with column names and subsequently assign to  CNA summary table
  my_CNA_names <- c()
  for (my_chr in names(my_alterations[[sample]][[1]])){
    my_CNA_names <- append(my_CNA_names, c(paste0(my_chr, "_whole_loss"),
                                           paste0(my_chr, "_whole_gain"),
                                           paste0(my_chr, "_partial_loss"),
                                           paste0(my_chr, "_partial_gain")))
  }  
  colnames(my_alterations_sum_table) <- my_CNA_names
  
  #Assign sampleIDs to CNA summary table
  row.names(my_alterations_sum_table) <- sampleIDs
  
  #Add number of cells per sample to CNA summary table
  my_alterations_sum_table$sample_length <- my_sample_length$sample_length
  
  #Calculate relative CNA gains and losses
  if (analysis.type == "relative"){
    #change to matrix to enable math
    my_alterations_sum_table <- as.matrix(my_alterations_sum_table)
    #normalize values by sample length
    my_alterations_sum_table <- my_alterations_sum_table/my_alterations_sum_table[,ncol(my_alterations_sum_table)]
    #remove last column (sample length) and change back to data.frame
    my_alterations_sum_table <- as.data.frame(subset(my_alterations_sum_table, select = -sample_length))
  }
  
  
  #Write CNA summary table to single excel
  message("=== Save CNA Summary table to .xlsx ===")
  write.xlsx(my_alterations_sum_table, file = paste0(output_dir,"/",Sys.Date(),"_" , "All_CNA-overview_", file_end,".xlsx") , sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE, showNA = FALSE)

  #Write CNA summary table containing only whole gains
  all_w_gains <- my_alterations_sum_table[,colnames(my_alterations_sum_table)[c(grep("whole_gain",colnames(my_alterations_sum_table)))]]
  write.xlsx(all_w_gains, file = paste0(output_dir,"/",Sys.Date(),"_" , "All_whole_gains_", file_end,".xlsx") , sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE, showNA = FALSE)
  
  #Write CNA summary table containing only whole losses
  all_w_losses <- my_alterations_sum_table[,colnames(my_alterations_sum_table)[c(grep("whole_loss",colnames(my_alterations_sum_table)))]]
  write.xlsx(all_w_losses, file = paste0(output_dir,"/",Sys.Date(),"_" , "All_whole_losses_", file_end,".xlsx") , sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE, showNA = FALSE)
  
  #Write CNA summary table containing only partial gains
  all_p_gains <- my_alterations_sum_table[,colnames(my_alterations_sum_table)[c(grep("partial_gain",colnames(my_alterations_sum_table)))]]
  write.xlsx(all_p_gains, file = paste0(output_dir,"/",Sys.Date(),"_" , "All_partial_gains_", file_end,".xlsx") , sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE, showNA = FALSE)
  
  #Write CNA summary table containing only partial losses
  all_p_losses <- my_alterations_sum_table[,colnames(my_alterations_sum_table)[c(grep("partial_loss",colnames(my_alterations_sum_table)))]]
  write.xlsx(all_p_losses, file = paste0(output_dir,"/",Sys.Date(),"_" , "All_partial_losses_", file_end,".xlsx") , sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE, showNA = FALSE)
  
  
  
  #create and write plots to single pdf
  message("=== Save CNA Counter Plots to PDF ===")
  grDevices::pdf(file=paste0(output_dir,"/",Sys.Date(),"_" ,"CNA_plots_", file_end,".pdf"), paper="a4")
  for (sample in sampleIDs) { 
    
    #create long format data table for ggplot2
    temp.df <- pivot_longer(my_alterations_sum[[sample]]$CNA_overview, names_to = "type", cols = c("sum_chr_whole_loss",  "sum_chr_partial_loss", "sum_chr_partial_gain","sum_chr_whole_gain" ))
    #sort factor for chromosomes from low to high
    temp.df$chr <- factor(temp.df$chr, levels = c(unique(temp.df$chr)))
    
    #create plot with CNA overview per chromosome
    p.CNA_overview <- ggplot(temp.df, aes(x= chr, y=value,  fill = factor(type,unique(type)))) +
      geom_bar(position = "stack", stat = "identity") +
      labs(title=paste0(sample," total number of CNA"), x="Chromosome", y = "Number of CNA", size = 10) +
      scale_fill_manual(values=c("#841859","#FFBFDE","#ABDFAC","#005600")) +
      labs(fill='CNA type')
      theme_minimal()
    
    print(p.CNA_overview)
  }
  grDevices::dev.off()
  
  #write used settings to .txt file
  writeLines(c(("=== CNA Counter Settings ==="),
               c(paste0("- Input directory = ",input_dir)),
               c(paste0("- Output directory = ",output_dir)),
               
               c(paste0("- Autosome copy number baseline = ",chr_n)),
               c(paste0("- Chr X copy number baseline = ",chrX_n)),
               c(paste0("- Minimum segment size = ",segment_min_size))),
    c(paste0(output_dir,"/",Sys.Date(),"_CNA-Counting-parameters_", file_end,".txt")))
  
  message("=== AneufinderCopyNumber script has finished ===")
}

