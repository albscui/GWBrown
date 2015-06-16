## Load packages
library(flowCore)
library(flowStats)
library(flowViz) # for flow data visualization
library(utils)

calcRatio <- function(fcs_filename, green_filter_name = "B510_20.A", 
                      red_filter_name = "YG582_15.A"){

  #Read FCS data and export to table
  frame <- read.FCS(fcs_filename, alter.names = T)
  frame <- exprs(frame)
  
  #Extract GFP & RFP data columns, calculate ratio
  frame_gfp <- frame[ , which(colnames(frame) == green_filter_name)]
  frame_rfp <- frame[ , which(colnames(frame) == red_filter_name)]
  frame_ratio <- frame_gfp/frame_rfp
  
  #Make dataframe
  out_table <- data.frame(CellID="",
                          GFP=frame_gfp, 
                          RFP=frame_rfp,
                          GFPRFPRatio=frame_ratio)
  #Remove all zeros
  out_table <- out_table[which(out_table$GFP > 0 &
                                   out_table$RFP > 0 &
                                   out_table$GFPRFPRatio > 0), ]
  out_table$CellID <- c(1:nrow(out_table))
  

  #Write CSV
  out_file_name <- substring(fcs_filename, 12, 14)
  write.csv(out_table, paste(out_file_name, ".csv", sep=""), row.names=F)
  print(paste("Progress ", 
              100*round(which(fcs_files == fcs_filename)/length(fcs_files), digits=3), 
              "%", collapse=""))
}

#Choose directory with FCS files and save filenames to a vector
fcs_dir <- choose.dir(getwd(), "Choose folder with FCS files")
fcs_files <- list.files(fcs_dir, pattern = "fcs$")

#Create a csv with GFP:RFP ratios across all fcs files
setwd(fcs_dir)
sapply(fcs.files, calcRatio)

