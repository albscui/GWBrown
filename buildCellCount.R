# Make a dataframe with cell counts for heatmap
buildCellCount <- function(filelist) {
    # Progress bar
    pb <- txtProgressBar(min = 0, max = 1, style = 3)
    count <- 0
    
    # empty data frame for storing values
    df_result <- data.frame()
    
    # loop through each file name in a list of files
    for (file in filelist) {
        
        # read into the csv file, initialize data frame
        classes <- c("NULL", "numeric", "numeric", "numeric")
        df <- read.csv(file, colClasses=classes)
        # loop through a string vector for every filter
        for (filter in c("GFP", "RFP", "GFPRFPRatio")) {
            
            colony <- getName(file)
            row <- substring(colony, 1, 1)
            col <- substring(colony, 2,  )
            gID <- getGeneID(colony)
            gGrp <- getGeneGroup(gID)
            cell_count <- nrow(df)
            summary <- summary(df[, filter])
            
            if (filter == "GFP") {
                # count how many cells have greater than 50,000 intensity
                cell_count_upper <- nrow(df[which(df$GFP > 50000), ])
            } else if (filter == "RFP") {
                # counts how many cells have greater than 5,000 intensity
                cell_count_upper <- nrow(df[which(df$RFP > 5000), ])
            } else if (filter == "GFPRFPRatio") {
                cell_count_upper <- nrow(df[which(df$GFPRFPRatio > 10), ])
            }
            
            # Put each parameter into a new row for final data frame
            df_result <- rbind(df_result,
                               data.frame(Filter=filter,
                                          Colony=colony,
                                          Row=row,
                                          Column=col,
                                          GeneID=gID,
                                          GeneGroup=gGrp,
                                          CellCount=cell_count,
                                          CellCount_Upper=cell_count_upper,
                                          Min=summary[1],
                                          FirstQu=summary[2],
                                          Median=summary[3],
                                          Mean=summary[4],
                                          ThirdQu=summary[5],
                                          Max=summary[6],
                                          row.names="Colony")
                               )
        
        }
        count <- count + 1
        Sys.sleep(time = 0.01)
        progress <- round(count/length(filelist), digits = 3)
        setTxtProgressBar(pb = pb, value = progress)
    }
    
    return(df_result)
}
