binData <- function(vector, f) {
    
    # cut sites depend on the vector
    # every set of cut sites amount to 101 numbers, equaling to 100 bins
    if (f == "GFP") {
        cutsites = c(seq(0, 10000, by=50.2), 80000)
    } else if (f == "RFP") {
        cutsites = c(seq(0, 5000, by= 25.1), 20000)
    } else if (f == "GFPRFPRatio") {
        cutsites = c(seq(0, 10, by=0.0502), 800)
    }
    
    bins <- cut(vector, breaks=cutsites)
    return(bins)
}


# Bin data for visualization of distributions
buildBins <- function(filelist) {
    df_out = data.frame()
    pb <- txtProgressBar(min = 0, max = 1, style = 3)
    count <- 0
    
    for (file in filelist) {
        classes <- c("NULL", "numeric", "numeric", "numeric")
        df_in <- read.csv(file, colClasses=classes)
        for(colName in colnames(df_in)) {
            colony <- getName(file)
            row <- substring(colony, 1, 1)
            col <- substring(colony, 2,  )
            gID <- getGeneID(colony)
            gGrp <- getGeneGroup(gID)

            Bins <- binData(df_in[, colName], colName)
            Bins_df <- data.frame(table(Bins))
            total <- sum(Bins_df$Freq)
            densities <- data.frame(Densities=round(Bins_df$Freq/total, digits = 5))
            Bins_df <- cbind(Bins_df, densities)
            
            newrows <- cbind(Bins_df,
                             data.frame(Colony=colony,
                                        GeneID=gID,
                                        GeneGroup=gGrp,
                                        Filter=colName)
                             )
            df_out <- rbind(df_out, newrows)
            
        }
        
        count <- count + 1
        Sys.sleep(time = 0.01)
        progress <- round(count/length(filelist), digits = 3)
        setTxtProgressBar(pb = pb, value = progress)
    }
    colnames(df_out)[2] <- "Bins_freq"
    close(pb)
    return(df_out)
    
}
