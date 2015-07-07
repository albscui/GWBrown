binData <- function(vector, f) {
    
    # cut sites depend on the vector
    # every set of cut sites amount to 101 numbers, equaling to 100 bins
    if (f == "GFP") {
        cutsites = c(seq(0,  131250, by = (131250/500)), 262500)
    } else if (f == "RFP") {
        cutsites = c(seq(0, 131050, by = (131050/500)), 262100)
    } else if (f == "GFPRFPRatio") {
        cutsites = c(seq(0, 595.5, by = (595.5/500)), 1191)
    }
    
    bins <- cut(vector, breaks=cutsites)
    levels(bins) <- cutsites[1:length(cutsites)-1]
    return(bins)
}


# Bin data for visualization of distributions
buildBins <- function(filelist) {
    df_out <- data.frame()
    pb <- txtProgressBar(min = 0, max = 1, style = 3)
    count <- 0
    
    for (file in filelist) {
        classes <- c("NULL", "numeric", "numeric", "numeric")
        df_in <- read.csv(file, colClasses=classes)
        for (colName in colnames(df_in)) {
            colony <- getName(file)
            row <- substring(colony, 1, 1)
            col <- substring(colony, 2, length(colony))
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

buildBins2 <- function(filelist, filter) {
    df_out <- data.frame()
    pb <- txtProgressBar(min = 0, max = 1, style = 3)
    count <- 0
    
    for (file in filelist) {
        classes <- c("NULL", "numeric", "numeric", "numeric")
        intensities <- read.csv(file, colClasses=classes)[, filter]
        
        colony <- getName(file)
        gID <- getGeneID(colony)
        gGrp <- getGeneGroup(gID)
        newRow <- data.frame(Colony = colony,
                             GeneID = gID,
                             GeneGroup = gGrp)
        bins <- binData(intensities, filter)
        bins_df <- data.frame(table(bins))
        bins_df$Dens <- round(bins_df$Freq/sum(bins_df$Freq), digits = 3)
        bins_df <- as.data.frame(bins_df[, 3])
        bins_df_row <- as.data.frame(t(bins_df))
        colnames(bins_df_row) <- seq(1, length(bins_df_row))
        newRow <- cbind(newRow, bins_df_row)
        row.names(newRow) <- NULL
        
        df_out <- rbind(df_out, newRow)
        
        count <- count + 1
        Sys.sleep(time = 0.01)
        progress <- round(count/length(filelist), digits = 3)
        setTxtProgressBar(pb = pb, value = progress)
    }
    close(pb)
    return(df_out)
}
