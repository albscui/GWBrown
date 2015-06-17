build_df_complete <- function(directory, filter = NULL, outputfile = NULL) {
    sourcedir <- paste(getwd())  # Remember source file directory
    
    # Retrieve list of full file directory names for reading
    setwd(paste(getwd(), "/", directory, sep = ""))
    fileList <- grep("[[:upper:][:digit:][:digit:]].csv$",
                     list.files(path = getwd()), value = TRUE)
    # Initialize empty dataframe for merging with extracted data
    df_final <- data.frame()
    
    # Progress bar
    pb <- txtProgressBar(min = 0, max = 1, style = 3)
    count <- 0  # set count to 1
    
    for (f in fileList) {
        # Progress bar specs
        Sys.sleep(0.01)
        progress <- round(count/length(fileList), digits = 3)
        setTxtProgressBar(pb, progress)
        
        # f is a filename in alist of files df is a dataframe initialized by
        # reading from file f
        classes <- sapply(read.csv(f, nrows = 10), class)
        df <- read.csv(f, colClasses = classes)
        
        # Column ID and corresponding values for final data frame
        plateNumb <- NULL
        colony <- getName(f)
        gID <- getGeneID(colony)
        gGrp <- geneGroups[which(geneGroups$GeneID == gID), "GeneGroup"]
        
        # Data vectors extracted from data source files.  Outliers are trimmed
        # according to specs below.
        if (!is.null(filter)) {
            datavector <- df[which(df[, filter] <= median(df[, filter]) * 
                                       4), filter]
            df_final <- rbind(df_final, data.frame(colonyLocation = colony,
                                                   geneID = gID,
                                                   geneGroup = gGrp,
                                                   Intensity = datavector,
                                                   Filter = filter))
        } else {
            # Generating the final data frame using the rbind function
            df_final <- rbind(df_final, 
                              data.frame(colonyLocation = colony, 
                                         geneID = gID, 
                                         geneGroup = gGrp, 
                                         Intensity = vgfp, 
                                         Filter = "GFP"), 
                              data.frame(colonyLocation = colony, 
                                         geneID = gID, 
                                         geneGroup = gGrp, 
                                         Intensity = vrfp, 
                                         Filter = "RFP"),
                              data.frame(colonyLocation = colony, 
                                         geneID = gID, 
                                         geneGroup = gGrp, 
                                         Intensity = vratio, 
                                         Filter = "Ratio"))
        }
        
        count <- count + 1
    }
    close(pb)
    setwd(sourcedir)
    return(df_final)
}
