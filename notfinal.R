<<<<<<< HEAD
library(ggplot2)
library(grid)
library(gridExtra)
#
#
=======
>>>>>>> cdc5559921eac3f55b2ff6c559d00c902cb473d3


# ------------------------- Create plate map matrix ---------------------------

<<<<<<< HEAD

# Optimize the filename.
=======
plateMap <- function(directory, csvfile) {
    
    ## 'csvfile' is a comma seperated file containing the map of the plate
    ## creates a 16x24 matrix containing the identifiers for each colony
    ## each identifier is a string character
    
    setwd(directory)
    mapdata <- read.csv(csvfile, head=FALSE, colClasses="character")
    rownames <- toupper(letters[1:16])
    colnames <- as.character(seq(1:24))
    m <- matrix(nrow=16, ncol=24, dimnames=list(rownames, colnames))
    i <- 1
    for (v in mapdata) {
        m[, i] <- v
        i <- i + 1
    }
    
    return(m)
}

plateMatrix <- plateMap("map", "platemap_simple.csv") # A matrix for platemap

#-------------------------- Main Code -----------------------------------------


getName <- function(file) {
    # Optimize the filename.
    
    name <- gsub(basename(file), pattern='.csv', replacement='')
    if(substr(name, 2,2) == 0) {
        name <- paste(substr(name, 1,1), substr(name, 3, 3), sep='')
    }
    return(name)
}


getGeneID <- function(name) {
    # Extract ID from plateMatrix
    
    r <- substring(name, 1, 1)
    c <- substring(name, 2, nchar(name))
    
    return(plateMatrix[r,c])  
}


plateMap <- function(directory, csvfile) {
    
    ## 'csvfile' is a comma seperated file containing the map of the plate
    ## creates a 16x24 matrix containing the identifiers for each colony
    ## each identifier is a string character
    
    setwd(directory)
    mapdata <- read.csv(csvfile, head=FALSE, colClasses="character")
    rownames <- toupper(letters[1:16])
    colnames <- as.character(seq(1:24))
    m <- matrix(nrow=16, ncol=24, dimnames=list(rownames, colnames))
    i <- 1
    for (v in mapdata) {
        m[, i] <- v
        i <- i + 1
    }
    
    return(m)
}

plateMatrix <- plateMap("map", "platemap_simple.csv") # A matrix for platemap

#-------------------------- Main Code -----------------------------------------


# Takes in a file, lists, -> modifies the lists.
binData <- function(file, colindex) {
    # Establish variables for use inside this functions
    gfpmax <- 262143.00
    gfpmin <- 1.74
    rfpmax <- 12754.86
    rfpmin <- 0.94
    ratiomax <- 4.226582e+02
    ratiomin <- 6.540862e-03
    data <- read.csv(file, head = T, colClasses = classes)
    colname <- colnames(df)[colindex]
    data <- log2(df[, colindex])
    bins <- cut(seq(6.540862e-03, 4.226582e+02, length.out=1001), 1000)
    # Bin the data
    if(colname == "GFP") {
        
        histobj <- hist(data$GFP,
                  breaks=seq(gfpmin, gfpmax, length.out=1001),
                  plot=FALSE)
        counts <- histobj$counts
        result <- data.frame(Bins=levels(bins), Counts=counts)
        
    } else if(colname == "RFP") {
        
        histobj <- hist(data$RFP,
                        breaks=seq(rfpmin, rfpmax, length.out=1001),
                        plot=FALSE)
        counts <- histobj$counts
        result <- data.frame(Bins=levels(bins), Counts=counts)
        
    } else if(colname == "GFPRFPRatio") {
        histobj <- hist(data$GFPRFPRatio,
                        breaks=seq(ratiomin, ratiomax, length.out=1001),
                        plot=FALSE)
        counts <- histobj$counts
        result <- data.frame(Bins=levels(bins), Counts=counts)
    }
    
    return(bins)

binData <- function(file, colindex) {
    # Takes in a file, lists, -> modifies the lists.
    # Establish variables for use inside this functions
    
    gfpmax <- 262143
    gfpmin <- 0
    rfpmax <- 262144
    rfpmin <- 0
    ratiomax <- 1192
    ratiomin <- 0
    df_old <- read.csv(file, head = T)
    colname <- colnames(df_old)[colindex]
    # Bin the data
    
    if(colname == "GFP") {
        
        cutsites <- seq(gfpmin, gfpmax, length.out=1001)
        bins <- levels(cut(cutsites, 1000))
        histobj <- hist(df_old$GFP,
                  breaks=cutsites,
                  plot=FALSE)
        counts <- histobj$counts
        result <- data.frame(Bins=bins, Counts=counts)
        
    } else if(colname == "RFP") {
        
        cutsites <- seq(rfpmin, rfpmax, length.out=1001)
        bins <- levels(cut(cutsites, 1000))
        histobj <- hist(df_old$RFP,
                        breaks=cutsites,
                        plot=FALSE)
        counts <- histobj$counts
        result <- data.frame(Bins=bins, Counts=counts)
        
    } else if(colname == "GFPRFPRatio") {
        
        print(max(df_old$GFPRFPRatio))
        cutsites <- seq(ratiomin, ratiomax, length.out=1001)
        bins <- levels(cut(cutsites, 1000))
        histobj <- hist(df_old$GFPRFPRatio,
                        breaks=cutsites,
                        plot=FALSE)
        counts <- histobj$counts
        result <- data.frame(Bins=bins, Counts=counts)
    }
    
    return(result)
}

binData_log2 <- function(file, colindex) {
    # Takes in a file, lists, -> modifies the lists.
    # Establish variables for use inside this functions
    
    gfpmax <- log2(262143)
    gfpmin <- -0.221004
    rfpmax <- log2(262144)
    rfpmin <- -0.08926734*1.5
    ratiomax <- log2(1192)
    ratiomin <- -9.093109
    df_old <- read.csv(file, head = T)
    colname <- colnames(df_old)[colindex]
    # Bin the data
    
    if(colname == "GFP") {
        print(min(log2(df_old$GFP)))
        cutsites <- seq(gfpmin, gfpmax, length.out=1001)
        bins <- levels(cut(cutsites, 1000))
        histobj <- hist(log2(df_old$GFP),
                        breaks=cutsites,
                        plot=FALSE)
        counts <- histobj$counts
        result <- data.frame(Bins=bins, Counts=counts)
        
    } else if(colname == "RFP") {
        print(min(log2(df_old$RFP)))
        cutsites <- seq(rfpmin, rfpmax, length.out=1001)
        bins <- levels(cut(cutsites, 1000))
        histobj <- hist(log2(df_old$RFP),
                        breaks=cutsites,
                        plot=FALSE)
        counts <- histobj$counts
        result <- data.frame(Bins=bins, Counts=counts)
        
    } else if(colname == "GFPRFPRatio") {
        print(min(log2(df_old$GFPRFPRatio)))
        cutsites <- seq(ratiomin, ratiomax, length.out=1001)
        bins <- levels(cut(cutsites, 1000))
        histobj <- hist(log2(df_old$GFPRFPRatio),
                        breaks=cutsites,
                        plot=FALSE)
        counts <- histobj$counts
        result <- data.frame(Bins=bins, Counts=counts)
    }
    
    return(result)

}


build <- function(files, colindex) {
    # Build data frame
    
    result <- data.frame()
    i <- 0
    pb <- txtProgressBar(min=0, max=1, style=3)
    
    for(f in filelst) {
        Sys.sleep(0.1)
        progress <- round(i/length(filelst), digits=3)
        setTxtProgressBar(pb, progress)
        

        coloc <- getName(f) # Colony location on plate e.g. "A1"
        gID <- extractGeneID(getName(f)) # Gene ID e.g. "RAD51"
        binsTable <- table(binData(f, columnIndex)) # Tabulated bin data

        colony <- getName(f) # Colony location on plate e.g. "A1"
        gID <- getGeneID(getName(f)) # Gene ID e.g. "RAD51"
        binsFrame <- binData(f, colindex) # Tabulated bin data
        
        combo <- data.frame(Colony = colony,
                            GeneID = gID,
                            binsFrame) # Data frame for specific file
        
        result <- rbind(result, combo) # Append to the final data frame
        i <- i + 1
    }
    close(pb)
    
    return(result)
}

columnIndex <- menu(colnames(read.csv(csvFiles[1], head=T)), graphics=T)
final_dataframe <- build(csvFiles)

# ---------------------------- Heat Map --------------------------------------


build_log2 <- function(files, colindex) {
    # Build data frame
    
    result <- data.frame()
    i <- 0
    pb <- txtProgressBar(min=0, max=1, style=3)
    
    for(f in files) {
        Sys.sleep(0.1)
        progress <- round(i/length(files), digits=3)
        setTxtProgressBar(pb, progress)
        
        colony <- getName(f) # Colony location on plate e.g. "A1"
        gID <- getGeneID(colony) # Gene ID e.g. "RAD51"
        binsFrame <- binData_log2(f, colindex) # Tabulated bin data
        
        combo <- data.frame(Colony = colony,
                            GeneID = gID,
                            binsFrame) # Data frame for specific file
        
        result <- rbind(result, combo) # Append to the final data frame
        i <- i + 1
    }
    close(pb)
    
    return(result)
}

getMedian <- function(file, col_index) {
    colonydata <- read.csv(file, head=T)
    vector_select <- colonydata[, col_index]
    median <- median(log2(vector_select))
    return(median)
}

buildDF_medians<- function(files, col_index) {
    ignoreLst_r <- sapply(ignoreLst, substring, 1, 1)
    ignoreLst_c <- sapply(ignoreLst, substring, 2, 3)
    ignoreLst_gID <- sapply(ignoreLst, getGeneID)
    result <- data.frame(R=ignoreLst_r, C=ignoreLst_c, GeneID=ignoreLst_gID, Median=0)
    i <- 0
    pb <- txtProgressBar(min=0, max=1, style=3)
    
    for (f in files) {
        Sys.sleep(0.1)
        progress <- round(i/length(files), digits=3)
        setTxtProgressBar(pb, progress)
        
        colony <- getName(f)
        geneID <- getGeneID(colony)
        row <- substring(colony, 1, 1)
        col <- substring(colony, 2,)
        median <- getMedian(f, col_index)
        df <- data.frame(R=row, C=col, GeneID=geneID, Median=median)
        result <- rbind(result, df)
        i <- i + 1
    }
    close(pb)
    return(result)
}


# ----------------------- Practice Code Here --------------------------------

filelst <- list.files("fcsdata", pattern='.csv', full.names=TRUE)
selectcol <- menu(colnames(read.csv(filelst[1], head=T)), graphics=T)
medians_filtered <- sapply(filelst, getMedian, col_index=selectcol)
df_medians_ratios <- buildDF_medians(filelst, selectcol)
df_medians_ratios_sorted <- df_medians_ratios[order("R", "C"), ]
ignoreLst <- c("E3", "E14", "G5", "G7", "G16", "G18", "J7", "J19", "J22", "L3", "L14", "N5", "N7", "N16", "N18", "N21", "O18")
