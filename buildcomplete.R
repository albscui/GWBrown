## Analysis of Flow Cytometry data Functions generate final dataframe
## for ploting and statistical analysis


## 



## Helper Functions
plateMap <- function(filename) {
    
    ## 'csvfile' is a comma seperated file containing the map of the plate
    ## creates a 16x24 matrix containing the identifiers for each colony
    ## each identifier is a string character
    
    mapdata <- read.csv(filename, head = FALSE, colClasses = "character")
    rownames <- toupper(letters[1:16])
    colnames <- as.character(seq(1:24))
    m <- matrix(nrow = 16, ncol = 24, dimnames = list(rownames, colnames))
    i <- 1
    
    for (v in mapdata) {
        m[, i] <- v
        i <- i + 1
    }
    
    return(m)
}


getName <- function(file) {
    # Optimize the filename.
    
    name <- gsub(basename(file), pattern = ".csv", replacement = "")
    if (substr(name, 2, 2) == 0) {
        name <- paste(substr(name, 1, 1), substr(name, 3, 3), sep = "")
    }
    return(name)
}


getGeneID <- function(name) {
    # Extract ID from plateMatrix
    
    r <- substring(name, 1, 1)
    c <- substring(name, 2, nchar(name))
    
    return(plateMatrix[r, c])
}

# Main
# Function----------------------------------------------------------------

build_df_complete <- function(directory, filter = NULL, outputfile = NULL) {
    sourcedir <- paste(getwd())  # Remember source file directory
    
    # Retrieve list of full file directory names for reading
    setwd(paste(getwd(), "/", directory, sep = ""))
    fileList <- grep("[[:upper:][:digit:][:digit:]].csv$",
                     list.files(path = getwd()), value = TRUE)
    # Initialize empty dataframe for merging with extracted data
    df_final <- data.frame()
    count <- 1  # set count to 1
    
    # Progress bar
    pb <- txtProgressBar(min = 0, max = 1, style = 3)
    
    for (f in fileList) {
        # f is a filename in alist of files df is a dataframe initialized by
        # reading from file f
        classes <- sapply(read.csv(f, nrows = 10), class)
        df <- read.csv(f, colClasses = classes)
        
        # Column ID and corresponding values for final data frame
        plateNumb <- NULL
        colony <- getName(f)
        gID <- getGeneID(colony)
        gGrp <- geneGroup[which(geneGroup$GeneID == gID), "GeneGroup"]
        
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
        # Progress bar specs
        Sys.sleep(0.01)
        progress <- round(count/length(fileList), digits = 3)
        setTxtProgressBar(pb, progress)
        count <- count + 1
    }
    close(pb)
    setwd(sourcedir)
    return(df_final)
}

# Make a dataframe with cell counts for heatmap
build_cell_count_df <- function(filelist) {
    # empty data frame for storing values
    df_result <- data.frame()
    
    for (file in filelist) {
        for (filter in c("GFP", "RFP", "GFPRFPRatio")) {
            df <- read.csv(file, colClasses=classes)
            colony <- getName(file)
            row <- substring(colony, 1, 1)
            col <- substring(colony, 2,  )
            gID <- getGeneID(colony)
            cell_count <- nrow(df)
            if (filter == "GFP") {
                cell_count_upper <- nrow(df[which(df$GFP > 50000), ])
            } else if (filter == "RFP") {
                cell_count_upper <- nrow(df[which(df$RFP > 5000), ])
            }
        }
    }
}


# Execute Code below -----------------------------------------------------

# Set working directory
setwd(choose.dir(getwd(), "choose wd"))
# Load plate map as matrix
plateMatrix <- plateMap("platemap_simple.csv")
# Load data frame for gene groups
geneGroup <- read.csv("genegroups.csv", stringsAsFactors = F)
# Specify the classes for combined dataframe
classes <- c("NULL", "factor", "factor", "factor", "numeric", "factor")
