##
##
## Helper Functions
setPlateMap <- function(filename) {
    
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
    
    return(plateMap[r, c])
}

getGeneGroup <- function(gene) {
    geneGroup <- geneGroups[which(geneGroups$GeneID == gene), "GeneGroup"]
    return(geneGroup)
}

subBinsDF <- function (df, colony="ALL", filter, plateNumber=NULL, ignore=NULL) {
    if (!is.null(ignore)) {
        subset <- df[which(!(df$Colony %in% ignore)), ]
    }
    
    subset <- df[which(df$Filter == filter), ]
    return(subset)
}

