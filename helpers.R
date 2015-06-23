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

subBinsDF <- function (df, filter=NULL, ignore=NULL, 
                       colony=NULL, plateNumbers=NULL,
                       gene=NULL) {
    
    subset <- NULL
    
    if (!is.null(ignore) & is.null(subset)) {
        if (is.null(subset)) {
            subset <- df[which(!(df$Colony %in% ignore)), ]
        } else {
            subset <- subset[which(!(subset$Colony %in% ignore)), ]
        }
    }
    
    if (!is.null(plateNumbers)) {
        coordinates <- generate96Coord(plateNumbers)
        if (is.null(subset)) {
            subset <- df[which(df$Colony %in% coordinates), ]
        } else {
            subset <- subset[which(subset$Colony %in% coordinates), ]
        }
    }
    
    if (!is.null(colony)) {
        if (is.null(subset)) {
            subset <- df[which(df$Colony %in% colony), ]
        } else {
            subset <- subset[which(subset$Colony %in% colony), ]
        }
        
    }
    
    if (!is.null(filter)) {
        if (is.null(subset)) {
            subset <- df[which(df$Filter == filter), ]
        } else {
            subset <- subset[which(subset$Filter == filter), ]
        }
        
    }
    
    if (!is.null(gene)) {
        if (is.null(subset)) {
            subset <- df[which(df$GeneID == gene), ]
        } else {
            subset <- subset[which(subset$GeneID == gene), ]
        }
        
    }
    
    subset <- droplevels(subset)
    return(subset)
}

generate96Coord <- function (plateNumbers) {
    lettersA <- sort(rep(toupper(letters[seq(1, 16, 2)]), times=12))
    lettersB <- sort(rep(toupper(letters[seq(2, 16, 2)]), times=12))
    numbers1 <- rep(seq(1, 24, 2), times=8)
    numbers2 <- rep(seq(2, 24, 2), times=8)
    
    coordinates <- c()
    for (n in plateNumbers)
        if (n == 1) {
            coordinates <- c(coordinates, paste(lettersA, numbers1, sep=""))
        } else if (n == 2 ) {
            coordinates <- c(coordinates, paste(lettersA, numbers2, sep=""))
        } else if (n == 3) {
            coordinates <- c(coordinates, paste(lettersB, numbers1, sep=""))
        } else if (n == 4) {
            coordinates <- c(coordinates, paste(lettersB, numbers2, sep=""))
        } else {
            print("Please provide a plate number (integer from 1 to 4 as a parameter.")
        }
    return(coordinates)
}

