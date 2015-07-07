source("D:/Home/Research/GWBrown/Project/FCS_project/FCS_scripts/GWBrown/helpers.R")


# load map of plate
plateMap <- setPlateMap("platemap_simple.csv")
# load gene groups and geneID
geneList <- read.csv("geneList.csv", colClasses="character")
# generate 96 well plate locations
plate1 <- generate96Coord(1)
plate2 <- generate96Coord(2)
plate3 <- generate96Coord(3)
plate4 <- generate96Coord(4)

# list of colonies to ignore
ignore <- read.table("ignore.txt")[, 1]

# list the files in the present working directory
files <- grep("[[:upper:][:digit:][:digit:]].csv$", list.files(path = getwd()), value = TRUE)



temp <- function() {
    rows <- toupper(letters[1:16])
    cols <- as.character(seq(1,24))
    
    result <- data.frame()
    for (r in rows) {
        for (c in cols) {
            gID <- plateMap[r, c]
            print(gID)
            gORF <- getGeneORF(gID)
            print(gORF)
            newRow <- data.frame(Row = r, Column = c, GeneID = gID, ORF = gORF)
            result <- rbind(result, newRow)
        }
    }
    
    return(result)
}

gfp <- buildVariance(Filter = "GFP", rep=1000)
rfp <- buildVariance(Filter = "RFP", rep=1000)
ratios <- buildVariance(Filter = "GFPRFPRatio", rep=1000)
