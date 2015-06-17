source("D:/Home/Research/GWBrown/Project/FCS_project/FCS_scripts/GWBrown/helpers.R")
source("D:/Home/Research/GWBrown/Project/FCS_project/FCS_scripts/GWBrown/buildCellCount.R")

# load map of plate
plateMap <- setPlateMap("platemap_simple.csv")
# load gene groups and geneID
geneGroups <- read.csv("genegroups.csv", colClasses="character")


# list the files in the present working directory
files <- grep("[[:upper:][:digit:][:digit:]].csv$", 
              list.files(path = getwd()), value = TRUE)

cellcounts <- buildCellCount(files)
