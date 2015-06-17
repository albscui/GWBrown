source("D:/Home/Research/GWBrown/Project/FCS_project/FCS_scripts/GWBrown/helpers.R")
source("D:/Home/Research/GWBrown/Project/FCS_project/FCS_scripts/GWBrown/buildCellCount.R")
source("D:/Home/Research/GWBrown/Project/FCS_project/FCS_scripts/GWBrown/buildBins.R")
source("D:/Home/Research/GWBrown/Project/FCS_project/FCS_scripts/GWBrown/generate96Coord.R")

# load map of plate
plateMap <- setPlateMap("platemap_simple.csv")
# load gene groups and geneID
geneGroups <- read.csv("genegroups.csv", colClasses="character")
# generate 96 well plate locations
plate1 <- generate96Coord(1)
plate2 <- generate96Coord(2)
plate3 <- generate96Coord(3)
plate4 <- generate96Coord(4)

# list the files in the present working directory
files <- grep("[[:upper:][:digit:][:digit:]].csv$", 
              list.files(path = getwd()), value = TRUE)

temp <- buildBins(files)


