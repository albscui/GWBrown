library(ggplot2)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    # Multiple plot function
    #
    # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
    # - cols:   Number of columns in layout
    # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
    #
    # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
    # then plot 1 will go in the upper left, 2 will go in the upper right, and
    # 3 will go all the way across the bottom.
    #
    
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}

mapCellCounts <- function(df, colour) {
    filters <- c("GFP", "RFP", "GFP/RFP Ratios")
    names(filters) <- c("green", "red", "yellow")
    
    heatmap <- ggplot(df, aes(x=Column, y=Row, fill = Median))
    
    heatmap +
        theme_bw() + theme(panel.border=element_blank(),
                           panel.grid.major=element_blank(),
                           panel.grid.minor=element_blank(),
                           axis.text.x=element_text(size=12),
                           axis.text.y=element_text(size=12)) +
        geom_tile() +
        geom_text(aes(label=GeneID, colour=factor(GeneID)), size=2) +
        guides(colour=FALSE) +
        scale_fill_gradient(low="black", high=colour) +
        scale_x_continuous(expand=c(0,0), breaks=seq(1:24)) +
        scale_y_discrete(expand=c(0,0), limits=rev(levels(df$R))) +
        ggtitle(paste("Median for ", filters[colour]))
}

plothist <- function(df, filter=NULL, by=NULL, geneGrp=NULL, gene=NULL, xlimit=NULL) {
    
    if (!is.null(geneGrp)) {
        plotobj <- ggplot(df[which(df$geneGroup == geneGrp | df$geneID == "HIS3"), ], 
                          aes(x=Intensity, colour=geneID))
        
    } else if (!is.null(gene)) {
        plotobj <- ggplot(df[which(df$geneID == gene), ], 
                          aes(x=Intensity, colour=colonyLocation))
        
    } else {
        plotobj <- ggplot(df, aes(x=Intensity, colour=geneGroup))
    }
    
    plotobj + geom_density() + theme_bw() +
        scale_x_continuous(limits=c(0, xlimit)) +
        xlab(paste("Intensity (", filter, ")", sep='')) +
        ylab("Density of Cell Counts") +
        ggtitle(paste("Distributions of RNR3 Expression by ", by, sep=''))
}



plotHeatMap <- function(df, colour) {
    filters <- c("GFP", "RFP", "GFP/RFP Ratios")
    names(filters) <- c("green", "red", "yellow")

    df$Bins <- as.character(df$Bins)
    df$Bins <- factor(df$Bins, levels=unique(df$Bins))
    df$Colony <- as.character(df$Colony)
    df$Colony <- factor(df$Colony, levels=unique(df$Colony))

    plot <- ggplot(df, aes(x=Colony, y=Bins, fill=Densities)) + geom_raster() +
        scale_fill_gradient(low="black", high=colour)

    axis_x_breaks <- levels(df$Colony)[seq(1, length(levels(df$Colony)), by=5)]
    axis_y_breaks <- levels(df$Bins)[seq(1, 201, by=1)]
    
    return(plot + theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5),
                 axis.text.y = element_blank(), axis.ticks=element_blank()) + 
        scale_x_discrete(expand=c(0,0), breaks=axis_x_breaks) +
        scale_y_discrete(expand=c(0,0), breaks=axis_y_breaks) + 
        ggtitle(paste("Distributions of cell densities for", filters[colour])))
}

gfp <- subBinsDF(bindata, filter="GFP", ignore=ignore, plateNumbers = 4)
rfp <- subBinsDF(bindata, filter="RFP", ignore=ignore, plateNumbers = 4)
ratio <- subBinsDF(bindata, filter="GFPRFPRatio", ignore=ignore, plateNumbers = 4)

multiplot(plotHeatMap(gfp, "green"),
          plotHeatMap(rfp, "red"), 
          plotHeatMap(ratio, "yellow"), cols=1)
