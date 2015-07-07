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

mapCellCounts <- function(df, fill, colour, ignore=NULL, plates=c(1,2,3,4), gene=NULL) {
    filters <- c("GFP", "RFP", "GFPRFPRatio")
    names(filters) <- c("green", "red", "yellow")
    f <- filters[colour]
    DF <- df
    DF$Column <- as.character(DF$Column)
    DF$Column <- factor(DF$Column, levels=unique(DF$Column))
    DF$Row <- as.character(DF$Row)
    DF$Row <- factor(DF$Row, levels=unique(DF$Row))
    DF <- subBinsDF(DF, filter=f, ignore=ignore, plateNumbers = plates, gene = gene)
    # Create the plot below
    heatmap <- ggplot(DF, aes_string(x = "Column", y = "Row", fill = fill))
    heatmap +
        theme_bw() + theme(panel.border=element_blank(),
                           panel.grid.major=element_blank(),
                           panel.grid.minor=element_blank(),
                           axis.title.x=element_blank(),
                           axis.title.y=element_blank()) +
        geom_tile() +
        geom_text(aes(label=GeneID, colour=factor(GeneID)), size=2) +
        guides(colour=FALSE) +
        scale_fill_gradient(low="black", high=colour) +
        scale_x_discrete(expand=c(0,0), limits=levels(DF$Column)) +
        scale_y_discrete(expand=c(0,0), limits=rev(levels(DF$Row))) +
        ggtitle(paste(fill, "plotted for ", filters[colour], "Plate:", plates))
}

plothist <- function(df, filter=NULL, by=NULL, geneGrp=NULL, gene=NULL, xlimit=NULL) {
    
    if (!is.null(geneGrp)) {
        plotobj <- ggplot(df[which(df$GeneGroup == geneGrp | df$GeneID == "HIS3"), ], 
                          aes(x=Intensity, colour=GeneID, group=GeneID))
        
    } else if (!is.null(gene)) {
        plotobj <- ggplot(df[which(df$GeneID == gene), ], 
                          aes(x=Intensity, colour=Colony, group=Colony))
        
    } else {
        plotobj <- ggplot(df, aes(x=Intensity, colour=GeneGroup, group=GeneGroup))
    }
    
    plotobj +
        stat_bin(binwidth=100, aes(y = ..density..), geom="line", position="identity") +
        
        theme_bw() +
        scale_x_continuous(limits=c(0, xlimit)) +
        xlab(paste("Intensity (", filter, ")", sep='')) +
        ylab("Density") +
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
    
    plot + theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5),
                 axis.text.y = element_blank(), axis.ticks=element_blank()) + 
        scale_x_discrete(expand=c(0,0), breaks=axis_x_breaks) +
        scale_y_discrete(expand=c(0,0), breaks=axis_y_breaks) + 
        ggtitle(paste("Distributions of cell densities for", filters[colour]))
}


plotScatter <- function(df) {
    ggplot(df, aes(x=Colony, y=CellCount)) +
        geom_point() +
        scale_x_discrete(expand=c(0,0), limits=c(plate1,plate2,plate3,plate4))
}


plotVariance <- function(df, filter="GFP") {
    colours <- c("green", "red", "yellow")
    names(colours) <- c("GFP", "RFP", "GFP/RFP Ratio")
    plot <- ggplot(data = df, aes(nSample, SE))
    plot + 
        geom_point() +
        geom_hline(aes(yintercept = median(SE)),
                   alpha=0.5,
                   size=2, 
                   colour=colours[filter]) +
        theme_bw() +
        labs(title = paste("Vairance based on sample size", filter),
             x = "Number of Samples")
}




plothist(dataGFP, filter="GFP", by="Gene Group", xlimit=50000)

# Plot cell count against intensity
ggplot(countsData[which(countsData$Filter == c("GFP", "RFP")), ], aes(x=CellCount, y=Median, colour=Filter)) + 
    geom_point() + 
    geom_smooth(method=lm) + 
    facet_grid(. ~ Filter) +
    theme_bw() +
    scale_y_continuous(limits=c(0, 3000)) +
    scale_colour_manual(values=c("#00CC00", "#CC0000"))

cor(x = temp$CellCount, y = temp$Median)
cor(x = counts_rfp$CellCount, y = counts_rfp$Median)


