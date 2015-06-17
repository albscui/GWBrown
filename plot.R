library(ggplot2)
library(reshape2)


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

mapMedians <- function(df, colour) {
    heatmap <- ggplot(df, aes(x = C, y = R, fill = Median, label=GeneID))
    
    heatmap + 
        theme_bw() + theme(panel.border=element_blank(),
                           panel.grid.major=element_blank(),
                           panel.grid.minor=element_blank()) +
        geom_tile() + geom_text(aes(colour=factor(GeneID)), size=2) +
        scale_colour_discrete(na.value="grey50", l=20) +
        scale_fill_gradient(low="black", high=colour) +
        scale_x_continuous(expand=c(0,0), breaks=seq(1:24)) +
        scale_y_discrete(expand=c(0,0), limits=rev(levels(df$R))) +
        ggtitle(paste("Median Values for GFP/RFP Intensity Ratios"))
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
    df$Bins <- factor(df$Bins, levels = df$Bins)
    
    plot <- ggplot(df, aes(x=Colony, y=Bins, fill=Densities))
    
    plot + geom_raster()  +
        theme(axis.text.y = element_blank()) +
        scale_fill_gradient(low="black", high=colour) +
        ggtitle("Distributions of cell densities")
}
