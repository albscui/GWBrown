binData <- function(file) {
    df_in <- read.csv(file, colClasses=c("NULL", "numeric", "numeric", "numeric"))
    filters <- c("GFP", "RFP", "GFPRFPRatio")
    for (f in filters) {
        vector <- df_in[, f]
        if (f == "GFP") {
            cutsites = c(seq(0, 10000, by=100), 80000)
        } else if (f == "RFP") {
            cutsites = c(seq(0, 5000, by = 50), 20000)
        } else if (f == "GFPRFPRatio") {
            cutsites = c(seq(0, 10, by=0.1), 800)
        }
        
        bins <- cut(vector, breaks=cutsites)
        
    }
}


# Bin data for visualization of distributions
buid_bin_df <- function(filelist) {
    df_result = data.frame()
    pb <- txtProgressBar(min = 0, max = 1, style = 3)
    count <- 0
    for (file in filelist) {
        df_in <- read.csv(file)
        for(col in df_in[]) {
            bins <- cut
        }
    }
    close(pb)
    return(df_result)
    
}
