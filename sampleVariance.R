
buildVariance <- function(geneName="HIS3", Filter="GFP", N=300, rep=10) {
    
    # Set seed so that results are reproducible
    set.seed(1234)
    temp <- read.csv("all_cell_counts.csv", nrows=2)
    classes <- sapply(temp, class)
    classes[c(1,4,5,8,9,10,11,12,14,15)] <- "NULL"
    data <- read.csv("all_cell_counts.csv", colClasses = classes)
    data <- subBinsDF(data, gene = geneName, filter = Filter, ignore=ignore)
    
    result <- data.frame()
    
    pb <- txtProgressBar(min = 0, max = 1, style = 3)
    counter <- 0
    for (i in seq(1, rep)) {
        for(n in seq(3, N-1)) {
            nSample <- n
            randInts <- sample(1:N, n)
            print(randInts)
            variance <- sd(data[randInts, "Mean"])
            newrow <- data.frame(Rep = i, nSample = n, SE = variance)
            result <- rbind(result, newrow)
            # For progress bar
            counter <- counter + 1
            Sys.sleep(time = 0.01)
            outerLoop <- length(seq(1, rep))
            innerLoop <- length(seq(3, N-1))
            loopTotal <- outerLoop * innerLoop
            progress <- round(counter/loopTotal, digits = 3)
            setTxtProgressBar(pb = pb, value = progress)
        }
    }
    close(pb)
    return(result)
}
