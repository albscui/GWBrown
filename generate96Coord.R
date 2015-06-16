generate96Coord <- function (plateNumber) {
    lettersA <- sort(rep(toupper(letters[seq(1, 16, 2)]), times=12))
    lettersB <- sort(rep(toupper(letters[seq(2, 16, 2)]), times=12))
    numbers1 <- rep(seq(1, 24, 2), times=8)
    numbers2 <- rep(seq(2, 24, 2), times=8)
    
    if (plateNumber == 1) {
        return(paste(lettersA, numbers1, sep=""))
    } else if (plateNumber == 2 ) {
        return(paste(lettersA, numbers2, sep=""))
    } else if (plateNumber == 3) {
        return(paste(lettersB, numbers1, sep=""))
    } else if (plateNumber == 4) {
        return(paste(lettersB, numbers2, sep=""))
    } else {
        print("Please provide a plate number (integer from 1 to 4 as a parameter.")
    }
}
