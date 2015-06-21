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
