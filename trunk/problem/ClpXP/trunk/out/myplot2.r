library("lattice")
library("grid")

parainfo <- read.csv("control.csv", strip.white=TRUE)

datafile <- sprintf("sample.txt")
z <- read.table(datafile)
p2 <- histogram(z[,1], nint = 30,
                colorkey=FALSE, xlab="",
                ylab="",
                scales=list(cex=2))

print(p2)


