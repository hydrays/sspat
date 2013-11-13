library("lattice")
library("grid")

parainfo <- read.csv("control.csv", strip.white=TRUE)

datafile <- sprintf("mean.txt")
z <- read.table(datafile)
p1 <- xyplot(z[,2]~z[,1],
                colorkey=FALSE, xlab="",
                ylab="",
                scales=list(cex=2))

print(p1)


