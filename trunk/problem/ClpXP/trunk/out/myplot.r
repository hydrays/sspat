library("lattice")
library("grid")

parainfo <- read.csv("control.csv", strip.white=TRUE)

datafile <- sprintf("mean.txt")
z1 <- read.table(datafile)

datafile <- sprintf("trans.txt")
z2 <- read.table(datafile)

png('plot.png')

p1 <- xyplot(300/z1[,2]~z1[,1],
             colorkey=FALSE, xlab="Mean dwell time",
             ylab="Translocation and unfolding rate",
             xlim=c(0,5),
             ylim=c(0,1),
             type='b',
             scales=list(cex=2),
             col='red',
             panel=function(...){
                 panel.xyplot(...)
                 panel.lines(z2[,1], z2[,4], type='b', pch=2, col='blue')
             })

print(p1)

dev.off()
