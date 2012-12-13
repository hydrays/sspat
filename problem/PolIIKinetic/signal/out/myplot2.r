library("lattice")
library("grid")
jet.colors <- colorRampPalette(c("white", "red", "blue", "green"))

#parainfo <- read.csv("control.csv", strip.white=TRUE)
#L <- parainfo$VALUE[parainfo$PARAMETER=='L']
#H <- parainfo$VALUE[parainfo$PARAMETER=='H']
#.tend <- parainfo$VALUE[parainfo$PARAMETER=='tend']
#.tpinc <- parainfo$VALUE[parainfo$PARAMETER=='tpinc']
#.divide <- which(parainfo$PARAMETER=='useomp')
#ompinfo <- parainfo[.divide:nrow(parainfo),]
#parainfo <- parainfo[1:(.divide-1),]

NSample = 2000
.pwidth = 2048
.pheight = 576
Nplot = 10000

z <- matrix(scan('tr.dat', n=Nplot*3, quiet=TRUE),
            Nplot, 3, byrow=TRUE)
png("signal.png", width=.pwidth, height=.pheight)
plot(z[,1], z[,3]/max(z[,3]), col='red', type='b')
lines(z[,1], z[,2])
dev.off()
