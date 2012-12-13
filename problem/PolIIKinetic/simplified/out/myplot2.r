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
Nplot = 100000

z <- matrix(scan('tr.dat', n=Nplot*11, quiet=TRUE),
            Nplot, 11, byrow=TRUE)
png("tr.png", width=.pwidth, height=.pheight)
plot(c(0,Nplot*0.1),c(-10,300),'n')
lines(z[,1], 100+30*z[,2])
for (j in seq(Nplot)){
  points(z[j,1], z[j,4], col=(z[j,3]))
}
dev.off()
