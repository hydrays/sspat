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
plot(c(0,Nplot*0.1),c(0,600),'n')
lines(z[,1], 300+30*z[,2])
for (j in seq(Nplot)){
  lines(c(z[j,1], z[j,1]), c(0, z[j,9])/10, col = 1)
  lines(c(z[j,1], z[j,1]), c(z[j,9], z[j,9]+z[j,10])/10, col = 2)
  lines(c(z[j,1], z[j,1]), c(z[j,9]+z[j,10], z[j,9]+z[j,10]+z[j,11])/10, col = 3)    
}
dev.off()
