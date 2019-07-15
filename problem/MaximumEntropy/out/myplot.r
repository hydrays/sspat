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

N = 1000
NSample = 2000
#pL = 400
#pH = 60
.pwidth = 1024
.pheight = 480

cat("processing file ...[",N,"]\n")
for (i in seq(N)) {
  datafile <- sprintf("%s%05d%s", "m", i, ".dat")
  outfile <- sprintf("%s%05d%s", "slice", i, ".png")
  z <- matrix(scan(datafile, n=2000*4, quiet=TRUE),
              2000, 4, byrow=TRUE)
  png(outfile, width=.pwidth, height=.pheight)
  plot(c(0,2000),c(0,150),'n', main=list("Controling site 3", cex=2),
       xlab = list("Cell index", cex=2),
       ylab = list("Energy level (arbitary unit)", cex = 2))
  text(1000, 130, paste("time =", i), cex = 2)
  axis(1, cex=2)
  
  for (j in seq(NSample)){
    points(j, z[j,4], col=(z[j,1]))
  }
  
  output.str1 <- sprintf("%5d", i)
  if (i > 1) cat("\b\b\b\b\b")
  cat(output.str1)
  dev.off()
}
