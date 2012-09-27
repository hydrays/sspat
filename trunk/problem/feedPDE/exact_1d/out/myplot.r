library("lattice")
library("grid")

## parainfo <- read.csv("control.csv", strip.white=TRUE)
## L <- parainfo$VALUE[parainfo$PARAMETER=='L']
## H <- parainfo$VALUE[parainfo$PARAMETER=='H']
## .tend <- parainfo$VALUE[parainfo$PARAMETER=='tend']
## .tpinc <- parainfo$VALUE[parainfo$PARAMETER=='tpinc']
## .divide <- which(parainfo$PARAMETER=='useomp')
## ompinfo <- parainfo[.divide:nrow(parainfo),]
## parainfo <- parainfo[1:(.divide-1),]

N = 20000

cat("processing file ...[",N,"]\n")
i <- 0
datafile <- sprintf("%s%05d%s", "m", i, ".dat")
outfile <- sprintf("%s%05d%s", "slice", i, ".png")
png(outfile)
z <- matrix(scan(datafile, quiet=TRUE), ncol=2, byrow=TRUE)
my.label.time <- sprintf("%s%d%s", "t = ", i, " (day)")
plot(z[,1], type='b', ylim=c(0, max(z)))
lines(z[,2], type='l', col='green')

for (i in seq(N)) {
  datafile <- sprintf("%s%05d%s", "m", i, ".dat")
  outfile <- sprintf("%s%05d%s", "slice", i, ".png")
  png(outfile)
  z <- matrix(scan(datafile, quiet=TRUE), ncol=2, byrow=TRUE)
  plot(z[,1], type='b', ylim=c(0, max(z)))
  lines(z[,2], type='l', col='green')
  #my.label.time <- sprintf("%s%d%s", "t = ", as.integer(i*.tpinc), " (day)")

  output.str1 <- sprintf("%5d", i)
  if (i > 1) cat("\b\b\b\b\b")
  cat(output.str1)
  dev.off()
}
