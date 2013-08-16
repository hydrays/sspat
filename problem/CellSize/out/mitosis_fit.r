#######################################
## Asynchronized cell size distribution
#######################################
#outfile <- sprintf("%s%05d%s", "fa", i, ".png")
#png(outfile, width=.pwidth, height=.pheight)
#pdf(outfile, width=.pwidth, height=.pheight)
#pdf("fasyn.pdf", width=7, height=5)

par(mfrow=c(2,1))

i <- 5000
datafile <- sprintf("%s%05d%s", "m", i, ".dat")
z <- matrix(scan(datafile, n=NPool*L, quiet=TRUE),
            NPool, L, byrow=TRUE)
d1 <- density(z[,1])
plot(d1, xlim=c(0,3000), ylim=c(0, 0.0015))

SimResultA <- spline(d1, n=101, xmin=0, xmax=3000)
SimResultA$y <- pmax(0, SimResultA$y)
SimResultA$y <- SimResultA$y/(sum(SimResultA$y)*30)
lines(SimResultA$x, SimResultA$y, type='l', col='red')

ExpResultA <- read.csv('asyn_dist.csv')
lines(ExpResultA$x, ExpResultA$y, type='l', col='blue')
ErrorA <- norm(as.matrix(ExpResultA$y-SimResultA$y))

# Newborn part
datafile <- sprintf("CellNewborn.dat")
z <- matrix(scan(datafile, n=NPool*L, quiet=TRUE),
            NPool, L, byrow=TRUE)
d2 <- density(z[,1])
plot(d2, xlim=c(0,3000), ylim=c(0, 0.003))

SimResultB <- spline(d2, n=101, xmin=0, xmax=3000)
SimResultB$y <- pmax(0, SimResultB$y)
SimResultB$y <- SimResultB$y/(sum(SimResultB$y)*30)
lines(SimResultB$x, SimResultB$y, type='l', col='red')

ExpResultB <- read.csv('newborn_dist.csv')
lines(ExpResultB$x, ExpResultB$y, type='l', col='blue')
ErrorB <- norm(as.matrix(ExpResultB$y-SimResultB$y))

ErrorT <- ErrorA + ErrorB
