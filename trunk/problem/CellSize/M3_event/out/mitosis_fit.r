#######################################
## Asynchronized cell size distribution
#######################################

require(seewave)
ExpResultA <- read.csv('asyn_dist.csv')
ExpResultB <- read.csv('newborn_dist.csv')

parainfo <- read.csv("control.csv", strip.white=TRUE)
NPool <- parainfo$VALUE[parainfo$PARAMETER=='NPool']

## ## Use stable distribution as a reference
## i <- 10
## datafile <- sprintf("%s%05d%s", "m", i, ".dat")
## z <- matrix(scan(datafile, n=NPool*4, quiet=TRUE),
##             NPool, 4, byrow=TRUE)
## d1 <- density(z[,1])
## ExpResultA <- spline(d1, n=101, xmin=0, xmax=3000)
## ExpResultA$y <- pmax(0, ExpResultA$y)
## ExpResultA$y <- ExpResultA$y/(sum(ExpResultA$y)*30)
## datafile <- sprintf("%s%05d%s", "n", i, ".dat")
## z <- matrix(scan(datafile, n=NCollect*1, quiet=TRUE),
##             NCollect, 1, byrow=TRUE)
## d2 <- density(z[,1])
## ExpResultB <- spline(d2, n=101, xmin=0, xmax=3000)
## ExpResultB$y <- pmax(0, ExpResultB$y)
## ExpResultB$y <- ExpResultB$y/(sum(ExpResultB$y)*30)

#for (i in seq(100)) {

i <- 5
outfile <- sprintf("%s%05d%s", "asynew", i, ".png")
png(outfile, width=600, height=800)
#pdf(outfile, width=.pwidth, height=.pheight)
#pdf("fasyn.pdf", width=7, height=5)

par(mfrow=c(2,1))

datafile <- sprintf("%s%05d%s", "m", i, ".dat")
z <- matrix(scan(datafile, n=NPool*6, quiet=TRUE),
            NPool, 6, byrow=TRUE)
d1 <- density(z[,1], from=0, to=3000)
plot(d1, xlim=c(0,3000), ylim=c(0, 0.0015))

SimResultA <- spline(d1, n=101, xmin=0, xmax=3000)
SimResultA$y <- pmax(0, SimResultA$y)
SimResultA$y <- SimResultA$y/(sum(SimResultA$y)*30)
lines(SimResultA$x, SimResultA$y, type='l', col='red')

lines(ExpResultA$x, ExpResultA$y, type='l', col='blue')
## ErrorA[i] <- norm(as.matrix(ExpResultA$y-SimResultA$y), type='O')
tt1<-cbind(ExpResultA$x, ExpResultA$y)
tt2<-cbind(SimResultA$x, SimResultA$y)
ErrorA <- kl.dist(tt1,tt2)$D1

# Newborn part
datafile <- sprintf("%s%05d%s", "n", i, ".dat")
z <- matrix(scan(datafile, n=NPool*1, quiet=TRUE),
            NPool, 1, byrow=TRUE)
d2 <- density(z[,1], from=0, to=3000)
plot(d2, xlim=c(0,3000), ylim=c(0, 0.003))

SimResultB <- spline(d2, n=101, xmin=0, xmax=3000)
SimResultB$y <- pmax(0, SimResultB$y)
SimResultB$y <- SimResultB$y/(sum(SimResultB$y)*30)
lines(SimResultB$x, SimResultB$y, type='l', col='red')

lines(ExpResultB$x, ExpResultB$y, type='l', col='blue')
## ErrorB[i] <- norm(as.matrix(ExpResultB$y-SimResultB$y), type='O')
tt1<-cbind(ExpResultB$x, ExpResultB$y)
tt2<-cbind(SimResultB$x, SimResultB$y)
ErrorB <- kl.dist(tt1,tt2)$D1

ErrorT <- ErrorA + ErrorB
dev.off()
#}
