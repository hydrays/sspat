require('lattice')
require('grid')
#######################################
## Asynchronized cell size distribution
#######################################
ExpResultA <- read.csv('asyn_dist.csv')
ExpResultB <- read.csv('newborn_dist.csv')

parainfo <- read.csv("out_sizetime/control.csv", strip.white=TRUE)
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

i <- 1
outfile <- sprintf("%s%05d%s", "asynew", i, ".png")
#png(outfile, width=800, height=600)
#pdf(outfile, width=.pwidth, height=.pheight)
pdf("fig_dists_sizetime.pdf", width=7, height=3)

datafile <- sprintf("%s%05d%s", "out_sizetime/m", i, ".dat")
z <- matrix(scan(datafile, n=NPool*6, quiet=TRUE),
            NPool, 6, byrow=TRUE)
d1 <- density(z[,1], from=0, to=3000)
SimResultA <- spline(d1, n=101, xmin=0, xmax=3000)
SimResultA$y <- pmax(0, SimResultA$y)
SimResultA$y <- SimResultA$y/(sum(SimResultA$y)*30)
trellis.par.set(clip=list(panel = "off"))
y.lim = c(0, 1.5e-3)
y.at <- pretty(y.lim)
y.labels <- formatC(1000*y.at, format = "g")
p1 <- xyplot(d1$y~d1$x,
             xlim=c(0,3000),
             ylim=c(0, 0.0015),
             lwd=2,
             type='l',
             col='black',
             xlab=list('s (cell size, fl)', cex=1),
             ylab=list('density frequency', cex = 1),
             scales=list(cex=1, y=list(at = y.at,
                                      labels = y.labels), tck=0.5),
             panel=function(...){
                 panel.lines(SimResultA$y~SimResultA$x, lwd=2,
                             type='l', col='red')
                 panel.lines(ExpResultA$y~ExpResultA$x, lwd=2,
                             type='l', col='blue')                 
                  panel.axis(side = c("top"),
                             at = c(100),
                             labels = c(expression(x10^{-3})),
                             ticks=FALSE,
                             outside=TRUE,
                             rot=c(0,90),
                             text.cex=1)
             }
             )
print(p1, position=c(0, 0, 0.5, 1), more=TRUE)
#lines(SimResultA$x, SimResultA$y, type='l', col='red')
#lines(ExpResultA$x, ExpResultA$y, type='l', col='blue')
#ErrorA[i] <- norm(as.matrix(ExpResultA$y-SimResultA$y), type='O')
ErrorA <- norm(as.matrix(ExpResultA$y-SimResultA$y), type='O')
## tt1<-cbind(ExpResultA$x, ExpResultA$y)
## tt2<-cbind(SimResultA$x, SimResultA$y)
## ErrorA[i] <- kl.dist(tt1,tt2)$D1

# Newborn part
datafile <- sprintf("%s%05d%s", "out_sizetime/n", i, ".dat")
z <- matrix(scan(datafile, n=NPool*1, quiet=TRUE),
            NPool, 1, byrow=TRUE)
d2 <- density(z[,1], from=0, to=3000)

SimResultB <- spline(d2, n=101, xmin=0, xmax=3000)
SimResultB$y <- pmax(0, SimResultB$y)
SimResultB$y <- SimResultB$y/(sum(SimResultB$y)*30)
##lines(SimResultB$x, SimResultB$y, type='l', col='red')
y.lim = c(0, 3e-3)
y.at <- pretty(y.lim, n=4)
y.labels <- formatC(1000*y.at, format = "g")
p2 <- xyplot(d1$y~d1$x,
             xlim=c(0,3000),
             ylim=c(0, 0.003),
             lwd=2,
             type='l',
             col='black',
             xlab=list('s (cell size, fl)', cex=1),
             ##ylab=list('density frequency', cex = 1),
             ylab='',
             scales=list(cex=1, y=list(at = y.at,
                                      labels = y.labels), tck=0.5),
             panel=function(...){
                 panel.lines(SimResultB$y~SimResultB$x, lwd=2,
                             type='l', col='red')
                 panel.lines(ExpResultB$y~ExpResultB$x, lwd=2,
                             type='l', col='blue')                 
                  panel.axis(side = c("top"),
                             at = c(100),
                             labels = c(expression(x10^{-3})),
                             ticks=FALSE,
                             outside=TRUE,
                             rot=c(0,90),
                             text.cex=1)
             }
             )
print(p2, position=c(0.5, 0, 1, 1))
## lines(ExpResultB$x, ExpResultB$y, type='l', col='blue')
## ##ErrorB[i] <- norm(as.matrix(ExpResultB$y-SimResultB$y), type='O')
## ErrorB <- norm(as.matrix(ExpResultB$y-SimResultB$y), type='O')
## ## tt1<-cbind(ExpResultB$x, ExpResultB$y)
## ## tt2<-cbind(SimResultB$x, SimResultB$y)
## ## ErrorB[i] <- kl.dist(tt1,tt2)$D1

## ErrorT <- ErrorA + ErrorB
dev.off()
#}
