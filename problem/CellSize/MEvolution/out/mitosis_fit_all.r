require('lattice')
require('grid')
#######################################
## Asynchronized cell size distribution
#######################################
i <- 1
pdf("fig_dists_all.pdf", width=7, height=3)
trellis.par.set(clip=list(panel = "off"))

ExpResultA <- read.csv('asyn_dist.csv')
ExpResultB <- read.csv('newborn_dist.csv')

parainfo <- read.csv("out_timeonly/control.csv", strip.white=TRUE)
NPool <- parainfo$VALUE[parainfo$PARAMETER=='NPool']
datafile <- sprintf("%s%05d%s", "out_timeonly/m", i, ".dat")
z <- matrix(scan(datafile, n=NPool*6, quiet=TRUE),
            NPool, 6, byrow=TRUE)
td1 <- density(z[,1], from=0, to=3000)
tSimResultA <- spline(td1, n=101, xmin=0, xmax=3000)
tSimResultA$y <- pmax(0, tSimResultA$y)
tSimResultA$y <- tSimResultA$y/(sum(tSimResultA$y)*30)

parainfo <- read.csv("out_sizetime/control.csv", strip.white=TRUE)
NPool <- parainfo$VALUE[parainfo$PARAMETER=='NPool']
datafile <- sprintf("%s%05d%s", "out_sizetime/m", i, ".dat")
z <- matrix(scan(datafile, n=NPool*6, quiet=TRUE),
            NPool, 6, byrow=TRUE)
sd1 <- density(z[,1], from=0, to=3000)
sSimResultA <- spline(sd1, n=101, xmin=0, xmax=3000)
sSimResultA$y <- pmax(0, sSimResultA$y)
sSimResultA$y <- sSimResultA$y/(sum(sSimResultA$y)*30)

parainfo <- read.csv("out_event/control.csv", strip.white=TRUE)
NPool <- parainfo$VALUE[parainfo$PARAMETER=='NPool']
datafile <- sprintf("%s%05d%s", "out_event/m", i, ".dat")
z <- matrix(scan(datafile, n=NPool*6, quiet=TRUE),
            NPool, 6, byrow=TRUE)
ed1 <- density(z[,1], from=0, to=3000)
eSimResultA <- spline(ed1, n=101, xmin=0, xmax=3000)
eSimResultA$y <- pmax(0, eSimResultA$y)
eSimResultA$y <- eSimResultA$y/(sum(eSimResultA$y)*30)

y.lim = c(0, 1.5e-3)
y.at <- pretty(y.lim)
y.labels <- formatC(1000*y.at, format = "g")
p1 <- xyplot(ExpResultA$y~ExpResultA$x,
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
                 panel.lines(...)
                 panel.lines(eSimResultA$y~eSimResultA$x, lwd=1.5,
                             type='l', lty=1, col='blue')
                 panel.lines(tSimResultA$y~tSimResultA$x, lwd=1.5,
                             type='l', lty=2, col='brown')
                 panel.lines(sSimResultA$y~sSimResultA$x, lwd=1.5,
                             type='l', lty=4, col='red')                 
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

## Newborn
parainfo <- read.csv("out_timeonly/control.csv", strip.white=TRUE)
NPool <- parainfo$VALUE[parainfo$PARAMETER=='NPool']
datafile <- sprintf("%s%05d%s", "out_timeonly/n", i, ".dat")
z <- matrix(scan(datafile, n=NPool*1, quiet=TRUE),
            NPool, 1, byrow=TRUE)
td1 <- density(z[,1], from=0, to=3000)
tSimResultA <- spline(td1, n=101, xmin=0, xmax=3000)
tSimResultA$y <- pmax(0, tSimResultA$y)
tSimResultA$y <- tSimResultA$y/(sum(tSimResultA$y)*30)

parainfo <- read.csv("out_sizetime/control.csv", strip.white=TRUE)
NPool <- parainfo$VALUE[parainfo$PARAMETER=='NPool']
datafile <- sprintf("%s%05d%s", "out_sizetime/n", i, ".dat")
z <- matrix(scan(datafile, n=NPool*1, quiet=TRUE),
            NPool, 1, byrow=TRUE)
sd1 <- density(z[,1], from=0, to=3000)
sSimResultA <- spline(sd1, n=101, xmin=0, xmax=3000)
sSimResultA$y <- pmax(0, sSimResultA$y)
sSimResultA$y <- sSimResultA$y/(sum(sSimResultA$y)*30)

parainfo <- read.csv("out_event/control.csv", strip.white=TRUE)
NPool <- parainfo$VALUE[parainfo$PARAMETER=='NPool']
datafile <- sprintf("%s%05d%s", "out_event/n", i, ".dat")
z <- matrix(scan(datafile, n=NPool*1, quiet=TRUE),
            NPool, 1, byrow=TRUE)
ed1 <- density(z[,1], from=0, to=3000)
eSimResultA <- spline(ed1, n=101, xmin=0, xmax=3000)
eSimResultA$y <- pmax(0, eSimResultA$y)
eSimResultA$y <- eSimResultA$y/(sum(eSimResultA$y)*30)

y.lim = c(0, 3e-3)
y.at <- pretty(y.lim, n=4)
y.labels <- formatC(1000*y.at, format = "g")
p2 <- xyplot(ExpResultB$y~ExpResultB$x,
             xlim=c(0,3000),
             ylim=c(0, 0.003),
             lwd=2,
             type='l',
             col='black',
             xlab=list('s (cell size, fl)', cex=1),
             ylab='', #list('density frequency', cex = 1),
             scales=list(cex=1, y=list(at = y.at,
                                      labels = y.labels), tck=0.5),
             panel=function(...){
                 panel.lines(...)
                 panel.lines(eSimResultA$y~eSimResultA$x, lwd=1.5,
                             type='l', lty=1, col='blue')
                 panel.lines(tSimResultA$y~tSimResultA$x, lwd=1.5,
                             type='l', lty=2, col='brown')
                 panel.lines(sSimResultA$y~sSimResultA$x, lwd=1.5,
                             type='l', lty=4, col='red')                 
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
#lines(SimResultA$x, SimResultA$y, type='l', col='red')
#lines(ExpResultA$x, ExpResultA$y, type='l', col='blue')
#ErrorA[i] <- norm(as.matrix(ExpResultA$y-SimResultA$y), type='O')
ErrorA <- norm(as.matrix(ExpResultA$y-SimResultA$y), type='O')
## tt1<-cbind(ExpResultA$x, ExpResultA$y)
## tt2<-cbind(SimResultA$x, SimResultA$y)
## ErrorA[i] <- kl.dist(tt1,tt2)$D1

## ## Newborn part
## datafile <- sprintf("%s%05d%s", "out_event/n", i, ".dat")
## z <- matrix(scan(datafile, n=NPool*1, quiet=TRUE),
##             NPool, 1, byrow=TRUE)
## d2 <- density(z[,1], from=0, to=3000)

## SimResultB <- spline(d2, n=101, xmin=0, xmax=3000)
## SimResultB$y <- pmax(0, SimResultB$y)
## SimResultB$y <- SimResultB$y/(sum(SimResultB$y)*30)
## ##lines(SimResultB$x, SimResultB$y, type='l', col='red')
## y.lim = c(0, 3e-3)
## y.at <- pretty(y.lim, n=4)
## y.labels <- formatC(1000*y.at, format = "g")
## p2 <- xyplot(d1$y~d1$x,
##              xlim=c(0,3000),
##              ylim=c(0, 0.003),
##              lwd=2,
##              type='l',
##              col='black',
##              xlab=list('s (cell size, fl)', cex=1),
##              ##ylab=list('density frequency', cex = 1),
##              ylab='',
##              scales=list(cex=1, y=list(at = y.at,
##                                       labels = y.labels), tck=0.5),
##              panel=function(...){
##                  panel.lines(SimResultB$y~SimResultB$x, lwd=2,
##                              type='l', col='red')
##                  panel.lines(ExpResultB$y~ExpResultB$x, lwd=2,
##                              type='l', col='blue')                 
##                   panel.axis(side = c("top"),
##                              at = c(100),
##                              labels = c(expression(x10^{-3})),
##                              ticks=FALSE,
##                              outside=TRUE,
##                              rot=c(0,90),
##                              text.cex=1)
##              }
##              )
## print(p2, position=c(0.5, 0, 1, 1))
## ## lines(ExpResultB$x, ExpResultB$y, type='l', col='blue')
## ## ##ErrorB[i] <- norm(as.matrix(ExpResultB$y-SimResultB$y), type='O')
## ## ErrorB <- norm(as.matrix(ExpResultB$y-SimResultB$y), type='O')
## ## ## tt1<-cbind(ExpResultB$x, ExpResultB$y)
## ## ## tt2<-cbind(SimResultB$x, SimResultB$y)
## ## ## ErrorB[i] <- kl.dist(tt1,tt2)$D1

## ErrorT <- ErrorA + ErrorB
dev.off()
#}
