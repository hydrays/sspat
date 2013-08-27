#######################################
## Asynchronized cell size distribution
#######################################

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

nn <- 50

eErrorA <- seq(nn)
tErrorA <- seq(nn)
sErrorA <- seq(nn)

eErrorT <- seq(nn)
tErrorT <- seq(nn)
sErrorT <- seq(nn)

for (i in seq(nn)) {
datafile <- sprintf("%s%05d%s", "trout_event/m", i, ".dat")
z <- matrix(scan(datafile, n=NPool*6, quiet=TRUE),
            NPool, 6, byrow=TRUE)
ed1 <- density(z[,1], from=0, to=3000)
eSimResultA <- spline(ed1, n=101, xmin=0, xmax=3000)
eSimResultA$y <- pmax(0, eSimResultA$y)
eSimResultA$y <- eSimResultA$y/(sum(eSimResultA$y)*30)
eErrorA[i] <- norm(as.matrix(ExpResultA$y-eSimResultA$y), type='O')

datafile <- sprintf("%s%05d%s", "trout_timeonly/m", i, ".dat")
z <- matrix(scan(datafile, n=NPool*6, quiet=TRUE),
            NPool, 6, byrow=TRUE)
td1 <- density(z[,1], from=0, to=3000)
tSimResultA <- spline(td1, n=101, xmin=0, xmax=3000)
tSimResultA$y <- pmax(0, tSimResultA$y)
tSimResultA$y <- tSimResultA$y/(sum(tSimResultA$y)*30)
tErrorA[i] <- norm(as.matrix(ExpResultA$y-tSimResultA$y), type='O')

datafile <- sprintf("%s%05d%s", "trout_sizetime/m", i, ".dat")
z <- matrix(scan(datafile, n=NPool*6, quiet=TRUE),
            NPool, 6, byrow=TRUE)
sd1 <- density(z[,1], from=0, to=3000)
sSimResultA <- spline(sd1, n=101, xmin=0, xmax=3000)
sSimResultA$y <- pmax(0, sSimResultA$y)
sSimResultA$y <- sSimResultA$y/(sum(sSimResultA$y)*30)
sErrorA[i] <- norm(as.matrix(ExpResultA$y-sSimResultA$y), type='O')
}

eErrorT <- eErrorA
tErrorT <- tErrorA
sErrorT <- sErrorA

## Newborn
for (i in seq(nn)) {
datafile <- sprintf("%s%05d%s", "trout_event/n", i, ".dat")
z <- matrix(scan(datafile, n=NPool*1, quiet=TRUE),
            NPool, 1, byrow=TRUE)
ed1 <- density(z[,1], from=0, to=3000)
eSimResultA <- spline(ed1, n=101, xmin=0, xmax=3000)
eSimResultA$y <- pmax(0, eSimResultA$y)
eSimResultA$y <- eSimResultA$y/(sum(eSimResultA$y)*30)
eErrorA[i] <- norm(as.matrix(ExpResultB$y-eSimResultA$y), type='O')

datafile <- sprintf("%s%05d%s", "trout_timeonly/n", i, ".dat")
z <- matrix(scan(datafile, n=NPool*1, quiet=TRUE),
            NPool, 1, byrow=TRUE)
td1 <- density(z[,1], from=0, to=3000)
tSimResultA <- spline(td1, n=101, xmin=0, xmax=3000)
tSimResultA$y <- pmax(0, tSimResultA$y)
tSimResultA$y <- tSimResultA$y/(sum(tSimResultA$y)*30)
tErrorA[i] <- norm(as.matrix(ExpResultB$y-tSimResultA$y), type='O')

datafile <- sprintf("%s%05d%s", "trout_sizetime/n", i, ".dat")
z <- matrix(scan(datafile, n=NPool*1, quiet=TRUE),
            NPool, 1, byrow=TRUE)
sd1 <- density(z[,1], from=0, to=3000)
sSimResultA <- spline(sd1, n=101, xmin=0, xmax=3000)
sSimResultA$y <- pmax(0, sSimResultA$y)
sSimResultA$y <- sSimResultA$y/(sum(sSimResultA$y)*30)
sErrorA[i] <- norm(as.matrix(ExpResultB$y-sSimResultA$y), type='O')
}

eErrorT <- eErrorA + eErrorT
tErrorT <- tErrorA + tErrorT
sErrorT <- sErrorA + sErrorT

eErrorB <- eErrorT - eErrorA
tErrorB <- tErrorT - tErrorA
sErrorB <- sErrorT - sErrorA

require('lattice')
pdf("fig_conv.pdf", width=7, height=5)
xyplot(eErrorT~10*seq(1:50),
       #xlim=c(0,500),
       #ylim=c(0, 0.0015),
       lwd=2,
       type='b',
       col='blue',
       xlab=list('time (hour)', cex=2),
       ylab=list('L1 error', cex = 2),
       scale=list(cex=2),
       key = list(x = 0.45, y=0.9,
           border=FALSE,
           lines=list(
               pch=c(2,3,1),
               col=c("brown","red","blue"),
               type=c("b","b","b"),
               lty=c(1,1,1),
               lwd=c(2,2,2)),
               cex = 1,
               text = list(lab = c("time only","size and time",
                               "integration"), cex=1.5),
               columns = 1,
               title = NULL
               ),       
       panel=function(...){
           panel.lines(...)
           panel.lines(tErrorT~10*seq(1:50), lwd=2,
                       type='b', pch=2, col='brown')
           panel.lines(sErrorT~10*seq(1:50), lwd=2,
                             type='b', pch=3, col='red')
             }
             )
dev.off()

