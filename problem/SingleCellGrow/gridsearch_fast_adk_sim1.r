library('adk')
dyn.load('simulator1.so')
source('simulator1n.r')

# Simulated data
.T <- 8
.r1 <- 1
.r2 <- 0.6
.d1 <- 0.3
.d2 <- 0.3
.v <- 0.3
.w <- 0.4
.N <- 10000

mcell <- simulator1n(.T, .r1, .r2, .d1, .d2, .v, .w, .N)

## Nsample = 10000
## T = 8
## r2 = 0.6
## d2 = 0.3
## d <- double(1000000)
## pvalue <- double(1000000)
## x <- seq(Nsample)
## j <- 1
## for ( r1 in seq(0.5, 1.5, by=0.1) ){
##     for ( d1 in seq(0, 1, by=0.1) ){
##         for ( v in seq(0, 1, by=0.1) ){
##           for ( w in seq(0, 1, by=0.1) ){
##             x <- simulator1n(T, r1, r2, d1, d2, v, w, Nsample)
##             #dis <- adk.test(mcell, x)
##             #d[j] <- dis$adk[1,2]
##             #pvalue[j] <- dis$adk[2,2]

##             dis <- ks.test(mcell, x)
##             d[j] <- dis$statistic
##             pvalue[j] <- dis$p.value
##             cat(c(d[j], r1, r2, d1, d2, v, w),'\n')
##             j <- j+1
##           }
##       }
##   }
## }


## j <- 1
## n <- 0
## m <- 0
## dvsr1 <- seq(1000000)
## dvsd1 <- seq(1000000)
## dvsv <- seq(1000000)
## dvsw <- seq(1000000)
## pl <- seq(1000000)
## for ( r1 in seq(0.5, 1.5, by=0.1) ){
##     for ( d1 in seq(0, 1, by=0.1) ){
##         for ( v in seq(0, 1, by=0.1) ){
##           for ( w in seq(0, 1, by=0.1) ){
## 	    if ( (pvalue[j] > 0.8) && (pvalue[j] < 1) ){
##                 m <- m+1
##                 dvsr1[m] <- r1
##                 dvsd1[m] <- d1
##                 dvsv[m] <- v
##                 dvsw[m] <- w
##                 pl[m] <- pvalue[j]
##                 cat(m, j, '\n')
##               }
##             j <- j+1
##           }
##         }
##   }
## }

## par(mfrow=c(2,2))
## plot(dvsr1[1:m], pl[1:m])
## plot(dvsd1[1:m], pl[1:m])
## plot(dvsv[1:m], pl[1:m])
## plot(dvsw[1:m], pl[1:m])

## -----------------------------------
## Plot the fit
## -----------------------------------

# Simulated data
.T <- 8
.r1 <- 1
.r2 <- 0.6
.d1 <- 0.3
.d2 <- 0.3
.v <- 0.3
.w <- 0.4
.N <- 10000

mcell <- simulator1n(.T, .r1, .r2, .d1, .d2, .v, .w, .N)

Nsample = 10000
y <- simulator1n(.T, .r1, .r2, .d1, .d2, .v, .w, Nsample)

T = 8
r1 = 1.5
r2 = 0.6
d1 = 0.7
d2 = 0.3
v = 0.8
w = 0.9
x <- simulator1n(T, r1, r2, d1, d2, v, w, Nsample)

Fn <- ecdf(x)
Fe <- ecdf(mcell)
Fy <- ecdf(y)
plot(Fn, xlim=c(0,300),
     main = "ECDF: simulation vs data",
     ylab = "Cumulative probability",
     xlab = "8-day Cell number ")
lines(Fe, col='red')
lines(Fy, col='green')

#dis <- adk.test(mcell, x)
#cat(dis$adk[1,2])

## dis <- ks.test(mcell, x)
## cat(dis$statistic, dis$p.value, '\n')

## disy <- ks.test(mcell, y)
## cat(disy$statistic, disy$p.value, '\n')

dis1 <- adk.test(mcell, x)
dis2 <- adk.test(mcell, y)

mtext(paste("@", date()), side=1, line=4, adj=1.04, cex=.66)

#dev.copy(jpeg,'search6.jpg')
#dev.off()
