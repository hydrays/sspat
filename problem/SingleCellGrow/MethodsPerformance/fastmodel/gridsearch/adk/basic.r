library('adk')
dyn.load('simulator1.so')
source('simulator1n.r')

## # Simulated data 
## .T <- 8
## .r1 <- 1
## .r2 <- 0.6
## .d1 <- 0.4
## .d2 <- 0.3
## .v <- 0.3
## .w <- 0.7
## .N <- 100

## #mcell <- as.matrix(read.csv('data/simcella.csv'))
## mcell <- simulator1n(.T, .r1, .r2, .d1, .d2, .v, .w, .N)

## Nsample = 2000
## T = 8
## r2 = 0.6
## d2 = 0.3
## d <- double(1000000)
## pvalue <- double(1000000)
## x <- seq(Nsample)
## j <- 1
## for ( r1 in seq(0, 1.5, by=0.05) ){
##     for ( d1 in seq(0, 1, by=0.05) ){
##         for ( v in seq(0, 1, by=0.05) ){
##           for ( w in seq(0, 1, by=0.05) ){
##             x <- simulator1n(T, r1, r2, d1, d2, v, w, Nsample)
##             dis <- adk.test(mcell, x)
##             d[j] <- dis$adk[1,2]
##             pvalue[j] <- dis$adk[2,2]
##             cat(c(d[j], r1, r2, d1, d2, v, w),'\n')
##             j <- j+1
##           }
##       }
##   }
## }


j <- 1
n <- 0
m <- 0
dvsr1 <- seq(1000000)
dvsd1 <- seq(1000000)
dvsv <- seq(1000000)
dvsw <- seq(1000000)
pl <- seq(1000000)
for ( r1 in seq(0, 2, by=0.05) ){
    for ( d1 in seq(0, 1, by=0.05) ){
        for ( v in seq(0, 1, by=0.05) ){
          for ( w in seq(0, 1, by=0.05) ){
	    if ( (pvalue[j] > 0.69) && (pvalue[j] < 1) ){
                m <- m+1
                dvsr1[m] <- r1
                dvsd1[m] <- d1
                dvsv[m] <- v
                dvsw[m] <- w
                pl[m] <- pvalue[j]
                cat(m, j, '\n')
              }
            j <- j+1
          }
        }
  }
}

## par(mfrow=c(2,2))
## plot(dvsr1[1:m], pl[1:m])
## plot(dvsd1[1:m], pl[1:m])
## plot(dvsv[1:m], pl[1:m])
## plot(dvsw[1:m], pl[1:m])

## ## -----------------------------------
## ## Plot the fit
## ## -----------------------------------

## par(mfrow=c(1,1))
## Nsample = 10000
## T = 8
## r1 = 1.0
## r2 = 0.6
## d1 = 0.35
## d2 = 0.3
## v = 0.35
## w = 0.75
## x <- seq(Nsample)
## for (i in seq(Nsample)){
##   y <- Csimulator1(T, r1, r2, d1, d2, v, w)
##   x[i] <- sum(y)
## }

## data7 <- read.csv('data/8dayfast.csv')
## mcell <- as.numeric(data7[8,3:ncol(data7[1,])])
## mcell[is.na(mcell)] <- 0

## Fn <- ecdf(x)
## Fe <- ecdf(mcell)
## plot(Fn, xlim=c(0,300),
##      main = "ECDF: simulation vs data",
##      ylab = "Cumulative probability",
##      xlab = "8-day Cell number ")
## lines(Fe, col='red')

## text(200, 0.6, "Goodness of fit")
## text(200, 0.55, "Black: simulated using parameters from grid search")
## text(200, 0.5, "Red: experiment data of the mixed population")
## text(200, 0.45, "Result (r1=1, r2=0.6, d1=0.4, d2=0.3, v=0, w=0.3):")
## text(200, 0.4, "KS-distance: 0.08")
## text(200, 0.35, "p-value: 0.81")

## dis <- adk.test(mcell, x)
## cat(dis$adk[1,2])

## library('png')
## if( require(png) ) {
##   img <- readPNG("bear.png")
##   my.symbols( x = unit(0.8, "npc"), 0, ms.image, MoreArgs=list(img=img),
##              inches=0.2, symb.plots=TRUE, add=TRUE)
## }
## mtext(paste("@", date()), side=1, line=4, adj=1.04, cex=.66)

## dev.copy(jpeg,'search6.jpg')
## dev.off()
