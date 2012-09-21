library('adk')
dyn.load('simulator1.so')
source('simulator1n.r')


# Simulated data using r2 = 0.6, d2 = 0.3, 50 samples.
mcell <- as.matrix(read.csv('data/simcella.csv'))

## Nsample = 2000
## T = 8
## d <- matrix(0, 101, 101)
## pvalue <- matrix(0, 101, 101)
## x <- seq(Nsample)
## i <- 1
## for ( r2 in seq(0, 1, by=0.01) ){
##     j <- 1
##     for ( d2 in seq(0, 1, by=0.01) ){
##             x <- simulator1n(T, 0, r2, 0, d2, 0, 0, Nsample)
##             dis <- adk.test(mcell, x)
##             d[i, j] <- dis$adk[1,2]
##             pvalue[i, j] <- dis$adk[2,2]
##             cat(c(d[i, j], r2, d2, i, j),'\n')
##             j <- j + 1
##     }
##     i <- i + 1
## }

## -----------------------------------
## Plot the contour
## -----------------------------------

filled.contour(x = seq(0, 1, length.out=101),
		 y = seq(0, 1, length.out=101),
		 d,
		 color=terrain.colors,
		 plot.title = title(main = "KS-distance between ECDFs [Good Cells]",
		 xlab = "proliferation rate",
		 ylab = "death rate"),
		 asp = 1,
		 plot.axes={ axis(1); axis(2); points(0.6,0.3,pch=17) },
		 level=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.65, 0.7))

mtext(paste("@", date()), side=1, line=4, adj=1.04, cex=.66)

text(0.4, 1.3, "Grid search for best fit using one-species model.")
text(0.4, 1.25, "The best fit is located at")
text(0.4, 1.2, "r1 = 0.84, d1 = 0.44 (triangle)")

dev.copy(pdf,'search5.pdf')
dev.off()


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
