## data7 <- read.csv('8day.csv')
## mcell <- as.numeric(data7[8,3:(ncol(data7[1,]))])
## mcell[is.na(mcell)] <- 0

## Nsample = 1000
## T = 8
## x <- seq(Nsample)
## d <- matrix(ncol=101, nrow=101)
## pvalue <- matrix(ncol=101, nrow=101)
## i <- 1
## for ( r1 in seq(0, 2, by=0.02) ){
##     j <- 1
##       for ( d1 in seq(0, 2, by=0.02) ){
##               for (k in seq(Nsample)){
##                 y <- simulator1(T, r1, 0, d1, 0, 0, 1)
##                 x[k] <- sum(y)
## 		}
##               dis <- ks.test(mcell, x)
##               pvalue[i,j] <- dis$p.value
##               d[i,j] <- dis$statistic
##               cat(c(d[i,j], pvalue[i,j], r1, d1, i, j),'\n')
##             j <- j+1
##         }
## 	i <- i+1
## }


## -----------------------------------
## Plot the contour
## -----------------------------------

## filled.contour(x = seq(0, 1.4, length.out=71), 
## 		 y = seq(0, 1.4, length.out=71), 
## 		 d[1:71,1:71], 
## 		 color=terrain.colors,
## 		 plot.title = title(main = "KS-distance between ECDFs [Good Cells]",
## 		 xlab = "proliferation rate", 
## 		 ylab = "death rate"),
## 		 asp = 1, 
## 		 plot.axes={ axis(1); axis(2); points(0.84,0.44,pch=17) },
## 		 level=c(0, 0.12, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1))

## mtext(paste("@", date()), side=1, line=4, adj=1.04, cex=.66)

## text(0.4, 1.3, "Grid search for best fit using one-species model.")
## text(0.4, 1.25, "The best fit is located at")
## text(0.4, 1.2, "r1 = 0.84, d1 = 0.44 (triangle)")

## if( require(png) ) {
##   img <- readPNG("bear.png")
##   my.symbols( 1.04, -0.36, ms.image, MoreArgs=list(img=img),
##              inches=0.2, symb.plots=TRUE, add=TRUE)
## }

## dev.copy(pdf,'search5.pdf')
## dev.off()

## -----------------------------------
## Plot the fit
## -----------------------------------
Nsample = 10000
T = 8
r1 = 0.84
r2 = 0
d1 = 0.44
d2 = 0
p = 0
w = 1 # means all good cell
x <- seq(Nsample)
for (i in seq(Nsample)){
  y <- simulator1(T, r1, r2, d1, d2, p, w)
  x[i] <- sum(y)
}

data7 <- read.csv('data/8dayfast.csv')
mcell <- as.numeric(data7[8,3:ncol(data7[1,])])
mcell[is.na(mcell)] <- 0

Fn <- ecdf(x)
Fe <- ecdf(mcell)
plot(Fn, xlim=c(0,300),
     main = "ECDF: simulation vs data",
     ylab = "Cumulative probability",
     xlab = "8-day Cell number ")
lines(Fe, col='red')

text(200, 0.6, "Goodness of fit")
text(200, 0.55, "Black: simulated using parameters from grid search")
text(200, 0.5, "Red: experiment data of the mixed population")
text(200, 0.45, "Result:")
text(200, 0.4, "KS-distance: 0.11")
text(200, 0.35, "p-value: 0.55")

dis <- ks.test(mcell, x)
cat(dis$p.value)

if( require(png) ) {
  img <- readPNG("bear.png")
  my.symbols( 230, -0.22, ms.image, MoreArgs=list(img=img),
             inches=0.2, symb.plots=TRUE, add=TRUE)
}
mtext(paste("@", date()), side=1, line=4, adj=1.04, cex=.66)

#dev.copy(pdf,'search5.pdf')
dev.copy(jpeg,'search5_2.jpg')
dev.off()
