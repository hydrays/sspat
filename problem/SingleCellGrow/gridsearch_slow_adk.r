library('adk')
dyn.load('simulator1.so')
source('simulator1n.r')

mcell <- as.matrix(read.csv('data/day8slow.csv'))

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



## ## -----------------------------------
## ## Plot the contour
## ## -----------------------------------

## filled.contour(x = seq(0, 1, length.out=101),
## 		 y = seq(0, 1, length.out=101),
## 		 d,
## 		 color=terrain.colors,
## 		 plot.title = title(main = "KS-distance between ECDFs [Good Cells]",
## 		 xlab = "proliferation rate",
## 		 ylab = "death rate"),
## 		 asp = 1,
## 		 plot.axes={ axis(1); axis(2); points(0.65,0.33,pch=1) },
## 		 level=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.66, 0.7))

## mtext(paste("@", date()), side=1, line=4, adj=1.04, cex=.66)

## text(0.4, 1.3, "Grid search for best fit using one-species model.")
## text(0.4, 1.25, "The best fit is located at")
## text(0.4, 1.2, "r1 = 0.84, d1 = 0.44 (triangle)")

## dev.copy(pdf,'search5.pdf')
## dev.off()


## -----------------------------------
## Plot the fit
## -----------------------------------

Nsample = 10000
T = 8
r1 = 0
r2 = 0.6
d1 = 0
d2 = 0.3
v = 0
w = 0
x <- simulator1n(T, r1, r2, d1, d2, v, w, Nsample)

Fn <- ecdf(x)
Fe <- ecdf(mcell)
plot(Fn, xlim=c(0,100),
     main = "ECDF: simulation vs data",
     ylab = "Cumulative probability",
     xlab = "8-day Cell number ")
lines(Fe, col='red')

text(200, 0.6, "Goodness of fit")
text(200, 0.55, "Black: simulated using parameters from grid search")
text(200, 0.5, "Red: experiment data of the mixed population")
text(200, 0.45, "Result (r1=0, r2=0.65, d1=0, d2=0.33, v=0, w=0):")
text(200, 0.4, "KS-distance: 0.08")
text(200, 0.35, "p-value: 0.81")

dis <- adk.test(mcell, x)
cat(dis$adk[1,2])

dev.copy(jpeg,'search6.jpg')
dev.off()
