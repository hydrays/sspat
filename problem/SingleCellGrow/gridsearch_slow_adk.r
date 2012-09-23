library('adk')
dyn.load('simulator1.so')
source('simulator1n.r')

mcell <- as.matrix(read.csv('data/day8slow.csv'))

## Nsample = 10000
## T = 8
## d <- matrix(0, 101, 101)
## pvalue <- matrix(0, 101, 101)
## x <- seq(Nsample)
## i <- 1
## for ( r2 in seq(0.4, 0.8, by=0.004) ){
##     j <- 1
##     for ( d2 in seq(0.2, 0.6, by=0.004) ){
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

## res <- which(d == max(d), arr.ind=T)
## r2max <- 0.4 + 0.004*res[1]
## d2max <- 0.2 + 0.004*res[2]

## par(bg = 'white')
## filled.contour(x = seq(0, 1, length.out=101),
##                y = seq(0, 1, length.out=101),
##                d,
##                color=terrain.colors,
##                plot.title = title(main = "Contour-plot of p-value under AD test (2D)",
## 		 xlab = "proliferation rate",
## 		 ylab = "death rate"),
##                asp = 1,
##                ## plot.axes={ axis(1); axis(2); points(r2max,d2max,pch=17) },
##                level=c(0, 0.0001, 0.01, 0.1, 0.2, 0.3, 0.5, 0.6, 0.65, 0.7))

## mtext(paste("@", date()), side=1, line=4, adj=1.04, cex=.66)

## text(0, 0.9, adj = c(0,0), "The method of fitting parameters is based on")
## text(0, 0.85, adj = c(0,0), "Minimize the AD-distance (or maximize its p-Value).")

## text(0, 0.75, adj = c(0,0), "Sample A: experimental data (day 8 cell number).")
## text(0, 0.7, adj = c(0,0), "Sample B: simulated data from the population model.")

## text(0, 0.6, adj = c(0,0), "We do a grid-search on the whole parameter space.")
## text(0, 0.55, adj = c(0,0), "Each parameter set gives a sample B.")
## text(0, 0.5, adj = c(0,0), "The parameter set that minimize the distance")
## text(0, 0.45, adj = c(0,0), "(or maximize the p-value) between sample A-B")
## text(0, 0.4, adj = c(0,0), "is the winner.")

## arrows(x0=0.2, y0=0.4, x1=0.4, y1=0.25)

## dev.copy(pdf,'I_method.pdf')
## dev.off()


## -----------------------------------
## Plot the fit
## -----------------------------------

Nsample = 10000
T = 8
r1 = 0
r2 = 0.652
d1 = 0
d2 = 0.34
v = 0
w = 0
x <- simulator1n(T, r1, r2, d1, d2, v, w, Nsample)

#y <- simulator1n(T, r1, 0.6, d1, 0.3, v, w, Nsample)

Fn <- ecdf(x)
Fe <- ecdf(mcell)
plot(Fn, xlim=c(0,100),
     main = "ECDF: simulation vs data",
     ylab = "Cumulative probability",
     xlab = "8-day Cell number ")
lines(Fe, col='red')

#Fy <- ecdf(y)
#lines(Fy, col='blue')

text(20, 0.5, adj = c(0, 0), "Model parameter estimated by minimize AD-distance.")
text(20, 0.45, adj = c(0, 0), "Red: 8-day slow cell population")
text(20, 0.4, adj = c(0, 0), "Black: simulated data using parameters from grid search")
text(20, 0.35, adj = c(0, 0), "Two model parameters: ( l2=0.652, d2=0.33 )")
text(20, 0.3, adj = c(0, 0), "max p-value: 0.66")

dis <- adk.test(mcell, x)
cat(dis$adk[1,2])

mtext(paste("@", date()), side=1, line=4, adj=1.04, cex=.66)

dev.copy(jpeg,'I_slowcell.jpg')
dev.off()
