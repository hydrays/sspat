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

filled.contour(x = seq(0, 1.4, length.out=71), 
		 y = seq(0, 1.4, length.out=71), 
		 d[1:71,1:71], 
		 color=terrain.colors,
		 plot.title = title(main = "KS-distance between ECDFs [Good Cells]",
		 xlab = "proliferation rate", 
		 ylab = "death rate"),
		 asp = 1, 
		 plot.axes={ axis(1); axis(2); points(0.84,0.44,pch=17) },
		 level=c(0, 0.12, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1))

mtext(paste("@", date()), side=1, line=4, adj=1.04, cex=.66)

text(0.4, 1.3, "Grid search for best fit using one-species model.")
text(0.4, 1.25, "The best fit is located at")
text(0.4, 1.2, "r1 = 0.84, d1 = 0.44 (triangle)")

if( require(png) ) {
  img <- readPNG("bear2.png")
  my.symbols( 1.04, -0.36, ms.image, MoreArgs=list(img=img),
             inches=0.2, symb.plots=TRUE, add=TRUE)
}

dev.copy(pdf,'search5.pdf')
dev.off()
