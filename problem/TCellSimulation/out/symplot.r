require(lattice)
require(rasterVis)
require(latticeExtra)
library(viridisLite)
coul <- viridis(100)
require("RColorBrewer")
gray.cols <- colorRampPalette(c('black', 'white'))(n=100)

## ##color Scheme I
## chemokine.cols <- colorRampPalette(c('white', '#F8766D'))(n=100)
## fb.col <- "red"
## tcell.col <- "#1F78B4"

## Color Scheme II
##chemokine.cols <- colorRampPalette(c('white', '#1F78B4'))(n=100)
##fb.col <- "#0000FF"
##tcell.col <- "#F8766D"

## ## Color Scheme III
## chemokine.cols <- colorRampPalette(c('white', '#FFFF00'))(n=100)
## fb.col <- "#1F78B4"
## tcell.col <- "#F8766D"

## ## Color Scheme IV
## chemokine.cols <- colorRampPalette(c('white', '#7CAE00'))(n=100)
## fb.col <- "#1F78B4"
## tcell.col <- "#F8766D"

## ## Color Scheme V
## chemokine.cols <- colorRampPalette(c('white', '#F7931E'))(n=100)
## fb.col <- "#1F78B4"
## tcell.col <- "#F8766D"

## Color Scheme VI
chemokine.cols <- colorRampPalette(c('white', '#FFFF00'))(n=100)
fb.col <- "#DEF1FC"
tcell.col <- "#F8766D"

highsig.col <- "#AFDCE8"
weaksig.col <- "#DEF1FC"
dark.col <- "#9FA0A0"
white.col <- "#FFFFFF"

fbhigh.col <- "#1F78B6"
fbweak.col <- "#95DEEF"
mycol = c(fbweak.col, fbhigh.col)
##L <- 2048
L <- 1024
##L  <- 512
##L  <- 256
#L  <- 128
colors <- c("#A7A7A7",
            "dodgerblue",
            "firebrick",
            "forestgreen",
            "gold")
myAt <- seq(5) - 1.5

pic <- matrix(unlist(read.csv('../matrix.txt', sep='', header=F)), nrow=L)

for( i in seq(0, 2000, by=10) )
{
    cat(i, '\n')
    padded_i <- sprintf("%05d", i)
    png(paste("config_", padded_i, ".png", sep=''), height=1*600, width=1*900)
    zz1 = raster(pic, xmn=0, xmx=L, ymn=0, ymx=L)
    zz1 = zz1 - 1
    p1l1  <- levelplot(zz1, col.regions = c(highsig.col, weaksig.col),  margin=FALSE,
                     colorkey=FALSE)

    A <- matrix(unlist(read.csv(paste('c', padded_i, '.dat', sep=''), header=F)), nrow=L)
    A[A==2] <- 0
    z1 = raster(A, xmn=0, xmx=L, ymn=0, ymx=L)
    coords <- xyFromCell(z1, which(z1[]==1))
    p1l2 <- xyplot(coords[,2]~coords[,1], lwd=1, pch=8, cex=0.001, col="#C9CACA")

    coords <- xyFromCell(z1, which(z1[]==3))
    p1l3 <- xyplot(coords[,2]~coords[,1], lwd=1, cex=0.05, col=tcell.col)    
    p1 <- p1l1 + as.layer(p1l2) + as.layer(p1l3)

    z1[z1==2] <- 0
    z1[z1==3] <- 0
    p2l1 <- levelplot(z1, col.regions = c(white.col), margin=FALSE, colorkey=FALSE)
    coords <- xyFromCell(z1, which(z1[]==1))
    p2l2 <- xyplot(coords[,2]~coords[,1], lwd=1, cex=0.5, col=dark.col)    
    p2 <- p2l1 + as.layer(p2l2)

    ##B <- matrix(unlist(read.csv(paste('phi', padded_i, '.dat', sep=''), header=F)), nrow=L)
    B <- matrix(unlist(read.csv(paste('phi', padded_i, '.dat', sep=''), header=F)), nrow=L)
    ##coords <- xyFromCell(z, which(z[]==3))
    z2 = raster(B, xmn=0, xmx=L, ymn=0, ymx=L)
    z2[1] = 1
    zlog2 <- log(z2)    
    p3l1 <- levelplot(zlog2, col.regions = chemokine.cols, margin=FALSE)
    C <- matrix(unlist(read.csv(paste('lambda', padded_i, '.dat', sep=''), header=F)), nrow=L)
    z3 = raster(C, xmn=0, xmx=L, ymn=0, ymx=L)    
    coords <- xyFromCell(z3, which(z3[]>1e-2))
    ##p3l2 <- xyplot(coords[,2]~coords[,1], cex=0.001, col=c(fbweak.col, fbhigh.col))
    p3l2 <- xyplot(coords[,2]~coords[,1], cex=0.001, col=mycol[z3[z3[]>1e-2]])
    p3 <- p3l1 + as.layer(p3l2)
    
    print(p1, position=c(0, 0.0, .5, 1), more=TRUE)
    print(p2, position=c(0.5, 0.5, 1, 1), more=TRUE)
    print(p3, position=c(0.5, 0, 1, 0.5))
    dev.off()
}

