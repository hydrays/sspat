require(lattice)
require(rasterVis)
require(latticeExtra)
library(viridisLite)
coul <- viridis(100)
gray.cols <- colorRampPalette(c('black', 'white'))(n=100)

##color Scheme I
chemokine.cols <- colorRampPalette(c('white', '#F8766D'))(n=100)
fb.col <- "red"
tcell.col <- "#1F78B4"

L  <- 512
##L  <- 256
#L  <- 128
colors <- c("#A7A7A7",
            "dodgerblue",
            "firebrick",
            "forestgreen",
            "gold")
myAt <- seq(5) - 1.5

pic <- matrix(unlist(read.csv('f00000.dat', header=F)), nrow=L)

for( i in seq(0, 2000, by=10) )
{
    cat(i, '\n')
    padded_i <- sprintf("%05d", i)
    png(paste("config_", padded_i, ".png", sep=''), height=600, width=1200)
    A <- matrix(unlist(read.csv(paste('c', padded_i, '.dat', sep=''), header=F)), nrow=L)
    A[A==2] <- 0
    z1 = raster(A, xmn=0, xmx=L, ymn=0, ymx=L)
    zz1 = raster(pic, xmn=0, xmx=L, ymn=0, ymx=L)
    p1  <- levelplot(zz1, col.regions = gray.cols, margin=FALSE,
                     colorkey=TRUE)
    coords <- xyFromCell(z1, which(z1[]==3))
    pp3 <- xyplot(coords[,2]~coords[,1], lwd=1, cex=0.05, col=tcell.col)    
    p1 <- p1 + as.layer(pp3)
    ##p1 <- as.layer(pp1) + p1

    ##B <- matrix(unlist(read.csv(paste('phi', padded_i, '.dat', sep=''), header=F)), nrow=L)
    B <- matrix(unlist(read.csv(paste('phi', padded_i, '.dat', sep=''), header=F)), nrow=L)
    ##coords <- xyFromCell(z, which(z[]==3))
    z2 = raster(B, xmn=0, xmx=L, ymn=0, ymx=L)
    p2 <- levelplot(z2, col.regions = coul, margin=FALSE,
                    panel=function(..., at, contour, region) {
                        panel.levelplot(..., contour = FALSE)
                        panel.contourplot(..., at=c(1), contour = TRUE, col="white", lty=2, region = FALSE)
                    }
                    )
    ## C <- matrix(unlist(read.csv(paste('lambda', padded_i, '.dat', sep=''), header=F)), nrow=L)
    ## z3 = raster(C, xmn=0, xmx=256, ymn=0, ymx=256)    
    ## coords <- xyFromCell(z3, which(z3[]>1e-2))
    ## p3 <- xyplot(coords[,2]~coords[,1], lwd=0.5, pch=1, cex=0.5, col="red")
    ## p4 <- p2 + as.layer(p3)

    ##C <- matrix(unlist(read.csv(paste('phi', padded_i, '.dat', sep=''), header=F)), nrow=L)
    C <- matrix(unlist(read.csv(paste('lambda', padded_i, '.dat', sep=''), header=F)), nrow=L)
    z3 = raster(C, xmn=0, xmx=L, ymn=0, ymx=L)
    p3 <- levelplot(z3, col.regions = coul, margin=FALSE)
    ## z3 = raster(C, xmn=0, xmx=256, ymn=0, ymx=256)    
    ## coords <- xyFromCell(z3, which(z3[]>1e-2))
    ## p3 <- xyplot(coords[,2]~coords[,1], lwd=0.5, pch=1, cex=0.5, col="red")
    ## p4 <- p2 + as.layer(p3)
    
    print(p1, position=c(0, 0.0, .5, 1), more=TRUE)
    print(p2, position=c(0.5, 0.5, 1, 1), more=TRUE)
    print(p3, position=c(0.5, 0, 1, 0.5))
    dev.off()
}

