require(lattice)
require(rasterVis)
require(latticeExtra)
library(viridisLite)
coul <- viridis(100)
require("RColorBrewer")
mycols <- colorRampPalette(c('white', '#1F78B4'))(n=100)

L  <- 512
##L  <- 256
#L  <- 128
colors <- c("#A7A7A7",
            "dodgerblue",
            "firebrick",
            "forestgreen",
            "gold")
myAt <- seq(5) - 1.5

for( i in seq(0, 2000, by=1) )
{
    cat(i, '\n')
    padded_i <- sprintf("%05d", i)
    png(paste("config_", padded_i, ".png", sep=''), height=600, width=1200)
    A <- matrix(unlist(read.csv(paste('c', padded_i, '.dat', sep=''), header=F)), nrow=L)
    A[A==2] <- sample(c(0,2), length(A[A==2]), replace = TRUE, prob=c(0.9, 0.1))
    z1 = raster(A, xmn=0, xmx=L, ymn=0, ymx=L)
    znull = z1
    znull[] = 0
    p1  <- levelplot(znull, at=myAt, col.regions = c('white', '#C9CACA', '#F8766D', '#F8766D'), margin=FALSE,
                     colorkey=FALSE)
    coords <- xyFromCell(z1, which(z1[]==1))
    pp1 <- xyplot(coords[,2]~coords[,1], lwd=1, pch=8, cex=0.5, col="#C9CACA")
    coords <- xyFromCell(z1, which(A[]==2))
    pp2 <- xyplot(coords[,2]~coords[,1], lwd=1, cex=0.25, col="#F8766D")
    coords <- xyFromCell(z1, which(z1[]==3))
    pp3 <- xyplot(coords[,2]~coords[,1], lwd=1, cex=0.25, col="#F8766D")    
    p1 <- p1 + as.layer(pp1) + as.layer(pp2) + as.layer(pp3)
    ##p1 <- as.layer(pp1) + p1
    
    ##B <- matrix(unlist(read.csv(paste('phi', padded_i, '.dat', sep=''), header=F)), nrow=L)
    B <- matrix(unlist(read.csv(paste('phi', padded_i, '.dat', sep=''), header=F)), nrow=L)
    ##coords <- xyFromCell(z, which(z[]==3))
    z2 = raster(B, xmn=0, xmx=L, ymn=0, ymx=L)
    zlog2 <- log(z2)
    p2 <- levelplot(z2, col.regions = mycols, margin=FALSE,
                    panel=function(..., at, contour, region) {
                        panel.levelplot(..., at, contour = FALSE)
                        panel.contourplot(..., at=c(0.5), contour = TRUE, col="green", lty=2, lwd=2, region = FALSE)
                        ##panel.dotplot(y~x)
                    }
                    )
    C <- matrix(unlist(read.csv(paste('lambda', padded_i, '.dat', sep=''), header=F)), nrow=L)
    z3 = raster(C, xmn=0, xmx=L, ymn=0, ymx=L)    
    coords <- xyFromCell(z3, which(z3[]>1e-2))
    p3 <- xyplot(coords[,2]~coords[,1], lwd=0.001, cex=0.001, col="red")
    p4 <- p2 + as.layer(p3)
    
    print(p1, position=c(0, 0.0, .5, 1), more=TRUE)
    print(p4, position=c(0.5, 0.0, 1, 1))
    dev.off()
}

