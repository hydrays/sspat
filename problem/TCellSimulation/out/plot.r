require(lattice)
require(rasterVis)
require(latticeExtra)
library(viridisLite)
coul <- viridis(100)

L  <- 256
colors <- c("#A7A7A7",
            "dodgerblue",
            "firebrick",
            "forestgreen",
            "gold")
myAt <- seq(5) - 1.5

for( i in seq(0, 2000, by=10) )
{
    cat(i, '\n')
    padded_i <- sprintf("%05d", i)
    png(paste("config_", padded_i, ".png", sep=''), height=600, width=1200)
    A <- matrix(unlist(read.csv(paste('c', padded_i, '.dat', sep=''), header=F)), nrow=L)
    A[A==2] <- sample(c(0,2), length(A[A==2]), replace = TRUE, prob=c(0.9, 0.1))
    z = raster(A, xmn=0, xmx=256, ymn=0, ymx=256)
    p1  <- levelplot(z, at=myAt, col.regions = c('white', '#C9CACA', '#F8766D', '#F8766D'), margin=FALSE,
                     colorkey=FALSE)

    B <- matrix(unlist(read.csv(paste('phi', padded_i, '.dat', sep=''), header=F)), nrow=L)
    coords <- xyFromCell(z, which(z[]==3))
    zz = raster(B, xmn=0, xmx=256, ymn=0, ymx=256)
    p2 <- levelplot(zz, col.regions = coul, margin=FALSE,
                    panel=function(..., at, contour, region) {
                        panel.levelplot(..., contour = FALSE)
                        panel.contourplot(..., at=c(0.5), contour = TRUE, col="white", lty=2, region = FALSE)
                        ##panel.dotplot(y~x)
                    }
                    )
    p3 <- xyplot(coords[,2]~coords[,1], lwd=0.5, pch=1, cex=0.5, col="red")
    p4 <- p2 + as.layer(p3)
    
    print(p1, position=c(0, 0.0, .5, 1), more=TRUE)
    print(p4, position=c(0.5, 0.0, 1, 1))
    dev.off()
}

