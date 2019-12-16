require(lattice)
require(rasterVis)
require(latticeExtra)
library(viridisLite)
coul <- viridis(100)

##L <- 1024
L  <- 512
##L  <- 256
#L  <- 128
colors <- c("#A7A7A7",
            "dodgerblue",
            "firebrick",
            "forestgreen",
            "gold")
myAt <- seq(5) - 1.5

C0 <- matrix(unlist(read.csv(paste('a00001.dat', sep=''), header=F)), nrow=L)
totalRateList <- NULL
totalDemandList <- NULL
for( i in seq(0, 2000, by=10) )
{
    cat(i, '\n')
    padded_i <- sprintf("%05d", i)
    png(paste("config_", padded_i, ".png", sep=''), height=600, width=900)
    A <- matrix(unlist(read.csv(paste('c', padded_i, '.dat', sep=''), header=F)), nrow=L)
    B <- matrix(unlist(read.csv(paste('phi', padded_i, '.dat', sep=''), header=F)), nrow=L)
    C <- matrix(unlist(read.csv(paste('a', padded_i, '.dat', sep=''), header=F)), nrow=L)
    ##A[A==2] <- sample(c(0,2), length(A[A==2]), replace = TRUE, prob=c(0.9, 0.1))

    
    totalRate <- sum( (A==0)*C )
    totalDemand <- sum( (A==0)*(1+B)*C0 )    
    totalRateList <- c(totalRateList, totalRate)
    totalDemandList <- c(totalDemandList, totalDemand)    
    pp2 <- xyplot(totalDemandList ~ seq(length(totalDemandList)))
    p2 <- xyplot(totalRateList ~ seq(length(totalRateList)), col='red')
    ## z2 = raster(B, xmn=0, xmx=L, ymn=0, ymx=L)
    ## p2 <- levelplot(z2, col.regions = coul, margin=FALSE,
    ##                 panel=function(..., at, contour, region) {
    ##                     panel.levelplot(..., contour = FALSE)
    ##                     panel.contourplot(..., at=c(1), contour = TRUE, col="white", lty=2, region = FALSE)
    ##                 }
    ##                 )
    A[A==2] = 0
    z1 = raster(A, xmn=0, xmx=L, ymn=0, ymx=L)
    p1  <- levelplot(z1, at=myAt, col.regions = c('white', '#C9CACA', '#F8766D', '#F8766D'), margin=FALSE,
                     colorkey=FALSE)

    p2 <- pp2 + as.layer(p2)

    z3 = raster(C, xmn=0, xmx=L, ymn=0, ymx=L)
    p3 <- levelplot(z3, col.regions = coul, margin=FALSE,
                    panel=function(..., at, contour, region) {
                        panel.levelplot(..., contour = FALSE)
                        panel.contourplot(..., at=c(0.5), contour = TRUE, col="white", lty=2, region = FALSE)
                        ##panel.dotplot(y~x)
                    }
                    )
    ## C <- matrix(unlist(read.csv(paste('lambda', padded_i, '.dat', sep=''), header=F)), nrow=L)
    ## z3 = raster(C, xmn=0, xmx=256, ymn=0, ymx=256)    
    ## coords <- xyFromCell(z3, which(z3[]>1e-2))
    ## p3 <- xyplot(coords[,2]~coords[,1], lwd=0.5, pch=1, cex=0.5, col="red")
    ## p4 <- p2 + as.layer(p3)
    
    print(p1, position=c(0, 0.0, .5, 1), more=TRUE)
    print(p2, position=c(0.5, 0.5, 1, 1), more=TRUE)
    print(p3, position=c(0.5, 0, 1, 0.5))
    dev.off()
}

