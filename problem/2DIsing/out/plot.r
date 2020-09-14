require(lattice)
require(rasterVis)
require(latticeExtra)
##library(viridisLite)
##coul <- viridis(100)
require("RColorBrewer")
##color Scheme I
chemokine.cols <- colorRampPalette(c('white', '#F8766D'))(n=100)
fb.col <- "red"
tcell.col <- "#1F78B4"

L  <- 64
##L  <- 256
#L  <- 128
colors <- c("#A7A7A7",
            "dodgerblue",
            "firebrick",
            "forestgreen",
            "gold")

for( i in seq(0, 2000, by=1) )
##i = 0
{
    cat(i, '\n')
    padded_i <- sprintf("%05d", i)
    png(paste("config_", padded_i, ".png", sep=''), height=600, width=1200)
    ##pdf(paste("config_", padded_i, ".pdf", sep=''), height=6, width=12)
    A <- matrix(unlist(read.csv(paste('c', padded_i, '.dat', sep=''), header=F)), nrow=L)
    z1 = raster(A, xmn=0, xmx=L, ymn=0, ymx=L)
    p1  <- levelplot(z1, margin=FALSE, colorkey=FALSE)
    print(p1)
    dev.off()
}

