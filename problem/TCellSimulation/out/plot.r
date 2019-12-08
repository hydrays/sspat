require(lattice)
L  <- 256
colors <- c("#A7A7A7",
            "dodgerblue",
            "firebrick",
            "forestgreen",
            "gold")
myAt <- seq(5) - 1.5

for( i in seq(0, 2000) )
{
    cat(i, '\n')
    padded_i <- sprintf("%05d", i)
    png(paste("config_", padded_i, ".png", sep=''), height=600, width=1200)
    A <- matrix(unlist(read.csv(paste('c', padded_i, '.dat', sep=''), header=F)), nrow=L)
    p1  <- levelplot(A, at=myAt, col.regions = c('white', 'red', 'yellow', 'green'))

    B <- matrix(unlist(read.csv(paste('p', padded_i, '.dat', sep=''), header=F)), nrow=L)
    ##B  <- B
    ##p2  <- levelplot(B, seq(0, 1, length.out=100))
    p2  <- levelplot(B)

    print(p1, position=c(0, 0.0, .5, 1), more=TRUE)
    print(p2, position=c(0.5, 0.0, 1, 1))
    dev.off()
}

