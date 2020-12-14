colList <- c('red', 'green', 'blue', 'black', 'yellow')

dist <- NULL
png(paste('stat/dist.png', sep=''), height=600, width=1000)
plot(0, 0, type='n', xlim=c(0, 16), ylim=c(0, 0.5))
for ( i in seq(3) )
{
    padded_i <- sprintf("%02d", i)
    resultfilename <- paste('stat/stat_size', padded_i, '.csv', sep='')
    A <- read.csv(resultfilename)
    h <- hist(log2(A[,2]), breaks=8, plot=FALSE)
    lines(h$mids, h$counts/sum(h$counts), col=i)
}
legend("topright", col=seq(i), lty=1, legend=c("time point 1", "time point 2", "time point 3", "time point 4", "time point 5"))
dev.off()
