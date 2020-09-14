L  <- 64
dH <- 10
dW <- 10
pad <- 7

xData  <- NULL
yData  <- NULL
nCounter  <- 0

for( fileid in seq(10) )
{
    filename_x = sprintf('%s%05d%s', 'c', fileid, '.dat')
    filename_y = sprintf('%s%05d%s', 'c', fileid+1, '.dat')
    BW_x <- as.matrix(read.csv(filename_x, header=F))
    BW_y <- as.matrix(read.csv(filename_y, header=F))

    LH <- dim(BW_x)
    if ( LH[1] != L | LH[2] != L)
    {
        cat('wrong dimension\n')
        break
    }
    
    xlist = seq( pad+1, LH[1]-dW-pad+1, by=dW)
    ylist = seq( pad+1, LH[2]-dH-pad+1, by=dH)    

    for( i in xlist )
    {
        for( j in ylist )
        {
            nCounter = nCounter + 1
            B_x = BW_x[ (i-pad) : (i+dW+pad-1), (j-pad) : (j+dH+pad-1)]
            B_y = BW_y[ i : (i+dW-1), j : (j+dH-1)]
            xData <- rbind(xData, as.vector(t(B_x)))
            yData <- rbind(yData, as.vector(t(B_y)))
            ## montage({B_x, B_y});
            ## outname_x = sprintf('%s%08d%s', 'output/x_', figureid, '.png');
            ## outname_y = sprintf('%s%08d%s', 'output/y_', figureid, '.png');
            ## imwrite(B_x, outname_x);
            ## imwrite(B_y, outname_y);
            cat(fileid, i, j, nCounter, '\n')
        }
    }
}
write.table(xData, file='output/xData.csv', sep=',', row.names=F, col.names=F)
write.table(yData, file='output/yData.csv', sep=',', row.names=F, col.names=F)
