require(matlab)
require(raster)

##L <- 1024
L  <- 512
##L  <- 256
##L  <- 128

output_file_prefix <- "sample01"
output_path = "out"

cmd0 <- "rm -r stat"
system(cmd0)

baseDir <- "stat"
if ( !dir.exists(baseDir) )
{
    dir.create(baseDir)
}

tList <- c(90, 170, 250)
##tList <- c(60, 90, 120, 150, 260)
##for( i in seq(0, 2000, by=10) )
for ( i in seq(length(tList)) )
##i = 400
{
    padded_i <- sprintf("%02d", i)
    outDir <- paste(baseDir, '/timepoint', padded_i, sep='')
    if ( !dir.exists(outDir) )
    {
        dir.create(outDir)
    }
    for( rep in seq(20) )
    {
        cat("\n ***********run********* ", rep, "\n")
        sourceFolderName <- paste("../case_", output_file_prefix, "_", as.character(rep), '/out/', sep='')
        ##sourceFolderName <- '../out/'

        padded_i <- sprintf("%05d", tList[i])
        A <- matrix(unlist(read.csv(paste(sourceFolderName, 'c', padded_i, '.dat', sep=''), header=F)), nrow=L)
        A[A>1] <- 0
        A[1:10, ] <- 1
        A[(L-10):L, ] <- 1
        A[,1:10] <- 1
        A[,(L-10):L] <- 1
        z1 = raster(A, xmn=0, xmx=L, ymn=0, ymx=L)
        z = aggregate(z1, 4, fun=max)
        z = 1-z

        padded_rep <- sprintf("%04d", rep)
        outputFileName <- paste(outDir, '/brick_', padded_rep, '.tif', sep='')
        writeRaster(z, filename=outputFileName, format="GTiff", overwrite=TRUE)
        ##write.csv(Tcell, 'Tcell.csv', row.names=FALSE)
        outputFileName <- paste(outDir, '/brick_', padded_rep, '.png', sep='')
        png(outputFileName)
        plot(z)
        dev.off()
    }
}

cat('begin matlab ...\n')
cmd1 <- "matlab -nodisplay -nosplash -nodesktop -r \"run('Stat.m'); exit;\" "
system(cmd1)
cat('end matlab ...\n')

cat('begin drawing ...\n')
cmd2 <- "Rscript Draw.r"
system(cmd2)
cat('end drawing ...\n')
