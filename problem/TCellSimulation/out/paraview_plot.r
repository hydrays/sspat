require(matlab)
require(raster)

##L <- 1024
L  <- 512
##L  <- 256
##L  <- 128
H <- 100

w = meshgrid(c(0, L/2, L), c(0, L/2, L), c(0, H/2, H), nargout=3)
g <- cbind(w$x, w$y, w$z)
colnames(g) = c("x", "y", "z")
write.csv(g, 'domain.csv', row.names=FALSE)  

sourceDIR = '/home/hydra/Tcell/sspat/problem/TCellSimulation/case_result3_1/out/'
baseDIR = paste(getwd(), '/', sep='')

#for( i in seq(0, 2000, by=10) )
i = 170
{
    cat(i, '\n')
    padded_i <- sprintf("%05d", i)
    A <- matrix(unlist(read.csv(paste(sourceDIR, 'c', padded_i, '.dat', sep=''), header=F)), nrow=L)
    B <- matrix(unlist(read.csv(paste(sourceDIR, 'phi', padded_i, '.dat', sep=''), header=F)), nrow=L)
    D <- matrix(unlist(read.csv(paste(sourceDIR, 'lambda', padded_i, '.dat', sep=''), header=F)), nrow=L)
    A[A==2] <- sample(c(0,2), length(A[A==2]), replace = TRUE, prob=c(0.9, 0.1))
    z1 = raster(A, xmn=0, xmx=L, ymn=0, ymx=L)
    z2 = raster(B, xmn=0, xmx=L, ymn=0, ymx=L)
    ##z2[z2[] < 1e-2] <- 1e-2
    chemokine <- cbind(xyFromCell(z2, which(z2[]>-100000)), 200.0*(z2[]), log10(z2[]))
    colnames(chemokine) = c("x", "y", "z", "h")
    TMcell <- cbind(xyFromCell(z1, which(z1[]==1 | z1[]==3)), H*seq(1, 1, length.out=length(z1[z1[]==1 | z1[]==3])), z1[z1[]==1 | z1[]==3])
    colnames(TMcell) = c("x", "y", "z", "t")
    inactiveTcell <- cbind(xyFromCell(z1, which(z1[]==2)), H*seq(1, 1, length.out=length(z1[z1[]==2])))
    ##runif(length(z1[z1[]==2]), min=90, max=100)*seq(1, 1, length.out=length(z1[z1[]==2])))
    colnames(inactiveTcell) = c("x", "y", "z")
    inactiveTcell <- inactiveTcell[!(inactiveTcell[,1] > L-10 | inactiveTcell[,1] < 10 | inactiveTcell[,2] > L-10 | inactiveTcell[,2] < 10),]
    ##upperPlane <- cbind(xyFromCell(z1, which(z1[]>-1)), 100*seq(1, 1, length.out=length(z1[])))
    ##colnames(upperPlane) = c("x", "y", "z")
    if (sum(sum(D)) > 1e-2)
    {
        z4 = raster(D, xmn=0, xmx=L, ymn=0, ymx=L)    
        FBcell <- cbind(xyFromCell(z4, which(z4[]>1e-4)), 4.0*log10(z2[z4[]>1e-4]), z4[z4[]>1e-4])
        colnames(FBcell) = c("x", "y", "z", "t")
        write.csv(FBcell, 'FBcell.csv', row.names=FALSE)
    }
    write.csv(chemokine, 'chemokine.csv', row.names=FALSE)    
    write.csv(TMcell, 'TMcell.csv', row.names=FALSE)
    write.csv(inactiveTcell, 'inactiveTcell.csv', row.names=FALSE)
    ##write.csv(upperPlane, 'upperPlane.csv', row.names=FALSE)        

    cmdpv = paste("python draw_script.py ", baseDIR, sep='')
    system(cmdpv)
    
    #outfilename <- paste("config_", padded_i, ".png", sep='')
    outfilename <- "nosignal.png"
    cmdmv = paste("mv pvoutput.png ", outfilename, sep='')
    system(cmdmv)
}

