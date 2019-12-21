library(png)
require(matlab)
require(raster)
require(lattice)
require(rasterVis)
require(latticeExtra)
library(viridisLite)
coul <- viridis(100)
require("RColorBrewer")
gray.cols <- colorRampPalette(c('black', 'white'))(n=100)

## ##color Scheme I
## chemokine.cols <- colorRampPalette(c('white', '#F8766D'))(n=100)
## fb.col <- "red"
## tcell.col <- "#1F78B4"

## Color Scheme II
##chemokine.cols <- colorRampPalette(c('white', '#1F78B4'))(n=100)
##fb.col <- "#0000FF"
##tcell.col <- "#F8766D"

## ## Color Scheme III
## chemokine.cols <- colorRampPalette(c('white', '#FFFF00'))(n=100)
## fb.col <- "#1F78B4"
## tcell.col <- "#F8766D"

## ## Color Scheme IV
## chemokine.cols <- colorRampPalette(c('white', '#7CAE00'))(n=100)
## fb.col <- "#1F78B4"
## tcell.col <- "#F8766D"

## ## Color Scheme V
## chemokine.cols <- colorRampPalette(c('white', '#F7931E'))(n=100)
## fb.col <- "#1F78B4"
## tcell.col <- "#F8766D"

## Color Scheme VI
chemokine.cols <- colorRampPalette(c('white', '#FFFF00'))(n=100)
fb.col <- "#DEF1FC"
tcell.col <- "#F8766D"

highsig.col <- "#AFDCE8"
weaksig.col <- "#DEF1FC"
dark.col <- "#9FA0A0"
white.col <- "#FFFFFF"

fbhigh.col <- "#1F78B6"
fbweak.col <- "#95DEEF"

L <- 1024
Lmin <- 300
Lmax <- 700
##L  <- 512
##L  <- 256
#L  <- 128
H <- 100

w = meshgrid(c(0, L/2, L), c(Lmin, (Lmin+Lmax)/2, Lmax), c(-18, H/2, H), nargout=3)
g <- cbind(w$x, w$y, w$z)
colnames(g) = c("x", "y", "z")
write.csv(g, 'domain.csv', row.names=FALSE)  

sourceDIR = '/home/hydra/Tcell/sspat/problem/TCellSimulation/case_xx03_2/out/'
baseDIR = paste(getwd(), '/', sep='')

A0 <- matrix(unlist(read.csv('f00000.dat', header=F)), nrow=L)
zz1 = raster(A0[Lmin:Lmax, ], xmn=0, xmx=L, ymn=0, ymx=L)

##for( i in seq(1, 2000, by=10) )
i = 160
{
    cat(i, '\n')
    padded_i <- sprintf("%05d", i)
    A <- matrix(unlist(read.csv(paste(sourceDIR, 'c', padded_i, '.dat', sep=''), header=F)), nrow=L)
    B <- matrix(unlist(read.csv(paste(sourceDIR, 'phi', padded_i, '.dat', sep=''), header=F)), nrow=L)
    D <- matrix(unlist(read.csv(paste(sourceDIR, 'lambda', padded_i, '.dat', sep=''), header=F)), nrow=L)
    A[A==2] <- sample(c(0,2), length(A[A==2]), replace = TRUE, prob=c(0.98, 0.02))
    z1 = raster(A[Lmin:Lmax, ], xmn=0, xmx=L, ymn=Lmin-1, ymx=Lmax)
    z2 = raster(B[Lmin:Lmax, ], xmn=0, xmx=L, ymn=Lmin-1, ymx=Lmax)
    z4 = raster(D[Lmin:Lmax, ], xmn=0, xmx=L, ymn=Lmin-1, ymx=Lmax)    
    z2[z2[] < 1e-2] <- 1e-2
    chemokine <- cbind(xyFromCell(z2, which(z2[]>-100000)), 4.0*log10(z2[]), log10(z2[]), zz1[])
    colnames(chemokine) = c("x", "y", "z", "h", "floor")
    TMcell <- cbind(xyFromCell(z1, which(z1[]==1 | z1[]==3)), H*seq(1, 1, length.out=length(z1[z1[]==1 | z1[]==3])), z1[z1[]==1 | z1[]==3])
    colnames(TMcell) = c("x", "y", "z", "t")
    inactiveTcell <- cbind(xyFromCell(z1, which(z1[]==2)), H*seq(1, 1, length.out=length(z1[z1[]==2])))
    ##runif(length(z1[z1[]==2]), min=90, max=100)*seq(1, 1, length.out=length(z1[z1[]==2])))
    colnames(inactiveTcell) = c("x", "y", "z")
    inactiveTcell <- inactiveTcell[!(inactiveTcell[,1] > L-10 | inactiveTcell[,1] < 10 | inactiveTcell[,2] > L-10 | inactiveTcell[,2] < 10),]
    ##upperPlane <- cbind(xyFromCell(z1, which(z1[]>-1)), 100*seq(1, 1, length.out=length(z1[])))
    ##colnames(upperPlane) = c("x", "y", "z")
    FBcell <- cbind(xyFromCell(z4, which(z4[]>1e-2)), z2[z4[]>1e-2], z4[z4[]>1e-2])
    colnames(FBcell) = c("x", "y", "z", "h")

    write.csv(chemokine, 'chemokine.csv', row.names=FALSE)    
    write.csv(TMcell, 'TMcell.csv', row.names=FALSE)
    write.csv(inactiveTcell, 'inactiveTcell.csv', row.names=FALSE)
    ##write.csv(upperPlane, 'upperPlane.csv', row.names=FALSE)        
    write.csv(FBcell, 'FBcell.csv', row.names=FALSE)

    z1[z1==2] <- 0
    z1[z1==3] <- 0
    p2 <- levelplot(z1, col.regions = c("white", "black"), scales=list(cex=2), colorkey=FALSE, margin=FALSE)
    coords <- xyFromCell(z1, which(z1[]==1))
    pp2 <- xyplot(coords[,2]~coords[,1], labels=list(cex=5), lwd=1, cex=1.5, col="black")    
    p2 <- p2 + as.layer(pp2)

    #p2l1 <- levelplot(z1, col.regions = c(white.col), margin=FALSE, colorkey=FALSE)
    #coords <- xyFromCell(z1, which(z1[]==1))
    #p2l2 <- xyplot(coords[,2]~coords[,1], lwd=1, cex=0.5, col=dark.col)    
    #p2 <- p2l1 + as.layer(p2l2)

    png("final_cut.png", width=1024, height=401)
    print(p2)
    dev.off()	

    cmdpv = paste("python sym_draw_script.py ", baseDIR, sep='')
    system(cmdpv)
    
    png("combine.png", width=800, height=800)
    plot(0:2000, 0:2000, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    rasterImage(readPNG(source="pvoutput.png"), 0, 800, 2000, 2000)
    rasterImage(readPNG(source="init_cut.png"), 0, 250, 600, 550)
    rasterImage(readPNG(source="final_cut.png"), 600, 0, 2000, 800)
    dev.off()

    #print(p1, position=c(0, 0.7, 1, 1), more=TRUE)
    #print(p2, position=c(0, 0.0, 0.5, 0.3), more=TRUE)
    #print(p2, position=c(0.5, 0.0, 1, 0.3)

    
    outfilename <- paste("config_", padded_i, ".png", sep='')
    #outfilename <- "bot_right.png"
    #cmdmv = paste("mv pvoutput.png ", outfilename, sep='')
    cmdmv = paste("mv combine.png ", outfilename, sep='')
    system(cmdmv)
}

