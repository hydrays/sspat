library("lattice")
library("grid")
#jet.colors <- colorRampPalette(c("blue", "yellow", "red"))
#jet.colors <- colorRampPalette()
jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

parainfo <- read.csv("control.csv", strip.white=TRUE)
NPool <- parainfo$VALUE[parainfo$PARAMETER=='NPool']
.tend <- parainfo$VALUE[parainfo$PARAMETER=='tend']
.tpinc <- parainfo$VALUE[parainfo$PARAMETER=='tpinc']
.tminc <- parainfo$VALUE[parainfo$PARAMETER=='tminc']

lambda1 <- 2000
gamma1 <- 1.0
lambda2 <- 0.25
gamma2 <- 0.15
k <- 0.5

L = 6
.pwidth = 800
.pheight = 600
i <- 5

########################
### Averaged Growth Rate
########################

datafile <- sprintf("%s%05d%s", "m", i, ".dat")
z <- matrix(scan(datafile, n=NPool*L, quiet=TRUE),
            NPool, L, byrow=TRUE)
##png("rate.png", width=.pwidth, height=.pheight)
trellis.par.set(clip=list(panel = "on"))
pdf("rate.pdf", width=7, height=6.3)
g <- (z[,3]-z[,2])/.tminc
g <- pmax(0, g)
index <- 25
Nbin <- 50
lx <- 0
ux <- 2525
dx <- (ux-lx)/Nbin
m <- seq(Nbin)
nn <- seq(Nbin)
va <- seq(Nbin)
x <- seq(Nbin)
for (i in seq(Nbin)){
  x[i] <- lx + (i-0.5)*dx
  m[i] <- 0
  va[i] <- 0
  nn[i] <- 0
}
for (i in seq(NPool)){
  s <- (z[i,3] + z[i,2])/2
  #s <- z[i,2]
  s <- min(ux, s)
  s <- max(lx, s)     
  index <- floor((s-lx)/dx)+1
  nn[index] <- nn[index] + 1
  va[index] <- va[index] + g[i]*g[i]
  m[index] <- m[index] + g[i]
}
for (i in seq(Nbin)){
  m[i] <- m[i] / nn[i]
  va[i] <- va[i] / nn[i] - m[i]*m[i]
}
p1 <- xyplot(m~x, xlim=c(0, 2500), grid=TRUE,
             lty = 1, type='b', col='black',
             lwd = 2,
             xlab=list("s (cell size, fl)", cex = 2),
             ylab=list("growth rate (fl/hour)", cex=2),
             scales=list(cex=2),
             panel = function(...){
               panel.xyplot(...)
               panel.lines(s.save, v, 
                           lwd=2, type='l',
                           lty=2, col='red', xlim=c(0, 2500))
             })
#panel.lines(s, v, 
#                         lwd=2, type='l', lty=2, col='black')
print(p1)
dev.off()

