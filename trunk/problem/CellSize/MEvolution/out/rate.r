library("lattice")
library("grid")

parainfo <- read.csv("control.csv", strip.white=TRUE)
NPool <- parainfo$VALUE[parainfo$PARAMETER=='NPool']
.tend <- parainfo$VALUE[parainfo$PARAMETER=='tend']
.tpinc <- parainfo$VALUE[parainfo$PARAMETER=='tpinc']

N = 200
.pwidth = 1064
.pheight = 600

i <- 190
datafile <- sprintf("%s%05d%s", "m", i, ".dat")
outfile <- sprintf("%s%05d%s", "rate", i, ".png")

png(outfile, width=340, height=300)
z <- matrix(scan(datafile, n=NPool*5, quiet=TRUE),
            NPool, 5, byrow=TRUE)
g <- 200*(z[,1]-z[,2])/max(z[,1]-z[,2])

index <- 0
Nbin <- 100
lx <- 50000
ux <- 250000
dx <- (ux-lx)/Nbin
m <- seq(Nbin)
nn <- seq(Nbin)
v <- seq(Nbin)
x <- seq(Nbin)
for (i in seq(Nbin)){
  x[i] <- lx + dx*(i - 0.5)
  m[i] <- 0
  v[i] <- 0
  nn[i] <- 0
}
for (i in seq(NPool)){
  z[i,1] = min(ux, z[i,1])
  z[i,1] = max(lx, z[i,1])     
  index <- floor((z[i,1]-lx)/dx)
  nn[index] <- nn[index] + 1
  v[index] <- v[index] + g[i]*g[i]
  m[index] <- m[index] + g[i]
}
for (i in seq(Nbin)){
  m[i] <- m[i] / nn[i]
  v[i] <- v[i] / nn[i] - m[i]*m[i]
}
p1 <- xyplot(m~x, xlim=c(0, ux), grid=TRUE)
print(p1)
dev.off()
