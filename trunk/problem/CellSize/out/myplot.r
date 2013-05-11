library("lattice")
library("grid")
jet.colors <- colorRampPalette(c("white", "red", "blue", "green"))

parainfo <- read.csv("control.csv", strip.white=TRUE)
NPool <- parainfo$VALUE[parainfo$PARAMETER=='NPool']
NCollect <- parainfo$VALUE[parainfo$PARAMETER=='NCollect']
NCollect2 <- parainfo$VALUE[parainfo$PARAMETER=='NCollect2']
.tend <- parainfo$VALUE[parainfo$PARAMETER=='tend']
.tpinc <- parainfo$VALUE[parainfo$PARAMETER=='tpinc']

N = 200
pL = 400
pH = 60
.pwidth = 800
.pheight = 640

cat("processing file ...[",N,"]\n")
i <- 190
datafile <- sprintf("%s%05d%s", "m", i, ".dat")

outfile <- sprintf("%s%05d%s", "slice", i, ".png")
png(outfile, width=.pwidth, height=.pheight)
#png(outfile, width=340, height=300)
z <- matrix(scan(datafile, n=NPool*5, quiet=TRUE),
            NPool, 5, byrow=TRUE)
my.label.time <- sprintf("%s%d%s", "t = ", i, " (day)")
p1 <- xyplot(200*(z[,1]-z[,2])/max(z[,1]-z[,2])~z[,1], xlim=c(0, 250000), grid=TRUE)
print(p1)
dev.off()

outfile <- sprintf("%s%05d%s", "fa", i, ".png")
png(outfile, width=.pwidth, height=.pheight)
#png(outfile, width=1080, height=640)
my.label.time <- sprintf("%s%d%s", "t = ", i, " (day)")
p1 <- histogram(z[,1], nint=30)
print(p1)
dev.off()


outfile <- sprintf("%s%05d%s", "rate", i, ".png")
png(outfile, width=.pwidth, height=.pheight)
#png(outfile, width=340, height=300)
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
  index <- floor((z[i,1]-lx)/dx)+1
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

## Plot histogram for mitosis cells.
datafile <- sprintf("CellMitosis.dat")
outfile <- sprintf("mitosis.png")
png(outfile, width=.pwidth, height=.pheight)
z <- matrix(scan(datafile, n=NPool*5, quiet=TRUE),
            NPool, 5, byrow=TRUE)
p1 <- histogram(z[,1], nint=30)
print(p1)
dev.off()

## Plot histogram for newborn cells.
datafile <- sprintf("CellNewborn.dat")
outfile <- sprintf("Newborn.png")
png(outfile, width=.pwidth, height=.pheight)
z <- matrix(scan(datafile, n=NPool*5, quiet=TRUE),
            NPool, 5, byrow=TRUE)
p1 <- histogram(z[,1], nint=30)
print(p1)
dev.off()

## Plot trajectory of the newborn cells syncronized.
outfile <- sprintf("syn.png")
png(outfile, width=.pwidth, height=.pheight)
ls <- 0
us <- 560000
T <- 24
M <- 80
ds <- (us-ls)/M
data = matrix(0, T, M)

for (i in seq(T)) {
  datafile <- sprintf("%s%05d%s", "s", i, ".dat")
  z <- matrix(scan(datafile, n=NCollect*5, quiet=TRUE),
              NCollect, 5, byrow=TRUE)
  for (j in seq(NCollect)){
    z[j,1] = min(us, z[j,1])
    z[j,1] = max(ls, z[j,1])         
    Mindex <- floor( (z[j,1]-ls)/ds ) + 1
    #print(Mindex)
    data[i, Mindex] <- data[i, Mindex] + 1
  }
}

M <- floor((280000 - ls)/ds + 1)
p1 <- levelplot(data[1:T, 1:M], col.regions=jet.colors,
                colorkey=TRUE, xlab="",
                aspect="fill",
                ylab="",
                scales=list(cex=2))
print(p1)
dev.off()

## Make figure 4
outfile <- sprintf("figure4.png")
png(outfile, width=.pwidth, height=.pheight)
ls <- 140000
us <- 260000
M <- 20
ds <- (us-ls)/M
data1 <- seq(M)
data2 <- seq(M)
x <- seq(M)
num1 <- seq(M)
num2 <- seq(M)

datafile1 <- sprintf("%s%02d%s", "prodiv", 6, ".dat")
datafile2 <- sprintf("%s%02d%s", "prodiv", 7, ".dat")
z1 <- matrix(scan(datafile1, n=NCollect2*1, quiet=TRUE),
            NCollect2, 1, byrow=TRUE)
z2 <- matrix(scan(datafile2, n=NCollect2*1, quiet=TRUE),
            NCollect2, 1, byrow=TRUE)
for (i in seq(M)){
  num1[i] <- 0
  num2[i] <- 0
}

for (j in seq(NCollect2)){
    if ( z1[j,1] > 0 ) {
      z1[j,1]  = min(us, z1[j,1])
      z1[j,1] = max(ls, z1[j,1])         
      Mindex <- floor( (z1[j,1]-ls)/ds ) + 1
      num1[Mindex] = num1[Mindex] + 1
      if (z2[j,1] > 0 ) {
        num2[Mindex] = num2[Mindex] + 1
      }
    }
  }
for ( i in seq(M) ) {
  x[i] = ls + i*ds - 0.5*ds
  data1[i] <- 1 - num2[i]/num1[i]
}

datafile1 <- sprintf("%s%02d%s", "prodiv", 9, ".dat")
datafile2 <- sprintf("%s%02d%s", "prodiv", 10, ".dat")
z1 <- matrix(scan(datafile1, n=NCollect2*1, quiet=TRUE),
            NCollect2, 1, byrow=TRUE)
z2 <- matrix(scan(datafile2, n=NCollect2*1, quiet=TRUE),
            NCollect2, 1, byrow=TRUE)
for (i in seq(M)){
  num1[i] <- 0
  num2[i] <- 0
}
for (j in seq(NCollect2)){
    if ( z1[j,1] > 0 ) {
      z1[j,1]  = min(us, z1[j,1])
      z1[j,1] = max(ls, z1[j,1])         
      Mindex <- floor( (z1[j,1]-ls)/ds ) + 1
      num1[Mindex] = num1[Mindex] + 1
      if (z2[j,1] > 0 ) {
        num2[Mindex] = num2[Mindex] + 1
      }
    }
  }
for ( i in seq(M) ) {
  data2[i] <- 1 - num2[i]/num1[i]
}
plot(x, data1, type = 'b', ylim = c(0, 1))
lines(x, data2, type = 'b', col='red', ylim=c(0,1))
dev.off()


## Trajecotries in the second paper
outfile <- sprintf("trace.png")
png(outfile, width=.pwidth, height=.pheight)
z <- matrix(scan("trace.dat", n=25*1000, quiet=TRUE),
            25, 1000, byrow=TRUE)
x <- matrix(0, 25, 1)
y <- matrix(0, 25, 1)
plot(c(0, 25000), c(0, 4000))
for ( j in seq(1000) ) {
  if ( (z[1, j] - 10000)^2 < 400){
  for (i in seq(25)) {
    x[i] <- z[i, j]
    if ( i < 2 ) {
      y[i] <- 0
    }
    else {
      y[i] <- z[i, j] - z[i-1, j]
    }
  }
  lines(x, y, type="b", col=j)
}
}

  
## for (i in seq(N)) {

##   datafile <- sprintf("%s%05d%s", "m", i, ".dat")
##   outfile <- sprintf("%s%05d%s", "slice", i, ".png")
##   png(outfile, width=.pwidth, height=.pheight)
##   z <- matrix(scan(datafile, n=L*H, quiet=TRUE),
##               L, H, byrow=TRUE)
##   z <- z[1:pL, 1:pH]

##   my.label.time <- sprintf("%s%d%s", "t = ", as.integer(i*.tpinc), " (day)")
##   p1 <- levelplot(z, col.regions=jet.colors,
##             colorkey=FALSE, xlab="",
##             ylab="",
##             panel=function(...){
##               panel.levelplot(...)
##               grid.text(my.label.time,
##                         y = unit(0.9, "npc"), gp=gpar(fontsize=30))
##             },
##             scales=list(cex=2))

##   print(p1)

##   ## print(p2)
##   ## print(p1, position=c(0, 0.5, 1, 1), more=TRUE)
##   ## print(p1, position=c(0, 0, 1, 0.5))

##   ## --------------
##   ## draw the customized legend
##   ## --------------
##   ## .xleft <-
##   output.str1 <- sprintf("%5d", i)
##   if (i > 1) cat("\b\b\b\b\b")
##   cat(output.str1)
##   dev.off()
## }
