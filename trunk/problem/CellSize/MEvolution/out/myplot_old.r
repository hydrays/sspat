library("lattice")
library("grid")
jet.colors <- colorRampPalette(c("white", "red", "blue", "green"))

parainfo <- read.csv("control.csv", strip.white=TRUE)
NPool <- parainfo$VALUE[parainfo$PARAMETER=='NPool']
NCollect <- parainfo$VALUE[parainfo$PARAMETER=='NCollect']
NCollect2 <- parainfo$VALUE[parainfo$PARAMETER=='NCollect2']
.tend <- parainfo$VALUE[parainfo$PARAMETER=='tend']
.tpinc <- parainfo$VALUE[parainfo$PARAMETER=='tpinc']
.tminc <- parainfo$VALUE[parainfo$PARAMETER=='tminc']

N = 200
L = 6
pL = 400
pH = 60
.pwidth = 800
.pheight = 640

cat("processing file ...[",N,"]\n")
i <- 190
datafile <- sprintf("%s%05d%s", "m", i, ".dat")

# Generate the growth curve of mRNA and Ribosome
y0 <- c(0, 1000)
T <- 15.0
N <- 1000
dt <- T/N
y <- matrix(0, N, 2)
y1 <- matrix(0, N, 2)
y[1,] <- y0
a <- seq(N)
s <- seq(N)
z <- seq(N)
a[1] <- 0
s[1] = y[1, 2]
lambda1 <- 10000
gamma1 <- 5.0
lambda2 <- 0.25
gamma2 <- 0.15
k <- 1.0
for ( j in seq(2, N) ) {
  a[j] <- j*dt
  y1[j,] <- y[j-1,]
  y[j,1] <- y1[j,1] + dt * (lambda1*(k*a[j])^2/(1+(k*a[j])^2) - gamma1*y1[j,1])
  y[j,2] <- y1[j,2] + dt * (lambda2*min(y1[j,1], y1[j,2]) - gamma2*y1[j,2])
  z[j] <- min(y[j, 1], y[j, 2])
  s[j] <- max(s[j-1], y[j,2])
}
pdf("ode1.pdf")
p1 <- xyplot(z~a, xlim=c(0, 12),
             ylim=c(0, 2500),
             xlab=list("cell age (hours)", cex = 1.5),
             ylab=list("amount of mRNA and Ribsome (arbitrary unit)", cex=1.5),
             type = "l",
             lwd = 15,
             lty = 1,
             col = "red",
             scales = list(cex = 1.5),
             panel=function(...){
               panel.xyplot(...)
               panel.lines(a, y[,1], lwd=4, type='l', lty = 4, col='black')
               panel.lines(a, y[,2], lwd=4, type='l', lty = 1, col='blue')
             },)
print(p1)
dev.off()

# Generate the statistics of growth rate vs cell size
y0 <- c(0, 1000)
T <- 15.0
N <- 1000
dt <- T/N
y <- matrix(0, N, 2)
y1 <- matrix(0, N, 2)
y[1,] <- y0
a <- seq(N)
s <- seq(N)
z <- seq(N)
a[1] <- 0
s[1] = y[1, 2]
lambda1 <- 10000
gamma1 <- 5.0
lambda2 <- 0.25
gamma2 <- 0.15
k <- 1.0
for ( j in seq(2, N) ) {
  a[j] <- j*dt
  y1[j,] <- y[j-1,]
  y[j,1] <- y1[j,1] + dt * (lambda1*(k*a[j])^2/(1+(k*a[j])^2) - gamma1*y1[j,1])
  y[j,2] <- y1[j,2] + dt * (lambda2*min(y1[j,1], y1[j,2]) - gamma2*y1[j,2])
  z[j] <- min(y[j, 1], y[j, 2])
  s[j] <- max(s[j-1], y[j,2])
}
v <- matrix(0, N, 1)
v[1:N-1] <- s[2:N]
v <- v - s
v[N] <- v[N-1]
v <- v/dt
s[1:N/2] <- seq(0, 1000, 1000/N/2)
v[1:N/2] <- seq(0, 100, 100/N/2)
outfile <- sprintf("%s%05d%s", "slice", i, ".png")
png(outfile, width=.pwidth, height=.pheight)
#png(outfile, width=340, height=300)
z <- matrix(scan(datafile, n=NPool*L, quiet=TRUE),
            NPool, L, byrow=TRUE)
my.label.time <- sprintf("%s%d%s", "t = ", i, " (day)")
p1 <- xyplot((z[,3]-z[,2])/.tminc~(z[,3]+z[,2])/2, xlim=c(0, 2500), grid=TRUE,
             panel=function(...){
               panel.xyplot(...)
               panel.lines(s, v, 
                           lwd=2, type='l', lty=2, col='black')
               },)
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
g <- (z[,3]-z[,2])/.tminc
index <- 0
Nbin <- 100
lx <- 500
ux <- 2500
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
  s <- (z[i,3] + z[i,2])/2
  s <- min(ux, s)
  s <- max(lx, s)     
  index <- floor((s-lx)/dx)+1
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
z <- matrix(scan(datafile, n=NPool*L, quiet=TRUE),
            NPool, L, byrow=TRUE)
p1 <- histogram(z[,1], nint=30)
print(p1)
dev.off()

## Plot histogram for newborn cells.
datafile <- sprintf("CellNewborn.dat")
outfile <- sprintf("Newborn.png")
png(outfile, width=.pwidth, height=.pheight)
z <- matrix(scan(datafile, n=NPool*L, quiet=TRUE),
            NPool, L, byrow=TRUE)
p1 <- histogram(z[,1], nint=30)
print(p1)
dev.off()

## Plot trajectory of the newborn cells syncronized.
outfile <- sprintf("syn.png")
png(outfile, width=.pwidth, height=.pheight)
ls <- 0
us <- 5600
T <- 24
M <- 100
ds <- (us-ls)/M
data = matrix(0, T, M)

for (i in seq(T)) {
  datafile <- sprintf("%s%05d%s", "s", i, ".dat")
  z <- matrix(scan(datafile, n=NCollect*L, quiet=TRUE),
              NCollect, L, byrow=TRUE)
  for (j in seq(NCollect)){
    z[j,1] = min(us, z[j,1])
    z[j,1] = max(ls, z[j,1])         
    Mindex <- floor( (z[j,1]-ls)/ds ) + 1
    #print(Mindex)
    data[i, Mindex] <- data[i, Mindex] + 1
  }
}

M <- floor((2800 - ls)/ds + 1)
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
ls <- 1400
us <- 2600
M <- 20
ds <- (us-ls)/M
data1 <- seq(M)
data2 <- seq(M)
x <- seq(M)
num1 <- seq(M)
num2 <- seq(M)

datafile1 <- sprintf("%s%02d%s", "prodiv", 7, ".dat")
datafile2 <- sprintf("%s%02d%s", "prodiv", 8, ".dat")
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

datafile1 <- sprintf("%s%02d%s", "prodiv", 10, ".dat")
datafile2 <- sprintf("%s%02d%s", "prodiv", 11, ".dat")
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
z <- matrix(scan("trace.dat", n=25*1000/.tminc, quiet=TRUE),
            25/.tminc, 1000, byrow=TRUE)
x <- matrix(0, 25/.tminc, 1)
y <- matrix(0, 25/.tminc, 1)
plot(c(0, 2500), c(0, 200))
for ( j in seq(1000) ) {
  if ( (z[1, j] - 1000)^2 < 1000){
  for (i in seq(25/.tminc)) {
    x[i] <- z[i, j]
    if ( i < 2 ) {
      y[i] <- 0
    }
    else {
      y[i] <- (z[i, j] - z[i-1, j])/.tminc
    }
  }
  lines(x, y, type="b", col=j)
}
}
