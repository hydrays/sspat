library("lattice")
library("grid")
#jet.colors <- colorRampPalette(c("blue", "yellow", "red"))
#jet.colors <- colorRampPalette()
jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

parainfo <- read.csv("control.csv", strip.white=TRUE)
NPool <- parainfo$VALUE[parainfo$PARAMETER=='NPool']
NCollect <- parainfo$VALUE[parainfo$PARAMETER=='NCollect']
NCollect2 <- parainfo$VALUE[parainfo$PARAMETER=='NCollect2']
.tend <- parainfo$VALUE[parainfo$PARAMETER=='tend']
.tpinc <- parainfo$VALUE[parainfo$PARAMETER=='tpinc']
.tminc <- parainfo$VALUE[parainfo$PARAMETER=='tminc']

lambda1 <- 10000
gamma1 <- 5.0
lambda2 <- 0.25
gamma2 <- 0.15
k <- 1.0

L = 6
.pwidth = 800
.pheight = 600
i <- 490

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
index1 <- 0
index2 <- 0
for ( j in seq(2, N) ) {
  a[j] <- j*dt
  y1[j,] <- y[j-1,]
  y[j,1] <- y1[j,1] + dt * (lambda1*(k*a[j])^4/(1+(k*a[j])^4) - gamma1*y1[j,1])
  y[j,2] <- y1[j,2] + dt * (lambda2*min(y1[j,1], y1[j,2]) - gamma2*y1[j,2])
  z[j] <- min(y[j, 1], y[j, 2])
  s[j] <- max(s[j-1], y[j,2])
  if ( s[j] == y[j,2] && index1 == 0) {
    index1 <- j
  }
  if ( y[j,1] < y[j,2]  && index2 == 0 && index1 != 0 ){
    index2 <- j
  }
}
pdf("ode1.pdf", width=7, height=5)
p1 <- xyplot(z~a, xlim=c(0, 12),
             ylim=c(0, 2500),
             xlab=list("cell age (hours)", cex = 2),
             ylab=list("s (cell size, fl)", cex=2),
             type = "l",
             lwd = 8,
             lty = 6,
             col = "red",
             scales = list(cex = 2),
             key = list(x = 0.2, y=0.26,
               border=TRUE,
               lines=list(
                 col=c("black", "blue", 'red', 'green'),
                 type=c("l","l","l","l"),
                 lty=c(4,1,6,1),
                 lwd=c(4,4,4,4)),
               cex = 1.4,
               text = list(lab = c("mRNA","ribosome",
                             "loaded ribosome", 'cell size')),
               columns = 1,
               #space = ""
               title = NULL
               ),
             panel=function(...){
               panel.xyplot(...)
               auto.key = TRUE
               panel.lines(a, s, lwd=4, type='l', lty=1, col='green')
               panel.lines(a, y[,1], lwd=4, type='l', lty = 4, col='black')
               panel.lines(a, y[,2], lwd=4, type='l', lty = 1, col='blue')
               panel.lines(c(a[index1], a[index1]),
                           c(0, s[index1]), lwd=4, type='l',
                           lty = 2, col='grey')
               panel.lines(c(a[index2], a[index2]),
                           c(0, s[index2]), lwd=4, type='l',
                           lty = 2, col='grey')
               grid.text('I',
                         just="left",
                         x = unit(0.09, "npc"),
                         y = unit(0.08, "npc"),
                         gp=gpar(fontsize=20) )
               grid.text('II',
                         just="left",
                         x = unit(0.43, "npc"),
                         y = unit(0.42, "npc"),
                         gp=gpar(fontsize=20) )
               grid.text('III',
                         just="left",
                         x = unit(0.84, "npc"),
                         y = unit(0.6, "npc"),
                         gp=gpar(fontsize=20) )                              
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
for ( j in seq(2, N) ) {
  a[j] <- j*dt
  y1[j,] <- y[j-1,]
  y[j,1] <- y1[j,1] + dt * (lambda1*(k*a[j])^4/(1+(k*a[j])^4) - gamma1*y1[j,1])
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
datafile <- sprintf("%s%05d%s", "m", i, ".dat")
outfile <- sprintf("%s%05d%s", "slice", i, ".png")
png(outfile, width=.pwidth, height=.pheight)
#pdf(outfile, width=.pwidth, height=.pheight)
z <- matrix(scan(datafile, n=NPool*L, quiet=TRUE),
            NPool, L, byrow=TRUE)
p1 <- xyplot((z[,3]-z[,2])/.tminc~(z[,3]+z[,2])/2, xlim=c(0, 2500), grid=TRUE,
             panel=function(...){
               panel.xyplot(...)
               panel.lines(s, v, 
                           lwd=2, type='l', lty=2, col='black')
               },)
print(p1)
dev.off()
s.save <- s

#######################################
## Asynchronized cell size distribution
#######################################
#outfile <- sprintf("%s%05d%s", "fa", i, ".png")
#png(outfile, width=.pwidth, height=.pheight)
#pdf(outfile, width=.pwidth, height=.pheight)
pdf("fasyn.pdf", width=7, height=5)
datafile <- sprintf("%s%05d%s", "m", i, ".dat")
z <- matrix(scan(datafile, n=NPool*L, quiet=TRUE),
            NPool, L, byrow=TRUE)

y.lim = c(c(0, 1.5e-3))
y.at <- pretty(y.lim)
y.labels <- formatC(1000*y.at, format = "g")
#y.labels[length(y.at)] <- paste(y.labels[length(y.at)], 'x10-3')

trellis.par.set(clip=list(panel = "off")) 
p1 <- histogram(z[,1], nint=30, xlim=c(250, 3000),
                ylim = y.lim,
                type = c("density"),
                xlab=list("s (cell size, fl)", cex = 1.5),
                ylab=list("proportion of cells", cex=1.5),
                scales=list(cex=1.5, y=list(at = y.at,
                                       labels = y.labels)),
                #mtext("10-3")
                panel = function(x, ...){
                  panel.histogram(x, ...)
                  panel.densityplot(x, ...,
                                    lwd = 4, col="red")
                  panel.axis(side = c("top"),
                             at = c(350),
                             labels = c(expression(x10^{-3})),
                             ticks=FALSE,
                             outside=TRUE,
                             rot=c(0,90),
                             text.cex=1.4)
                })
print(p1)
dev.off()


########################
### Averaged Growth Rate
########################
#png("rate.png", width=.pwidth, height=.pheight)
trellis.par.set(clip=list(panel = "on"))
pdf("rate.pdf", width=7, height=6)
g <- (z[,3]-z[,2])/.tminc
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
             xlab=list("s (cell size, fl)", cex = 1.5),
             ylab=list("growth rate (fl/hour)", cex=1.5),
             scales=list(cex=1.5),
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

trellis.par.set(clip=list(panel = "off")) 
####################################
## Plot histogram for mitosis cells.
####################################
pdf("mitosis.pdf", width=7, height=5.5)
#png("mitosis.png", width=.pwidth, height=.pheight)
datafile <- sprintf("CellMitosis.dat")
z <- matrix(scan(datafile, n=NPool*L, quiet=TRUE),
            NPool, L, byrow=TRUE)
y.lim = c(0, 1.4e-3)
y.at <- pretty(y.lim)
y.labels <- formatC(1000*y.at, format = "g")
p1 <- histogram(z[,1], nint=30, xlim=c(300, 3000),
                ylim=y.lim,
                type = c("density"),
                xlab=list("s (cell size, fl)", cex = 1.5),
                ylab=list("proportion of cells", cex=1.5),
                scales=list(cex=1.5, y=list(at = y.at,
                                       labels = y.labels)),
                panel = function(x, ...){
                  panel.histogram(x, ...)
                  panel.densityplot(x, ...,
                                    lwd = 4, col="red")
                  panel.axis(side = c("top"),
                             at = c(400),
                             labels = c(expression(x10^{-3})),
                             ticks=FALSE,
                             outside=TRUE,
                             rot=c(0,90),
                             text.cex=1.4)
                })
print(p1)
dev.off()

####################################
## Plot histogram for newborn cells.
####################################
pdf("newborn.pdf", width=7, height=5)
#png(outfile, width=.pwidth, height=.pheight)
datafile <- sprintf("CellNewborn.dat")
z <- matrix(scan(datafile, n=NPool*L, quiet=TRUE),
            NPool, L, byrow=TRUE)
y.lim = c(0, 3e-3)
y.at <- pretty(y.lim)
y.labels <- formatC(1000*y.at, format = "g")
p1 <- histogram(z[,1], nint=30, xlim=c(250, 3000),
                ylim=y.lim,
                type = c("density"),
                xlab=list("s (cell size, fl)", cex = 1.5),
                ylab=list("proportion of cells", cex=1.5),
                scales=list(cex=1.5, y=list(at = y.at,
                                       labels = y.labels)),
                panel = function(x, ...){
                  panel.histogram(x, ...)
                  panel.densityplot(x, ...,
                                    lwd = 4, col="red")
                  panel.axis(side = c("top"),
                             at = c(350),
                             labels = c(expression(x10^{-3})),
                             ticks=FALSE,
                             outside=TRUE,
                             rot=c(0,90),
                             text.cex=1.4)
                })
print(p1)
dev.off()

####################################
## Plot trajectory of the newborn cells syncronized.
####################################
outfile <- sprintf("syn.png")
pdf("syn.pdf", width=7, height=5.5)
#png(outfile, width=.pwidth, height=.pheight)
ls <- 0
us <- 5600
T <- 24
M <- 160
ds <- (us-ls)/M
data = matrix(0, T, M)

for (i in seq(T)) {
  print(i)
  datafile <- sprintf("%s%05d%s", "s", i, ".dat")
  z <- matrix(scan(datafile, n=NCollect*L, quiet=TRUE),
              NCollect, L, byrow=TRUE)
  for (j in seq(NCollect)){
    z[j,1] = min(us, z[j,1])
    z[j,1] = max(ls, z[j,1])         
    Mindex <- floor( (z[j,1]-ls)/ds ) + 1
    data[i, Mindex] <- data[i, Mindex] + 1
  }
}

M <- floor((2800 - ls)/ds + 1)
lM <- floor((400 - ls)/ds) 
p1 <- levelplot(data[1:T, lM:M]/NCollect,
                  col.regions=jet.colors(1000),
                cuts=1000,
                useRaster=TRUE,
                region=FALSE,
                colorkey=FALSE,
                xlab=list("time (hours)", cex=2),
                ylab=list("s (cell size, fl)", cex=2),
                aspect="fill",
                scales=list(cex=1.5,
                  y = list(relation="sliced",
                    at=list(c(M-60-lM, M-40-lM, M-20-lM)),
                    label=list(c(700, 1400, 2100))),
                x = list(relation="sliced",
                  at=list(c(6, 12, 18))))
                )
print(p1)
dev.off()

## Make figure 4
#png("fig4.png", width=.pwidth, height=.pheight)
pdf("fig4.pdf", width=7, height=7)
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

p1 <- xyplot(data1[1:(M-3)]~x[1:(M-3)],
             xlab = list("s (cell size, fl)", cex=2),
             ylab = list("Proportion of divisions", cex=2),
             scales=list(cex=2),
             type = 'b', lwd = 2, col="blue",
             ylim = c(0, 1),
             xlim = c(1400, 2600),
             grid=TRUE,
             cex = 1.5,
             pch = 20,
             panel = function(...){
               panel.xyplot(...)
               panel.lines(x[1:(M-3)], data2[1:(M-3)], 
                         lwd=2, type='b', col='red', cex=1.5, pch=20)
             }
             )
print(p1)
dev.off()


## ## Trajecotries in the second paper
## outfile <- sprintf("trace.png")
## png(outfile, width=.pwidth, height=.pheight)
## z <- matrix(scan("trace.dat", n=25*1000/.tminc, quiet=TRUE),
##             25/.tminc, 1000, byrow=TRUE)
## x <- matrix(0, 25/.tminc, 1)
## y <- matrix(0, 25/.tminc, 1)
## plot(c(0, 2500), c(0, 200))
## for ( j in seq(1000) ) {
##   if ( (z[1, j] - 1000)^2 < 1000){
##   for (i in seq(25/.tminc)) {
##     x[i] <- z[i, j]
##     if ( i < 2 ) {
##       y[i] <- 0
##     }
##     else {
##       y[i] <- (z[i, j] - z[i-1, j])/.tminc
##     }
##   }
##   lines(x, y, type="b", col=j)
## }
## }

## ## Generate lifespan
## outfile <- sprintf("lifespan.png")
## png(outfile, width=.pwidth, height=.pheight)
## z <- matrix(scan("lifespan.dat", n=NCollect2*1, quiet=TRUE),
##               NCollect2, 1, byrow=TRUE)
## p1 <- histogram(z, nint=30, scales=list(cex=2))
## print(p1)
## dev.off()
