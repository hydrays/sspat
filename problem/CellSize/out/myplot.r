library("lattice")
library("grid")

parainfo <- read.csv("control.csv", strip.white=TRUE)
NPool <- parainfo$VALUE[parainfo$PARAMETER=='NPool']
.tend <- parainfo$VALUE[parainfo$PARAMETER=='tend']
.tpinc <- parainfo$VALUE[parainfo$PARAMETER=='tpinc']

N = 200
pL = 400
pH = 60
.pwidth = 1064
.pheight = 600

cat("processing file ...[",N,"]\n")
i <- 190
datafile <- sprintf("%s%05d%s", "m", i, ".dat")

outfile <- sprintf("%s%05d%s", "slice", i, ".png")
png(outfile, width=340, height=300)
z <- matrix(scan(datafile, n=NPool*5, quiet=TRUE),
            NPool, 5, byrow=TRUE)
my.label.time <- sprintf("%s%d%s", "t = ", i, " (day)")
p1 <- xyplot(200*(z[,1]-z[,2])/max(z[,1]-z[,2])~z[,1], xlim=c(0, 250000), grid=TRUE)
print(p1)
dev.off()

outfile <- sprintf("%s%05d%s", "fa", i, ".png")
png(outfile, width=1080, height=640)
my.label.time <- sprintf("%s%d%s", "t = ", i, " (day)")
p1 <- histogram(z[,1], nint=30)
print(p1)
dev.off()


outfile <- sprintf("%s%05d%s", "rate", i, ".png")
png(outfile, width=340, height=300)
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
