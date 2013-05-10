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
png(outfile, width=.pwidth, height=.pheight)
z <- matrix(scan(datafile, n=NPool*5, quiet=TRUE),
            NPool, 5, byrow=TRUE)
my.label.time <- sprintf("%s%d%s", "t = ", i, " (day)")
p1 <- xyplot(z[,1]-z[,2]~z[,1])
print(p1)
dev.off()

outfile <- sprintf("%s%05d%s", "fa", i, ".png")
png(outfile, width=.pwidth, height=.pheight)
my.label.time <- sprintf("%s%d%s", "t = ", i, " (day)")
p1 <- histogram(z[,1], nint=30)
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
