library("lattice")
library("grid")
jet.colors <- colorRampPalette(c("white", "red", "blue", "green"))

N = 1000
L = 1000
H = 200
pL = 600
pH = 50
.pwidth = 2048
.pheight = 576
for (i in seq(N)) {

  datafile <- sprintf("%s%05d%s", "m", i, ".dat")
  outfile <- sprintf("%s%05d%s", "slice", i, ".png")
  png(outfile, width=.pwidth, height=.pheight)
  z <- matrix(scan(datafile, n=L*H),
              L, H, byrow=TRUE)
  z <- z[1:pL, 1:pH]

  my.label.time <- sprintf("%s%d%s", "t = ", i, " (day)")
  p1 <- levelplot(z, col.regions=jet.colors,
            colorkey=FALSE, xlab="",
            ylab="",
            panel=function(...){
              panel.levelplot(...)
              grid.text(my.label.time,
                        y = unit(0.9, "npc"), gp=gpar(fontsize=30))
            },
            scales=list(cex=2))
  
  print(p1)
  ## print(p2)
  ## print(p1, position=c(0, 0.5, 1, 1), more=TRUE)
  ## print(p1, position=c(0, 0, 1, 0.5))

  ## --------------
  ## draw the customized legend
  ## --------------
  ## .xleft <- 
  cat(i, datafile, outfile, "\n")
  dev.off()
}
