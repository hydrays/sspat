library("lattice")
library("grid")
jet.colors <- colorRampPalette(c("white", "red", "blue", "green"))

parainfo <- read.csv("control.csv", strip.white=TRUE)
.tend <- parainfo$VALUE[parainfo$PARAMETER=='tend']
.tpinc <- parainfo$VALUE[parainfo$PARAMETER=='tpinc']
.divide <- which(parainfo$PARAMETER=='useomp')
ompinfo <- parainfo[.divide:nrow(parainfo),]
parainfo <- parainfo[1:(.divide-1),]
  
N = 200
L = 1000
H = 200
pL = 600
pH = 100
.pwidth = 2048
.pheight = 576

cat("processing file ...[",N,"]\n")
i <- 0
datafile <- sprintf("%s%05d%s", "m", i, ".dat")
outfile <- sprintf("%s%05d%s", "slice", i, ".png")
png(outfile, width=.pwidth, height=.pheight)
z <- matrix(scan(datafile, n=L*H, quiet=TRUE),
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
                  grid.text(paste("------ # Parameter Setting",
                                  "( figure generated at",
                                   date(), " ) # ------"),
                            x = unit(0.01, "npc"), y = unit(0.8, "npc"),
                            just="left",
                            gp=gpar(fontsize=20) )                  
                  for ( j in seq(nrow(parainfo)) ){
                    grid.text(paste("# ", parainfo[j,1], " = ", parainfo[j,2]),
                              just="left",
                              x = unit(0.01 + 0.1*((j-1)%%4), "npc"),
                              y = unit(0.7-0.1*as.integer((j-1)/4), "npc"),
                              gp=gpar(fontsize=20) )
                  }

                  if (ompinfo[1,2]==1){
                    grid.text("------ # OMP in use! # ------",
                              x = unit(0.7, "npc"), y = unit(0.8, "npc"),
                              just="left",
                              gp=gpar(fontsize=20) )                  
                    for ( j in seq(nrow(ompinfo)) ){
                      grid.text(paste("# ", ompinfo[j,1], " = ", ompinfo[j,2]),
                                just="left",
                                x = unit(0.7 + 0.1*((j-1)%%2), "npc"),
                                y = unit(0.7-0.1*as.integer((j-1)/2), "npc"),
                                gp=gpar(fontsize=20) )
                    }
                  }
                },
                scales=list(cex=2))
print(p1)
#dev.off()

for (i in seq(N)) {

  datafile <- sprintf("%s%05d%s", "m", i, ".dat")
  outfile <- sprintf("%s%05d%s", "slice", i, ".png")
  png(outfile, width=.pwidth, height=.pheight)
  z <- matrix(scan(datafile, n=L*H, quiet=TRUE),
              L, H, byrow=TRUE)
  z <- z[1:pL, 1:pH]

  my.label.time <- sprintf("%s%d%s", "t = ", as.integer(i*.tpinc), " (day)")
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
  output.str1 <- sprintf("%5d", i)
  if (i > 1) cat("\b\b\b\b\b")
  cat(output.str1)
  dev.off()
}
