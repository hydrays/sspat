library("lattice")
library("grid")
jet.colors <- colorRampPalette(c("white", "red", "blue", "green", "black"))

parainfo <- read.csv("control.csv", strip.white=TRUE)
L <- parainfo$VALUE[parainfo$PARAMETER=='L']
H <- parainfo$VALUE[parainfo$PARAMETER=='H']
.tend <- parainfo$VALUE[parainfo$PARAMETER=='tend']
.tpinc <- parainfo$VALUE[parainfo$PARAMETER=='tpinc']
.divide <- which(parainfo$PARAMETER=='useomp')
ompinfo <- parainfo[.divide:nrow(parainfo),]
parainfo <- parainfo[1:(.divide-1),]

N = 2000
pL = 400
pH = 80
#.pwidth = 2048/2
#.pheight = 288
.pwidth = 1028
.pheight = 256

dataH <- H + 8

cat("processing file ...[",N,"]\n")
i <- 0
datafile <- sprintf("%s%05d%s", "m", i, ".dat")
outfile <- sprintf("%s%05d%s", "slice", i, ".png")
png(outfile, width=.pwidth, height=.pheight)
z <- matrix(scan(datafile, n=L*dataH, quiet=TRUE),
            L, dataH, byrow=TRUE)
z <- z[200:pL-200, 1:pH]
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
  #outfile <- sprintf("%s%05d%s", "slice", i, ".png")
  #png(outfile, width=.pwidth, height=.pheight)
  outfile <- sprintf("%s%05d%s", "slice", i, ".pdf")
  pdf(outfile, width=.pwidth/60, height=.pheight/60)
  
  z <- matrix(scan(datafile, n=L*dataH, quiet=TRUE),
              L, dataH, byrow=TRUE)
  nutri <- 50*z[, H+1]
  p0 <- 50*z[,H+7]
  v0 <- 50*z[,H+8]  
                                        #z <- z[200:pL-200, 1:pH]
  z <- z[1:pL, 1:pH]
  if(max(z)<4){
    z[200, pH] <- 4
    #z[pL, pH] <- 4
  }

  my.label.time <- sprintf("%s%d%s", "t = ", as.integer(i*.tpinc), " (day)")
  p1 <- levelplot(z, col.regions=jet.colors,
                  colorkey=FALSE, xlab="",
                  ylab="",
                  ylim = c(0,60),
                  draw=FALSE,
            panel=function(...){
              panel.levelplot(...)
              panel.lines(seq(pL), nutri, lwd=4, type='l', col=colors()[100])
              panel.lines(seq(pL), p0, lwd=4, type='l', col=colors()[450])
              panel.lines(seq(pL), v0, lwd=4, type='l', col=colors()[68])    
              #panel.lines(seq(pL), 5, lwd=2, lty=2, col='yellow')
              #grid.text(my.label.time,
              #           x = unit(0.85, "npc"),
              #           y = unit(0.85, "npc"), gp=gpar(fontsize=30))
            },
            scales=list(cex=2, y=list(relation="free",
                                 at=list(c(0, 25, 50)),
                                 labels=list(c(0, 0.5, 1)))))

  #p2 <- xyplot(nutri~seq(pL))
  print(p1)

  ## p2 <- xyplot(c(0,0)~c(400, 60),
  ##                 ylab="",
  ##                 ylim = c(0,60),
  ##              xlim = c(0,400),
  ##                 draw=FALSE,
  ##           panel=function(...){
  ##             panel.lines(seq(pL), nutri, lwd=4, type='l', col=colors()[100])
  ##             panel.lines(seq(pL), p0, lwd=4, type='l', col=colors()[450])
  ##             panel.lines(seq(pL), v0, lwd=4, type='l', col=colors()[68])    
  ##             #panel.lines(seq(pL), 5, lwd=2, lty=2, col='yellow')
  ##             grid.text(my.label.time,
  ##                       x = unit(0.85, "npc"),
  ##                       y = unit(0.85, "npc"), gp=gpar(fontsize=30))
  ##           },
  ##           scales=list(cex=2, y=list(relation="free",
  ##                                at=list(c(0, 25, 50)),
  ##                                labels=list(c(0, 0.5, 1)))))

  
  ## print(p2)
  #print(p1, position=c(0, 0.5, 1, 1), more=TRUE)
  #print(p2, position=c(0, 0, 1, 0.5))

  ## --------------
  ## draw the customized legend
  ## --------------
  ## .xleft <-
  output.str1 <- sprintf("%5d", i)
  if (i > 1) cat("\b\b\b\b\b")
  cat(output.str1)
  dev.off()
}
