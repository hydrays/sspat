library("lattice")
library("grid")

parainfo <- read.csv("parameter.dat", strip.white=TRUE)
x <- matrix(scan("x.dat"))
.tend <- parainfo$VALUE[parainfo$PARAMETER=='tfinal']
.tpinc <- parainfo$VALUE[parainfo$PARAMETER=='tpinc']
.pm <- parainfo$VALUE[parainfo$PARAMETER=='p_m']
.vm <- parainfo$VALUE[parainfo$PARAMETER=='v_m']
MSC <- parainfo$VALUE[parainfo$PARAMETER=='alpha1']
MMC <- parainfo$VALUE[parainfo$PARAMETER=='alpha3']

my.label.diff1 <- sprintf("%s%7.2f", "M_SC = ", MSC)
my.label.diff2 <- sprintf("%s%6.3f", "M_MC = ", MMC)

N = 20000
.pwidth = 1024
.pheight = 768

NumCol <- 6
cat("processing file ...[",N,"]\n")
for (i in seq(N)) {

  datafile <- sprintf("%s%05d%s", "m", i, ".dat")
  outfile <- sprintf("%s%05d%s", "slice", i, ".png")
  png(outfile, width=.pwidth, height=.pheight)
  z <- matrix(scan(datafile, quiet=TRUE), ncol=NumCol, byrow=TRUE)
  my.label.time <- sprintf("%s%d%s", "t = ", as.integer(i*.tpinc), " (day)")

  p1 <- xyplot(c(0,2)~c(0, max(x)),
               type='n', xlab="", ylab="",
               scales=list(cex=2),
               key = list(x = 0.5, y=0.8,
               border=FALSE,
               lines=list(
               pch=c(1,1,NA,NA,NA,NA,NA,NA),
               col=c("red","black","blue","red","black",
               "grey","green","red","black"),
               type=c("b","b","l","l","l","l","l","l","l"),
               lty=c(1,1,1,2,2,3,3,3,3),
               lwd=c(1,1,4,4,4,6,0,0,0)),
               cex = 2,
               text = list(lab = c("SC","TC","P0","Net growth rate of SC",
                           "Net growth rate of MC","Pressure","","","")),
               columns = 1,
               #space = ""
               title = NULL
               ),
               panel=function(...){
                   panel.lines(x, z[,1], type='b', col='red')
                   panel.lines(x, z[,2], type='b', col='black')
                   panel.lines(x, z[,3], type='l',
                               lwd = 4, col='blue')
                   panel.lines(x, z[,4], type='l', lwd = 4,
                               lty = "dotted", col='red')
                   panel.lines(x, z[,5], type='l', lwd = 4,
                               lty = "dotted", col='black')
                   panel.lines(x, z[,6], type="l",
                               lwd = 6, lty="dotted", col='grey')
                   grid.text(my.label.time,
                             y = unit(0.9, "npc"), gp=gpar(fontsize=30))
                   grid.text(my.label.diff1,
                             y = unit(0.75, "npc"),
                             x = unit(0.3, "npc"),
                             gp=gpar(fontsize=24))
                   grid.text(my.label.diff2,
                             y = unit(0.7, "npc"),
                             x = unit(0.3, "npc"),
                             gp=gpar(fontsize=24))                   
                   
               }
               )
  print(p1)
  
  output.str1 <- sprintf("%5d", i)
  if (i > 1) cat("\b\b\b\b\b")
  cat(output.str1)
  dev.off()
}
