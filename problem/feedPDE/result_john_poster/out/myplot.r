library("lattice")
library("grid")

parainfo <- read.csv("parameter.dat", strip.white=TRUE)
x <- matrix(scan("x.dat"))
.tend <- parainfo$VALUE[parainfo$PARAMETER=='tfinal']
.tpinc <- parainfo$VALUE[parainfo$PARAMETER=='tpinc']
.pm <- parainfo$VALUE[parainfo$PARAMETER=='p_m']
.vm <- parainfo$VALUE[parainfo$PARAMETER=='v_m']

N = 20000
.pwidth = 1024
.pheight = 768

NumCol <- 8
cat("processing file ...[",N,"]\n")
for (i in seq(N)) {

  datafile <- sprintf("%s%05d%s", "m", i, ".dat")
  outfile <- sprintf("%s%05d%s", "slice", i, ".png")
  png(outfile, width=.pwidth, height=.pheight)
  z <- matrix(scan(datafile, quiet=TRUE), ncol=NumCol, byrow=TRUE)
  my.label.time <- sprintf("%s%d%s", "t = ", 400+as.integer(i*.tpinc), " (day)")

  p1 <- xyplot(c(0,2)~c(0, max(x)),
               type='n', xlab="", ylab="",
               scales=list(cex=2),
               key = list(x = 0.8, y=0.9,
               border=FALSE,
               lines=list(
               pch=c(1,1,1,NA,NA,NA,NA,NA),
               col=c("red","green","black","pink","grey",
               "grey","green","red","black"),
               type=c("b","b","b","l","l","l","l","l","l"),
               lty=c(1,1,1,2,2,3,3,3,3),
               lwd=c(1,1,1,4,6,0,0,0,0)),
               cex = 2,
               text = list(lab = c("SC","TC","MC","TGFb",
                           "pressure","","","","")),
               columns = 1,
               #space = ""
               title = NULL
               ),
               panel=function(...){
                   panel.lines(x, z[,1], type='b', col='red')
                   panel.lines(x, z[,2], type='b', col='green')
                   panel.lines(x, z[,3], type='b', col='black')
                   panel.lines(x, z[,4], type='l',
                               lwd = 4, lty = 2, col='pink')
                   #panel.lines(x, z[,5], type='l',
                   #            lwd = 1, lty=1, col='blue')
                   #panel.lines(x, z[,6], type='l',
                   #            lwd = 2, lty=2, col='red')
                   #panel.lines(x, z[,7], type='l',
                   #            lwd = 2, lty=2, col='black')
                   panel.lines(x, z[,8], type="l",
                               lwd = 6, lty="dotted", col='grey')
                   grid.text(my.label.time,
                             y = unit(0.9, "npc"), gp=gpar(fontsize=30))
               }
               )
  print(p1)
  
  output.str1 <- sprintf("%5d", i)
  if (i > 1) cat("\b\b\b\b\b")
  cat(output.str1)
  dev.off()
}
