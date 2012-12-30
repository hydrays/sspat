library("lattice")

Nplot = 500

z <- matrix(scan('tr.dat', n=Nplot*3, quiet=TRUE),
            Nplot, 3, byrow=TRUE)

p1 <- xyplot(c(0,max(z[,3]/800))~c(0, max(25)),
             type='o',
             xlab="",
             ylab=list(label='transcription level',cex=2),
             xlim=c(0, 50),
             #ylim=c(-0.1, 0.8),
             scales=list(cex=2),
             panel=function(...){
               panel.lines(z[,1], z[,3]/1000, type='s', cex=1.2, lwd=2, col='red')
             },
             )

p2 <- xyplot(c(0,max(1))~c(0, max(25)),
             type='o',
             xlab=list(label='t (hour)',cex=2),
             ylab=list(label="signal", cex=2),
             xlim=c(0, 25),
             ylim=c(-0.2, 1.2),
             scales=list(y=list(relation="free",
               cex=2, at=c(0,1), labels=c("off","on")),
               x=list(cex=2)),
             panel=function(...){
               panel.lines(z[,1], z[,2], type='s', lwd=2, col='black')
             },
             )



print(p1, position=c(0, 0.4, 1, 1), more=TRUE)
print(p2, position=c(0, 0, 1, 0.4))

dev.copy(pdf, height=6, width=8, 'fig_signal_s1.pdf')
dev.off()
