require('lattice')
require('grid')

lambda1 <- 2000
gamma1 <- 1.0
lambda2 <- 0.25
gamma2 <- 0.15
k <- 0.5

y0 <- c(0, 1000)
T <- 15.0
N <- 1000
dt <- T/N
y <- matrix(0, N, 2)
y1 <- matrix(0, N, 2)
y[1,] <- y0
a <- seq(N)
z <- seq(N)
z[1] <- min(y[1, 1], y[1, 2])
a[1] <- 0
index1 <- 0
index2 <- 0
for ( j in seq(2, N) ) {
  a[j] <- j*dt
  y1[j,] <- y[j-1,]
  y[j,1] <- y1[j,1] + dt * (lambda1*(k*a[j])^4/(1+(k*a[j])^4) - gamma1*y1[j,1])
  y[j,2] <- y1[j,2] +
      dt * max(0, (lambda2*min(y1[j,1], y1[j,2]) - gamma2*y1[j,2]))
  z[j] <- min(y[j, 1], y[j, 2])
  if ( y[j,1] > y[j,2] && index1 == 0) {
    index1 <- j-1
  }
  if ( y[j,1] < y[j,2]  && index2 == 0 && index1 != 0 ){
    index2 <- j
  }
}
pdf("ode1.pdf", width=7, height=5)
p1 <- xyplot(z~a, xlim=c(0, 13),
             ylim=c(0, 3000),
             xlab=list("cell age (hours)", cex = 1.5),
             ylab=list("s (cell size, fl)", cex=1.5),
             type = "l",
             lwd = 8,
             lty = 6,
             col = "red",
             scales = list(cex = 1.5, x=list(at=c(2,4,6,8,10,12,14))),
             key = list(x = 0.1, y=0.95,
               border=TRUE,
               lines=list(
                 col=c("black", "blue", 'red'),
                 type=c("l","l","l"),
                 lty=c(4,1,6),
                 lwd=c(4,4,4)),
               cex = 1.4,
               text = list(lab = c("mRNA","cell size",
                             "min(m, s)")),
               columns = 1,
               #space = ""
               title = NULL
               ),
             panel=function(...){
               panel.xyplot(...)
               auto.key = TRUE
               panel.lines(a, y[,1], lwd=4, type='l', lty = 4, col='black')
               panel.lines(a, y[,2], lwd=4, type='l', lty = 1, col='blue')
               panel.lines(c(a[index1], a[index1]),
                           c(0, y[index1,2]), lwd=4, type='l',
                           lty = 2, col='grey')
               panel.lines(c(a[index2], a[index2]),
                            c(0, y[index2,2]), lwd=4, type='l',
                            lty = 2, col='grey')
               grid.text('I',
                         just="left",
                         x = unit(0.1, "npc"),
                         y = unit(0.2, "npc"),
                         gp=gpar(fontsize=20) )
               grid.text('II',
                         just="left",
                         x = unit(0.45, "npc"),
                         y = unit(0.3, "npc"),
                         gp=gpar(fontsize=20) )
               grid.text('III',
                         just="left",
                         x = unit(0.84, "npc"),
                         y = unit(0.4, "npc"),
                         gp=gpar(fontsize=20) )                              
             },)
print(p1)
dev.off()
