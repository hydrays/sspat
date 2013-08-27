library("lattice")
library("grid")

y0 <- c(0, 1100)
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
lambda1 <- 20000
gamma1 <- 10.0
lambda2 <- 0.25
gamma2 <- 0.15
k <- 1.0
for ( i in seq(2, N) ) {
  a[i] <- i*dt
  y1[i,] <- y[i-1,]
  y[i,1] <- y1[i,1] + dt * (lambda1*(k*a[i])^4/(1+(k*a[i])^4) - gamma1*y1[i,1])
  y[i,2] <- y1[i,2] + dt * (lambda2*min(y1[i,1], y1[i,2]) - gamma2*y1[i,2])
  z[i] <- min(y[i, 1], y[i, 2])
  s[i] <- max(s[i-1], y[i,2])
}

pdf("ode1.pdf")
p1 <- xyplot(z~a, xlim=c(0, 12),
             ylim=c(0, 2500),
             xlab=list("cell age (hours)", cex = 1.5),
             ylab=list("amount of mRNA and Ribsome (arbitrary unit)", cex=1.5),
             type = "l",
             lwd = 15,
             lty = 1,
             col = "red",
             scales = list(cex = 1.5),
             grid=TRUE,
             panel=function(...){
               panel.xyplot(...)
               panel.lines(a, y[,1], lwd=4, type='l', lty = 4, col='black')
               panel.lines(a, y[,2], lwd=4, type='l', lty = 1, col='blue')
               panel.lines(a, y[,2]*y[,2]/(1000*a), lwd=4, type='l', lty = 1, col='green')
             },)
print(p1)
dev.off()
#lines(a, y[,2], col="red")

v <- matrix(0, N, 1)
v[1:N-1] <- s[2:N]
v <- v - s
v[N] <- v[N-1]
v <- v/dt
pdf("ode2.pdf")
p2 <- xyplot(v[10:N]~s[10:N], xlim=c(0, 2500),
             xlab=list("s (cell size, fl)", cex = 1.5),
             ylab=list("growth rate (fl/hour)", cex=1.5),
             type = "l",
             lwd = 4,
             lty = 1,
             col = "red",
             scales = list(cex = 1.5),
             panel=function(...){
               panel.xyplot(...)
               panel.lines(seq(0, 1000, 1000/50), seq(0, 100, 100/50),
                           lwd=2, type='l', lty=2, col='black')
             },)
print(p2)
dev.off()
