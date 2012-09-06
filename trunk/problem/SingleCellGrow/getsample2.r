Nsample = 1000
T = 8
r1 = 1
r2 = 0.4
d1 = 0
d2 = 0.0
p = 0.6
w = 1
x <- seq(Nsample)
for (i in seq(Nsample)){
  y <- twotypes(T, r1, r2, d1, d2, p, w)
  x[i] <- sum(y)
}
xx = c(x, -x)
yy <- density(xx)
yy$y <- 2*yy$y[yy$x>=0]
yy$x <- yy$x[yy$x>=0]
plot(yy)
