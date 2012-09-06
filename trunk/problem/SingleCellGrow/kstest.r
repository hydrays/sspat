Nsample = 10000
T = 8
r1 = 0.8
r2 = 0.3
d1 = 0.1
d2 = 0.1
p = 0.0
w = 0.8
x <- seq(Nsample)
for (i in seq(Nsample)){
  y <- twotypes(T, r1, r2, d1, d2, p, w)
  x[i] <- sum(y)
}

data7 <- read.csv('data7.csv')
mcell <- as.numeric(data7[8,3:ncol(data7[1,])])
mcell[is.na(mcell)] <- 0

Fn <- ecdf(x)
Fe <- ecdf(mcell)

plot(Fn)
lines(Fe, col='red')

dis <- ks.test(mcell, x)
cat(dis$p.value)
