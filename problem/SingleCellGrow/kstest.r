Nsample = 1000
T = 8
r1 = 0.8
r2 = 0.3
p1 = c(0.2, 0.2, 0.2, 0.2, 0.2)
p2 = c(0.8, 0.2)
w = 1.0 # means all good cell
x <- seq(Nsample)
for (i in seq(Nsample)){
  y <- simulator2(T, r1, r2, p1, p2, w)
  x[i] <- sum(y)
}

data7 <- read.csv('8day.csv')
mcell <- as.numeric(data7[8,3:ncol(data7[1,])])
mcell[is.na(mcell)] <- 0
mcell <- mcell[mcell>=50]

Fn <- ecdf(x)
Fe <- ecdf(mcell)

plot(Fn)
lines(Fe, col='red')

dis <- ks.test(mcell, x)
cat(dis$p.value)
