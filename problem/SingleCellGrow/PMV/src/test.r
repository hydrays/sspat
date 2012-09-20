dyn.load("simulator1.so")
source('simulator1n.r')
x <- simulator1n(8, 1.15, .6, 0.45, 0.3, 0.3, 0.75, 20000)
#cat(x, '\n')
Fe <- ecdf(x)
plot(Fe)

source('simulator1s.r')
for (i in seq(20000)){
  y <- simulator1s(8, 1.15, .6, 0.45, 0.3, 0.3, 0.75)
  x[i] = sum(y)
}
cat(x, '\n')
Fe <- ecdf(x)
lines(Fe, col='red')

source('simulator1p.r')
for (i in seq(20000)){
  temp <- runif(1)
  if (temp < 0.75) {
    z <- c(1, 0)
  } else {
    z <- c(0, 1)
  }
  y <- simulator1p(8, 1.15, .6, 0.45, 0.3, 0.3, z)
  x[i] <- sum(y)
}
cat(x, '\n')
Fe <- ecdf(x)
lines(Fe, col='blue')
