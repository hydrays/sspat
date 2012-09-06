#jpeg('rplot.jpg')
Nsample = 2000
T = 8
r1 = 0
r2 = 0.6
d1 = 0
d2 = 0.3
p = 0
w = 0
x <- seq(Nsample)
for (i in seq(Nsample)){
  y <- simulator1(T, r1, r2, d1, d2, p, w)
  x[i] <- sum(y)
}

data7 <- read.csv('8dayslow.csv')
mcell <- as.numeric(data7[7,2:(ncol(data7[1,])-2)])
mcell[is.na(mcell)] <- 0

Fn <- ecdf(x)
Fe <- ecdf(mcell)
plot(Fn, pch='.', cex=5, xlim=c(0,300), col='green')
lines(Fe, col='black')

cat(mcell,'\n')
dis <- ks.test(mcell, x)
cat(dis$p.value)

T = 8
r1 = 0.75
r2 = 0.6 #Fix
d1 = 0.2
d2 = 0.3 #Fix
p = 0.2
w = 0.6
x <- seq(Nsample)
for (i in seq(Nsample)){
  y <- simulator1(T, r1, r2, d1, d2, p, w)
  x[i] <- sum(y)
}
data7 <- read.csv('8day.csv')
mcell <- as.numeric(data7[8,3:ncol(data7[1,])])
mcell[is.na(mcell)] <- 0

Fe <- ecdf(mcell)
Fn <- ecdf(x)
lines(Fn, pch='.', cex=5, xlim=c(0,300), col='blue')
lines(Fe, col='red')


cat(mcell,'\n')
dis <- ks.test(mcell, x)
cat(dis$p.value)

text(150, 0.5, 'step 1: fit the slow population:')
text(150, 0.4, 'r2 = 0.6, d2 = 0.3')

text(150, 0.3, 'step 2: fit the fast population using the above result:')
text(150, 0.2, 'r1 = 0.75, d1 = 0.2, p = 0.2, w = 0.6')

#dev.off()
