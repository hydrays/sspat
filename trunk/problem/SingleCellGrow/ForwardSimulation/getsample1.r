source('goodcell.r')
source('../plfit.r')
Nsample <-100000
T = 13
p = 0.4
x <- seq(Nsample)
for (i in seq(Nsample)) x[i] <- goodcell(T, p)

write.table(x, file='../out/xFor02.csv', sep = ",", row.names = FALSE, col.names=FALSE)
#para <- as.numeric(plfit(x))
write.table(para, file='../out/paraFor02.csv', sep = ",", row.names = FALSE, col.names=FALSE)


y <- x
yy <- c(y, -y)
dens <- density(yy)
dens$y <- 2*dens$y[dens$x >=0]
dens$x <- dens$x[dens$x >= 0]
plot(dens,col="red")

y <- rexp(10000)*mean(x)
yy <- c(y, -y)
dens <- density(yy)
dens$y <- 2*dens$y[dens$x >=0]
dens$x <- dens$x[dens$x >= 0]
lines(dens,col="gray")
