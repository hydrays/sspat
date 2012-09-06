source('goodcell.r')
Nsample <-100
T = 10
p = 0.2
x <- seq(Nsample)
for (i in seq(100)) x[i] <- goodcell(T, p)