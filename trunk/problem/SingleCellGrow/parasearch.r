data7 <- read.csv('data7.csv')
mcell <- as.numeric(data7[8,3:ncol(data7[1,])])
mcell[is.na(mcell)] <- 0

Nsample = 1000
T = 8
d <- seq(1000000)
x <- seq(Nsample)
j <- 1
for ( r1 in seq(0, 1, by=0.1) ){
  for ( r2 in seq(0, 0.5, by=0.1) ){
    for ( d1 in seq(0, 1, by=0.1) ){
      for ( d2 in seq(0, 0.5, by=0.1) ){
        for ( p in seq(0, 1, by=0.1) ){
          for ( w in seq(0, 1, by=0.1) ){
            for (i in seq(Nsample)){
              y <- twotypes(T, r1, r2, d1, d2, p, w)
              x[i] <- sum(y)
            }
            dis <- ks.test(mcell, x)
            d[j] <- dis$p.value
            cat(c(d[j], r1, r2, d1, d2, p, w),'\n')
            j <- j+1
          }
        }
      }
    }
  }
}

