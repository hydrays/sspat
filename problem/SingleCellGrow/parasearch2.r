data7 <- read.csv('data7.csv')
mcell <- as.numeric(data7[8,3:ncol(data7[1,])])
mcell[is.na(mcell)] <- 0

load('ksdistance.RData')

Nsample = 10000
T = 8
x <- seq(Nsample)
dd <- seq(1000000)
j <- 1
n <- 0

for ( r1 in seq(0, 1, by=0.1) ){
  for ( r2 in seq(0, 0.5, by=0.1) ){
    for ( d1 in seq(0, 1, by=0.1) ){
      for ( d2 in seq(0, 0.5, by=0.1) ){
        for ( p in seq(0, 1, by=0.1) ){
          for ( w in seq(0, 1, by=0.1) ){
	    if ( (d[j] > 0.1) && (d[j] < 1) ){
	      n <- n + 1
              for (i in seq(Nsample)){
                y <- twotypes(T, r1, r2, d1, d2, p, w)
                x[i] <- sum(y)
              }
              dis <- ks.test(mcell, x)
              dd[j] <- dis$p.value
              cat(c(d[j], dd[j], r1, r2, d1, d2, p, w),'\n')
	    }
            j <- j+1
          }
        }
      }
    }
  }
}

