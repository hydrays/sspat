data7 <- read.csv('data7.csv')
mcell <- as.numeric(data7[8,3:ncol(data7[1,])])
mcell[is.na(mcell)] <- 0

load('ksdistance.RData')
load('ksdist2.RData')
j <- 1
n <- 0
m <- 0
dvsr1 <- seq(1000000)
dvsr2 <- seq(1000000)
dvsd1 <- seq(1000000)
dvsd2 <- seq(1000000)
dvsp <- seq(1000000)
dvsw <- seq(1000000)
pl <- seq(1000000)
for ( r1 in seq(0, 1, by=0.1) ){
  for ( r2 in seq(0, 0.5, by=0.1) ){
    for ( d1 in seq(0, 1, by=0.1) ){
      for ( d2 in seq(0, 0.5, by=0.1) ){
        for ( p in seq(0, 1, by=0.1) ){
          for ( w in seq(0, 1, by=0.1) ){
	    if ( (d[j] > 0.1) && (d[j] < 1) ){
              n <- n+1
              if ( (d[j] > 0.1) && (d[j] < 1) ){
                m <- m+1
                dvsr1[m] <- r1
                dvsr2[m] <- r2
                dvsd1[m] <- d1
                dvsd2[m] <- d2
                dvsp[m] <- p
                dvsw[m] <- w
                pl[m] <- d[j]
                cat(m, j, '\n')
              }
            }
            j <- j+1
          }
        }
      }
    }
  }
}

par(mfrow=c(3,2))
plot(dvsr1[1:n], pl[1:m])
plot(dvsr2[1:m], pl[1:m])
plot(dvsd1[1:m], pl[1:m])
plot(dvsd2[1:m], pl[1:m])
plot(dvsp[1:m], pl[1:m])
plot(dvsw[1:m], pl[1:m])

