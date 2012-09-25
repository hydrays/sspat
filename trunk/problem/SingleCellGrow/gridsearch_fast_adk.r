library('adk')
dyn.load('simulator1.so')
source('simulator1n.r')

mcell <- as.matrix(read.csv('data/day8fast.csv'))

Nsample = 10000
T = 8
r2 = 0.652
d2 = 0.34
d <- double(1000000)
pvalue <- double(1000000)
x <- seq(Nsample)
j <- 1
for ( r1 in seq(0.5, 1.5, by=0.05) ){
    for ( d1 in seq(0, 1, by=0.05) ){
        for ( v in seq(0, 1, by=0.05) ){
          for ( w in seq(0, 1, by=0.05) ){
            x <- simulator1n(T, r1, r2, d1, d2, v, w, Nsample)
            dis <- adk.test(mcell, x)
            d[j] <- dis$adk[1,2]
            pvalue[j] <- dis$adk[2,2]
            cat(c(d[j], r1, r2, d1, d2, v, w),'\n')
            j <- j+1
          }
      }
  }
}


j <- 1
n <- 0
m <- 0
dvsr1 <- seq(1000000)
dvsd1 <- seq(1000000)
dvsv <- seq(1000000)
dvsw <- seq(1000000)
pl <- seq(1000000)
for ( r1 in seq(0.5, 1.5, by=0.05) ){
    for ( d1 in seq(0, 1, by=0.05) ){
        for ( v in seq(0, 1, by=0.05) ){
          for ( w in seq(0, 1, by=0.05) ){
	    if ( (d[j] > 0.54) && (d[j] < 1) ){
              n <- n+1
              if ( (d[j] > 0.1) && (d[j] < 1) ){
                m <- m+1
                dvsr1[m] <- r1
                dvsd1[m] <- d1
                dvsv[m] <- v
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

par(mfrow=c(2,2))
plot(dvsr1[1:n], pl[1:m])
plot(dvsd1[1:m], pl[1:m])
plot(dvsv[1:m], pl[1:m])
plot(dvsw[1:m], pl[1:m])

## ## -----------------------------------
## ## Plot the fit
## ## -----------------------------------

## par(mfrow=c(1,1))
## Nsample = 10000
## T = 8
## r1 = 1.1
## r2 = 0.652
## d1 = 0.4
## d2 = 0.34
## v = 0.25
## w = 0.4
## x <- simulator1n(T, r1, r2, d1, d2, v, w, Nsample)

## Fn <- ecdf(x)
## Fe <- ecdf(mcell)
## plot(Fn, xlim=c(0,300),
##      main = "ECDF: simulation vs data",
##      ylab = "Cumulative probability",
##      xlab = "8-day Cell number ")
## lines(Fe, col='red')

## text(200, 0.6, "Goodness of fit")
## text(200, 0.55, "Black: simulated using parameters from grid search")
## text(200, 0.5, "Red: experiment data of the mixed population")
## text(200, 0.45, "Result (r1=1, r2=0.6, d1=0.4, d2=0.3, v=0, w=0.3):")
## text(200, 0.4, "KS-distance: 0.08")
## text(200, 0.35, "p-value: 0.81")

## dis <- adk.test(mcell, x)
## cat(dis$adk[1,2])

## mtext(paste("@", date()), side=1, line=4, adj=1.04, cex=.66)

## dev.copy(jpeg,'search6.jpg')
## dev.off()
