dyn.load('simulator1.so')
source('simulator1s.r')

T <- 8
r1 <- 1
r2 <- 0.6
d1 <- 0.4
d2 <- 0.3
v <- 0.3

N <- 100000
w <- double(N)
n <- 0
for ( i in seq(N) ){
  x <- simulator1s(T, r1, r2, d1, d2, v, 1)
  ## cat(x,'\n')
  if ( sum(x) > 100 ){
    n <- n+1    
    w[n] <- x[1]/sum(x)
  }
}

den <- density(w[1:n])
plot(den)
