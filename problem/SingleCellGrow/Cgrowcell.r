dyn.load('simu1.so')
source('Csimulator1.r')

for (j in seq(7)){
  N <- 100
  m <- c(0, 50, 100, 200, 800, 2000, 10000)
  T <- 30
  x <- c(100, m[j])
  r <- seq(T+1)
  r[1] <- x[1]/sum(x)

  #model parameters
  l1 <- 1
  l2 <- 0.6
  d1 <- 0.3
  d2 <- 0.3
  u <- 0.3

  for (i in seq(T)){
    x <- Csimulator1(1, l1, l2, d1, d2, u, x)
    cat(i, x, r[i], '\n')
    r[i+1] <- x[1]/sum(x) 
  }

  if (j==1) {
    plot(0:T, r, type="o", pch=1, main="ratio vs time",
         xlim = c(0, T),
         ylim = c(0, 1),
         xlab="time (day)", ylab="ratio (good/total)")
    abline(h=c(((l1-d1)-(l2-d2)-u)/((l1-d1)-(l2-d2))))
  }
  else{
    lines(0:T, r, type="o", col=j, pch=j)
  }
}

text(8, 0.8, adj=0, "Parameters")
text(8, 0.75, adj=0, "r1=1, d1=0.4")
text(8, 0.7, adj=0, "r2=0.6, d2=0.3")
text(8, 0.65, adj=0, "p=0.3")
text(8, 0.6, adj=0, "w is changing")

dev.copy(jpeg, "T40.jpg")
dev.off()
