for (j in seq(7)){
  N <- 100
  m <- c(0, 50, 100, 200, 800, 1000, 10000)
  T <- 14
  x <- c(100, m[j])
  r <- seq(T+1)
  r[1] = x[1]/sum(x)
  for (i in seq(T)){
    tempx <- c(0, 0)
    for (k in seq(x[1])){
      y <- simulator1(1, 1, 0.6, 0.4, 0.3, 0.3, 1)
      tempx <- tempx + y
    }
    for (k in seq(x[2])){
      y <- simulator1(1, 1, 0.6, 0.4, 0.3, 0.3, 0)
      tempx <- tempx + y
    }
    x <- tempx
    cat(i, x, r[i], '\n')
    r[i+1] <- x[1]/sum(x) 
  }

  if (j==1) {
    plot(0:T, r, type="o", pch=1, main="ratio vs time",
         xlim = c(0, 14),
         ylim = c(0, 1),
         xlab="time (day)", ylab="ratio (good/total)")
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
