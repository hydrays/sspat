twotypes <- function(T, r1, r2, d1, d2, p, w){
  ## x is the system state
  x <- c(0, 0)

  ## define reaction system
  nr <- 5
  ns <- 2
  c <- c(r1, r2, d1, d2, p)
  nu <- matrix(0, ns, nr)
  nu[,1] = c(1, 0)
  nu[,2] = c(0, 1)
  nu[,3] = c(-1, 0)
  nu[,4] = c(0, -1)
  nu[,5] = c(-1,1)
  
  ## determing initial cell
  u <- runif(1)
  if ( u < w ) {
    x[1] = 1
  } else {
    x[2] = 1
  }

  ## Begin main loop
  .t <- 0
  while ( TRUE ){
    a <- getrate(x, c)
    if ( sum(a) == 0 ) break
    ## cat(a, '\n')
    ca <- cumsum(a)
    tau <- rexp(1, rate = 1)
    tau <- tau/ca[nr]
    .t <- .t + tau
    ## cat(.t, '\n')
    if ( .t > T ) {
      ## cat(.t, x, '\n')
      ## cat('breaked')
      break
    }
    
    u <- runif(1)
    u <- u*ca[nr]
    rindex <- -1
    for ( j in seq(nr) ){
      if ( u < ca[j] ) {
        rindex <- j
        break
      }
    }
    x <- x + nu[,rindex]
    ## cat(x, '\n')
  }
  return(x)
}

getrate <- function(x, c){
  a <- c(c[1]*x[1], c[2]*x[2], c[3]*x[1], c[4]*x[2], c[5]*x[1])
  return(a)
}
