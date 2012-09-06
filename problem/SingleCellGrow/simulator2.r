simulator2 <- function(T, r1, r2, p1, p2, w){
  ## x is the system state
  x <- c(0, 0)

  if ( length(p1)!=5 ) {
    cat('p1 is not set correctly. p1 = [p12, p11, p10, p1d, p1a] where \n')
    cat('p12 ===> Type1 -> Type1 + Type1, \n')
    cat('p11 ===> Type1 -> Type1 + Type2, \n')
    cat('p10 ===> Type1 -> Type2 + Type2, \n')
    cat('p1d ===> Type1 -> Type2, \n')
    cat('p1a ===> Type1 -> 0, \n')
    return(-1)
  }


  if ( length(p2)!=2 ) {
    cat('p2 is not set correctly. p2 = [p22, p21] where \n')
    cat('p22 ===> Type2 -> Type2 + Type2, \n')
    cat('p2a ===> Type2 -> 0, \n')
    return(-1)
  }
  
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
