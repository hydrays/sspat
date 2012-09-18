Csimulator1 <- function(T, l1, l2, d1, d2, v, w){
  if (FALSE){
    stop("argument wrong!")
  }
  out <- .C("simulator1", T = as.double(T), l1 = as.double(l1),
            l2 = as.double(l2), d1 = as.double(d1),
            d2 = as.double(d2), v = as.double(v),
            w = as.double(w), x = integer(2))
#  cat(summary(out))
  return(out$x)
}
