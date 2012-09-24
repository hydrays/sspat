simulator1p <- function(T, l1, l2, d1, d2, v, x){
  if (!is.numeric(T)){
    stop("argument wrong!")
  }
  out <- .C("simulator1p", T = as.double(T), l1 = as.double(l1),
            l2 = as.double(l2), d1 = as.double(d1),
            d2 = as.double(d2), v = as.double(v),
            x = as.integer(x))
#  cat(summary(out))
  return(out$x)
}
