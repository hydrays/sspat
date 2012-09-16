getrate <-
function(x, c){
  a <- c(c[1]*x[1], c[2]*x[1], c[3]*x[1], c[4]*x[1], c[5]*x[1],
         c[6]*x[2], c[7]*x[2])
  return(a)
}

