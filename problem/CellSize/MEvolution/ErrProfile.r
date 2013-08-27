source('MitoFun.r')
n <- 20
mp1 <- seq(0.1, 2, length.out=n)
mp2 <- seq(0.1, 2, length.out=n)
Err <- matrix(0, n, n)
for ( i in seq(n) ) {
    for ( j in seq(n) ) {
        Err[i, j] <- MitoFun(c(mp1[i], mp2[j]))
        #Err[i, j] <- mp1[i] + mp2[j]
    }
}
