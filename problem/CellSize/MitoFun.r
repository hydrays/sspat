#######################################
## Error function for optimazation
#######################################

MitoFun <- function(mpv){

  if (length(mpv) != 2 ) {
    cat('Number of parameters must equal to p\n')
  }

  ## Prepare control file
  mp1 <- mpv[1]
  mp2 <- mpv[2]

  Content <- readLines('control.txt')
  Content[24] <- paste('\tmp1 = ', mp1, ',')
  Content[25] <- paste('\tmp2 = ', mp2)
  writeLines(Content, 'control.txt')
  
  cat('Invoking Fortran program\n')
  system('run > log')

  parainfo <- read.csv("out/control.csv", strip.white=TRUE)
  NPool <- parainfo$VALUE[parainfo$PARAMETER=='NPool']
  L = 6
  
  ## Asyn part
  i <- 490
  datafile <- sprintf("%s%05d%s", "out/m", i, ".dat")
  z <- matrix(scan(datafile, n=NPool*L, quiet=TRUE),
              NPool, L, byrow=TRUE)
  d1 <- density(z[,1])
  SimResultA <- spline(d1, n=101, xmin=0, xmax=3000)
  SimResultA$y <- pmax(0, SimResultA$y)
  SimResultA$y <- SimResultA$y/(sum(SimResultA$y)*30)
  ExpResultA <- read.csv('asyn_dist.csv')
  ## L1 norm
  ## ErrorA <- norm(as.matrix(ExpResultA$y-SimResultA$y))
  ## KL divgence
  tt1<-cbind(ExpResultA$x, ExpResultA$y)
  tt2<-cbind(SimResultA$x, SimResultA$y)
  ErrorA <- kl.dist(tt1,tt2)$D1
  
  ## Newborn part
  datafile <- sprintf("out/CellNewborn.dat")
  z <- matrix(scan(datafile, n=NPool*L, quiet=TRUE),
              NPool, L, byrow=TRUE)
  d2 <- density(z[,1])
  SimResultB <- spline(d2, n=101, xmin=0, xmax=3000)
  SimResultB$y <- pmax(0, SimResultB$y)
  SimResultB$y <- SimResultB$y/(sum(SimResultB$y)*30)
  ExpResultB <- read.csv('newborn_dist.csv')
  ## L1 norm
  ## ErrorB <- norm(as.matrix(ExpResultB$y-SimResultB$y))
  ## KL divgence
  tt1<-cbind(ExpResultB$x, ExpResultB$y)
  tt2<-cbind(SimResultB$x, SimResultB$y)
  ErrorB <- kl.dist(tt1,tt2)$D1
  
  ## Total error
  ErrorT <- ErrorA + ErrorB

  ##cat('Para:', mp1, mp2, '-> Error:', ErrorT, '\n')
  print(c(mp1, mp2, ErrorT), digits=16)
  return(ErrorT)
}
