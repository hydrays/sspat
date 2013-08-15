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
  Content[25] <- paste('\tmp2 = ', mp2, ',')
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
  ErrorA <- norm(as.matrix(ExpResultA$y-SimResultA$y))

  ## Newborn part
  datafile <- sprintf("out/CellNewborn.dat")
  z <- matrix(scan(datafile, n=NPool*L, quiet=TRUE),
              NPool, L, byrow=TRUE)
  d2 <- density(z[,1])
  SimResultB <- spline(d2, n=101, xmin=0, xmax=3000)
  SimResultB$y <- pmax(0, SimResultB$y)
  SimResultB$y <- SimResultB$y/(sum(SimResultB$y)*30)
  ExpResultB <- read.csv('newborn_dist.csv')
  ErrorB <- norm(as.matrix(ExpResultB$y-SimResultB$y))

                                        # Total error
  ErrorT <- ErrorA + ErrorB

  return(ErrorT)
}
