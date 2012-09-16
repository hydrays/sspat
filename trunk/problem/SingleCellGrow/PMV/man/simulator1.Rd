\name{simulator1}
\alias{simulator1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simulator1
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simulator1(T, r1, r2, d1, d2, p, w)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T}{
%%     ~~Describe \code{T} here~~
}
  \item{r1}{
%%     ~~Describe \code{r1} here~~
}
  \item{r2}{
%%     ~~Describe \code{r2} here~~
}
  \item{d1}{
%%     ~~Describe \code{d1} here~~
}
  \item{d2}{
%%     ~~Describe \code{d2} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{w}{
%%     ~~Describe \code{w} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(T, r1, r2, d1, d2, p, w){
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
