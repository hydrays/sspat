\name{PMValidate}
\alias{PMValidate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{PMValidate
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PMValidate()
}
%- maybe also 'usage' for other objects documented here.
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
function(){

Nsample = 1000
T = 8
r1 = 0.8
r2 = 0.3
p1 = c(0.2, 0.2, 0.2, 0.2, 0.2)
p2 = c(0.8, 0.2)
w = 1.0 # means all good cell
x <- seq(Nsample)
for (i in seq(Nsample)){
  y <- simulator2(T, r1, r2, p1, p2, w)
  x[i] <- sum(y)
  }

data7 <- read.csv('data/8dayfast.csv')
mcell <- as.numeric(data7[8,3:ncol(data7[1,])])
mcell[is.na(mcell)] <- 0
mcell <- mcell[mcell>=50]

Fn <- ecdf(x)
Fe <- ecdf(mcell)

plot(Fn, xlim=c(0, 300))
lines(Fe, col='red')

dis <- ks.test(mcell, x)
cat(dis$p.value)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
