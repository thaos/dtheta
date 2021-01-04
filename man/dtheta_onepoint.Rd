\name{dtheta_onepoint}
\alias{dtheta_onepoint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compute the local dimension and the persistence of one point.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~

}
\usage{
dtheta_onepoint(dX, rho = 0.98, method = "Ferro", method.args = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dX}{
%%     ~~Describe \code{x} here~~
  A vector of distances between the point of interest and the others.
  }
  \item{rho}{
%%     ~~Describe \code{x} here~~
  The empirical quantile level used to compute the local dimension and the persitance indicators.
  }
  \item{method}{
%%     ~~Describe \code{x} here~~
  A string with name of the persistence estimator to use. Availble choices are "Ferro", "Sueveges" and "Caby"
  }
  \item{method.args}{
%%     ~~Describe \code{x} here~~
 A named list with additionnal arguments needed by the persistence estimator. The argument "m" for the Caby estimator.
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
  return a vector of two values: the estimates of local dimension and of the persistence.
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
n <- 1000
q <- 0.99
ar1 <- 0.5
X <- as.vector(arima.sim(list(order=c(1,0,0), ar=ar1), n=n)*sqrt(1-ar1^2))
dX <- -log(sqrt((X - X[1])^2))
dtheta_onepoint(dX, rho = 0.98, method = "Caby", method.args = list(m = 1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

