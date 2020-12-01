\name{Sueveges}
\alias{Sueveges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Sueveges's Estimator of the Extremal Index
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
 estimates the extremal index with the  Sueveges's Estimator
}
\usage{
Sueveges(X, rho)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{x} here~~
  a vector of numerical value.
}
  \item{rho}{
%%     ~~Describe \code{x} here~~
  the empirical quantile level used to define the time series of excessess over threshold.
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
  return the estimate of the extremal index
}
\references{
%% ~put references to the literature/web site here ~
Mária Süveges (2007). Likelihood estimation of the extremal index. Extremes. 10:41–55 DOI 10.1007/s10687-007-0034-2

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
X <- arima.sim(list(order=c(1,0,0), ar=ar1), n=n)*sqrt(1-ar1^2)
Sueveges(X, q)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

