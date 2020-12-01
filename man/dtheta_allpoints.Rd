\name{dtheta_allpoints}
\alias{dtheta_allpoints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compute the local dimension and the persistence for several points.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
   Compute the local dimension and the persistence for several points using the future API for parallel computing.
}
\usage{
dtheta_allpoints(X0, Xref, rho = 0.98, method = "Ferro", method.args = list(), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X0}{
%%     ~~Describe \code{x} here~~
  The N points for which the local dimension and persistence are to be estimated.
  It is a matrix of dimension NxD where D is the number of variable defining the state of the system.
}
  \item{Xref}{
%%     ~~Describe \code{x} here~~
  A set of points defining the trajectory of the dynamical system over T timesteps.
  It is a matrix of dimension TxD where D is the number of variable defining the state of the system.
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
  \item{...}{
%%     ~~Describe \code{x} here~~
 A named list with additionnal arguments for the parallel API future.
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
  return a matrix of of dimension Nx2 with the estimates of local dimension in the first colunm and of the persistence in the second.
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
future::plan("multisession")
dtheta_allpoints(X, X, rho = 0.98, method = "Caby", method.args = list(m = 1))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

