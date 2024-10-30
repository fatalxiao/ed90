\name{pava}
\alias{pava}
\alias{pava.sa}
\title{
Linear order isotonic regression.
}
\description{
The ``pool adjacent violators algorithm'' (PAVA) is applied to
calculate the isotonic regression of a set of data, with respect
to the usual increasing (or decreasing) linear ordering on the indices.
}
\usage{
pava(y, w, decreasing=FALSE, long.out=FALSE, stepfun=FALSE)
pava.sa(y, w, decreasing=FALSE, long.out=FALSE, stepfun=FALSE)
}
\arguments{
\item{y}{
Vector of data whose isotonic regression is to be calculated.
}
\item{w}{
Optional vector of weights to be used for calculating a
weighted isotonic regression; if w is not given, all weights
are taken to equal 1.
}
\item{decreasing}{Logical scalar; should the isotonic regression
be calculated with respect to \emph{decreasing} (rather than
increasing) order?}
\item{long.out}{
Logical argument controlling the nature of the value returned.
}
\item{stepfun}{
Logical scalar; if \code{TRUE} a step function representation
of the isotonic regression is returned.
}
}
\value{
If long.out is TRUE then the result returned consists of a list whose
components are:

\item{y}{
the fitted values
}
\item{w}{
the final weights
}
\item{tr}{
a set of indices made up of the smallest index in each level set,
which thus "keeps track" of the level sets.
}
\item{h}{
a step function which represents the results of the isotonic
regression.  This component is present \emph{only if} 
\code{stepfun} is \code{TRUE}.
}

If \code{long.out} is \code{FALSE} and \code{stepfun} is \code{TRUE}
then only the step function is returned.

If \code{long.out} and \code{stepfun} are both \code{FALSE} then only
the vector of fitted values is returned.

}
\details{
The function \code{pava()} uses dynamically loading of a fortran
subroutine "pava" to effect the computations.  The function \code{pava.sa()}
("sa" for "stand-alone") does all of the computations in raw R.  Thus
\code{pava.sa()} could be considerably slower for large data sets.

The \code{x} values for the step function returned by these functions
(if \code{stepfun} is \code{TRUE}) are thought of as being 1, 2,
\dots, \code{n=length(y)}.  The knots of the step function are the
\code{x} values (indices) \emph{following} changes in the \code{y}
values (i.e. the starting indices of the level sets, except for the
first level set).  The \code{y} value corresponding to the first level
set is the ``left hand'' value of \code{y} or \code{yleft}.
The step function is formed using the default arguments of \code{stepfun()}.
In particular it is \emph{right} continuous.
}
\author{Rolf Turner
\email{rolfturner@posteo.net}
}
\references{
Robertson, T., Wright, F. T. and Dykstra, R.  L.  (1988).
Order Restricted Statistical Inference.  Wiley, New York.
}
\seealso{\code{\link{ufit}()} \code{\link{stepfun}()} \code{\link{biviso}()} }
\examples{
# Increasing order:
y <- (1:20) + rnorm(20)
ystar <- pava(y)
plot(y)
lines(ystar,type='s')
# Decreasing order:
z <- NULL
for(i in 4:8) {
	z <- c(z,rep(8-i+1,i)+0.05*(0:(i-1)))
}
zstar <- pava(z,decreasing=TRUE)
plot(z)
lines(zstar,type='s')
# Using the stepfunction:
zstar <- pava(z,decreasing=TRUE,stepfun=TRUE)
plot(z)
plot(zstar,add=TRUE,verticals=FALSE,pch=20,col.points="red")
}
\keyword{regression}
\keyword{nonlinear}
