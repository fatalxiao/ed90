\name{ufit}
\alias{ufit}
\title{
Unimodal isotonic regression.
}
\description{
A "divide and conquer" algorithm is applied to calculate the isotonic
regression of a set of data, for a unimodal order.  If the mode of
the unimodal order is not specified, then the optimal (in terms of
minimizing the error sum of squares) unimodal fit is calculated.
}
\usage{
ufit(y, lmode=NULL, imode=NULL, x=NULL, w=NULL, lc=TRUE, rc=TRUE,
        type=c("raw","stepfun","both"))
}
\arguments{
\item{y}{
Vector of data whose isotonic regression is to be calculated.
}
\item{lmode}{
Numeric scalar specifiing the location of the mode.  It must be one
of the values of \code{x} (see below) otherwise an error is thrown.
}
\item{imode}{
Integer scalar specifying the index, amongst the values of \code{x}
(see below) of the location of the mode.  It must be one of the
indices from 1 to \code{n}, where \code{n} is the length of \code{y},
otherwise an error is thrown.

It is an error to specify both \code{lmode} and \code{imode}.

Note that if neither \code{lmode} nor \code{imode} is specified
then the function performs an exhaustive search among all possible
mode locations for the optimal (in terms of minimising the error
sum of squares) location.
}
\item{x}{
A somewhat notional vector of \eqn{x} values corresponding to the
data vector \code{y}; the value of the mode must be given, or will
be determined in terms of these \eqn{x} values.  Conceptually the
model is \code{y = m(x) + E}, where \code{m()} is a unimodal function
with mode at \code{lmode}, and where \code{E} is random "error".  If \code{x} is
not specified, it defaults to an equi-spaced sequence of length \code{n}
on [0,1] (where \code{n} is the length of \code{y}).
}
\item{w}{
Optional vector of weights to be used for calculating a weighted
isotonic regression; if \code{w} is not specified then all weights
are taken to equal 1.
}
\item{lc}{
Logical scalar; should the isotonization be left continuous?  If
\code{lc==FALSE} then the value of the isotonization just before the
mode is set to \code{NA}, which causes line plots to have a jump
discontinuity at (just to the left of) the mode.  The default is
\code{lc=TRUE}.}
\item{rc}{
Logical scalar; should the isotonization be right continuous?  If
\code{rc==FALSE} then the value of the isotonization just after the
mode is set to \code{NA}, which causes line plots to have a jump
discontinuity at (just to the right of) the mode.  The default is
\code{rc=TRUE}.}
\item{type}{String specifying the type of the output; see \bold{Value}.
May be abbreviated.}
}
\value{
If \code{type=="raw"} then the value is
a list with components:
\item{x}{
The argument \code{x} if this is specified, otherwise the default value.
}
\item{y}{
The fitted values.
}
\item{mode}{
The value of the location of the mode as determined by \code{lmode}
or \code{imode} if one of these was specified.  Otherwise it is
the value of the location of the mode which was found to minimize
the error sum of squares.
}
\item{mse}{ The mean squared error.}

If \code{type=="both"} then a component \code{h} which is the step function
representation of the isotonic regression is added to the foregoing list.

If \code{type=="stepfun"} then only the step function representation
\code{h} is returned.
}
\details{
This function dynamically loads fortran subroutines "pava", "ufit"
and "unimode" to do the actual work.
}
\author{Rolf Turner
\email{rolfturner@posteo.net}
}
\references{
Mureika, R. A., Turner, T. R. and Wollan, P. C.  (1992).
An algorithm for unimodal isotonic regression, with
application to locating a maximum.  University of New Brunswick
Department of Mathematics and Statistics Technical Report Number
92 -- 4.

Robertson, T., Wright, F. T. and Dykstra, R.  L.  (1988).
Order Restricted Statistical Inference.  Wiley, New York.

Shi, Ning-Zhong. (1988) A test of homogeneity for
umbrella alternatives and tables of the level probabilities.
Commun.  Statist. --- Theory Meth. vol. 17, pp. 657 -- 670.

Turner, T. R., and Wollan, P. C. (1997) Locating a maximum using
isotonic regression.  Computational Statistics and Data Analysis
vol. 25, pp. 305 -- 320.
}
\seealso{\code{\link{pava}()} \code{\link{biviso}()} }

\examples{
y <- c(0,1,2,3,3,2)
f1 <- ufit(y,lmode=0.4) # The third entry of the default
                        # value of x = c(0.0,0.2,0.4,0.6,0.8,1.0).
f2 <- ufit(y,imode=3)   # Identical to f1.
f3 <- ufit(y,lmode=3,x=1:6)   # Effectively the same as f1 and f2.
                              # But is different  in appearance.
f4 <- ufit(y,imode=3,x=1:6)   # Identical to f3.

\dontrun{
    ufit(y,lmode=3)     # Throws an error.
    ufit(y,imode=7)     # Throws an error.
}

x <- c(0.00,0.34,0.67,1.00,1.34,1.67,2.00,2.50,3.00,3.50,4.00,4.50,
       5.00,5.50,6.00,8.00,12.00,16.00,24.00)
y <- c(0.0,61.9,183.3,173.7,250.6,238.1,292.6,293.8,268.0,285.9,258.8,
       297.4,217.3,226.4,170.1,74.2,59.8,4.1,6.1)
z <- ufit(y,x=x,type="b")
plot(x,y)
lines(z,col="red")
plot(z$h,do.points=FALSE,col.hor="blue",col.vert="blue",add=TRUE)
abline(v=z$mode,col="green",lty=2)
}
\keyword{regression}
\keyword{nonlinear}
