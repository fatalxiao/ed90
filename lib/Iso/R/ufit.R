ufit <- function(y,lmode=NULL,imode=NULL,x=NULL,w=NULL,lc=TRUE, rc=TRUE,
                 type=c("raw","stepfun","both")) {
#
# Function `ufit'.  Calculates the isotonic unimodal fit to a data
# sequence y, with mode at ``lmode''.  If lmode==NULL, then ufit() determines
# the optimal (least squares) location for the mode, and the fit
# for this optimum value.  The optimum mode may, by virtue of the
# nature of the procedure, be taken to be one of the points x_i, i =
# 1, ..., n.  NOTE that the optimum will occur at one of the
# midpoints (x_i + x_{i+1})/2, i = 1, ..., n-1 AND at one of the two
# adjacent points, i.e. either at x_i or at x_{i+1}.  If x is null, x
# is taken to be an equispaced sequence on [0,1].
#
type <- match.arg(type)

n <- length(y)
if(is.null(w)) w <- rep(1,n)
if(is.null(x)) x <- seq(0,1,length=n)
nargtype <- 1 + (!is.null(lmode)) + (!is.null(imode))*2
# 1 <--> neither is specified.
# 2 <--> lmode is specified, imode is not
# 3 <--> imode is specified, lmode is not
# 4 <--> both are specified

switch(EXPR=nargtype,
    {
# 1; neither lmode nor imode has been specified.
        imode <- -1 # Triggers search for optimum.
    },
    {
# 2; lmode is specified, imode is not specified.
        if(!lmode %in% x) {
            whinge <- paste0("If \"lmode\" is specified, it must be an entry\n",
                             "  of \"x\" which defaults to",
                             " seq(0,1,length=length(y)).\n")
            stop(whinge)
        }
        imode <- which(x==lmode)
    },
# 3; lmode is not specified, imode is specified.
    {
        if(!imode %in% 1:n) {
            whinge <- paste0("If \"imode\" is specified, it must be an integer\n",
                             "  between 1 and n = length(y).\n")
            stop(whinge)
        }
        lmode <- x[imode]
    },
# 4; both lmode and imode are specified --- error.
    {
        stop("Specify at most one of \"lmode\" and \"imode\".\n")
    }
)

rslt <- .Fortran(
		"ufit",
		y=as.double(y),
		w=as.double(w),
		imode=as.double(imode),
		ymdf=double(n),
		wmdf=double(n),
		mse=double(1),
		y1=double(n),
		w1=double(n),
		y2=double(n),
		w2=double(n),
		ind=integer(n),
		kt=integer(n),
		n=as.integer(n),
		PACKAGE="Iso"
		)
if(nargtype==1) {
    imode <- rslt$imode
    lmode <- x[imode]
}
ystar <- rslt$ymdf
if(type%in%c("stepfun","both")) {
	kind <- 1+which(diff(ystar)!=0)
	if(!(n%in%kind)) kind <- c(kind,n)
	y0   <- c(ystar[1],ystar[kind])
	h    <- stepfun(x[kind],y0)
}
i     <- floor(imode)
if(!lc) ystar[i]   <- NA
if( (!rc) & (i < n) ) ystar[i+1] <- NA
switch(type,raw=list(x=x,y=ystar,mode=lmode,mse=rslt$mse),
            stepfun=h,
            both=list(x=x,y=ystar,mode=lmode,mse=rslt$mse,h=h))
}
