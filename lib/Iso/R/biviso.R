biviso <- function(y, w=NULL,eps=NULL,eps2=1e-9,ncycle=50000,
                   fatal=TRUE,warn=TRUE) {
#
# Function 'biviso'.  To perform bivariate isotonic regression for simple
# (increasing) linear ordering on both variables.  Uses Applied Statistics
# Algorithm AS 206 (Isotonic regression in two independent variables;
# Gordon Bril, Richard Dykstra, Carolyn Pillers, and Tim Robertson;
# Algorithm AS 206; JRSSC (Applied Statistics), vol. 33, no. 3, pp.
# 352-357, 1984.)

# Check that ncycle makes sense:
if(ncycle!=round(ncycle) | ncycle < 2)
    stop("Argument ncycle must be an integer with value at least 2.\n")

# Check that y is of the right shape:
if(!is.numeric(y) | !is.matrix(y))
	stop("Argument \"y\" must be a numeric matrix.\n")

if(is.null(w)) w <- matrix(1,nrow=nrow(y),ncol=ncol(y)) else {
	if(!isTRUE(all.equal(dim(y),dim(w))))
		stop("Arguments \"y\" and \"w\" must have the same dimension.\n")
}

# Set epsilon:
if(is.null(eps)) eps <- sqrt(.Machine$double.eps)

nr   <- nrow(y)
nc   <- ncol(y)
nd   <- max(nr,nc)
rslt <- .Fortran(
	"smooth",
	NROW=as.integer(nr),
	NCOL=as.integer(nc),
	NDIM=as.integer(nd),
	X=as.double(y),
	W=as.double(w),
	A=double(4*nr*nc),
	B=double(5*nd),
	NCYCLE=as.integer(ncycle),
	ICYCLE=integer(1),
	G=double(nr*nc),
	EPS1=as.double(eps),
	EPS2=as.double(eps2),
	IFAULT=integer(1),
        FX=double(nd),
        PW=double(nd),
        W1=double(nd),
        WT=double(nd),
        NW=integer(nd),
	PACKAGE="Iso"
)
if(rslt$IFAULT != 0) {
	if(rslt$ifault == 4 && warn) {
		warning(paste("A near zero weight less than delta=0.00001\n",
                              "was replaced by delta.\n",sep=""))
	} else if(fatal) {
		stop(paste("Failed with ifault = ",rslt$ifault,".\n",sep=""))
        } else if(warn) {
		warning(paste("Algorithm gave ifault = ",rslt$ifault,".\n",sep=""))
	}
}
m <- matrix(rslt$G,nrow=nr,ncol=nc)
attr(m,"icycle") <- rslt$ICYCLE
attr(m,"ifault") <- rslt$IFAULT
m
}
