subroutine ufit(y,w,imode,ymdf,wmdf,mse,y1,w1,y2,w2,ind,kt,n)
implicit double precision(a-h,o-z)
double precision imode, mse, imax
dimension y(n), w(n), ymdf(n), wmdf(n),y1(n), w1(n), y2(n), w2(n), ind(n), kt(n)

# Nude virgin of ufit --- 18/8/95.
# The changes are based upon Pete's observation that when we are
# seeking the OPTIMIUM mode (in terms of SSE) we need only search
# over the ``half-points'' --- 1.5, 2.5, ..., n-0.5.  If the optimum
# is at k, then the half-points k-0.5 and k+05 give SSEs that are
# at least as small as and hence are equal to the SSE at k.  This is
# because if a function is increasing on 1,...,k and decreasing on
# k,...n, then it is increasing on 1,...,(k-1)) and decreasing on
# k,...,n !!!  (And likewise for 1,...,k and (k+1),...n.)   Thus if
# there is an optimum at k then there are optima at k-0.5 and k+0.5.
# Of course if k=1 then k-0.5 is not considered and likewise if k=n
# then k+0.5 is not considered.
#
# Note also that if there is an optimum at the half-point k-0.5, then
# there is also a whole-point optimum either at k-1 or k.
#
# Explanation revised (corrected) 31/05/2015.
#
# Note that in this subroutine there is no "x" argument.  The
# *indices* of a conceptual "x" are dealt with. (22/07/2023)

if(imode < 0) {
	m      =  n-1
	tau    =  1.5d0
	imax   = -1.d0
	ssemin =  1.d200

	do i = 1,m {
		do j = 1,n {
			ymdf(j) = y(j)
			wmdf(j) = w(j)
		}
		call unimode(ymdf,wmdf,y1,w1,y2,w2,ind,kt,tau,n)
		sse = 0.d0
		do j = 1,n {
			sse = sse + (ymdf(j)-y(j))**2
		}
		if(sse < ssemin) {
			ssemin = sse
			imax   = tau
		}
		tau = tau+1.d0
	}
	k1 = int(imax-0.5d0)
	k2 = int(imax+0.5d0)
}
else imax = imode
do j = 1,n {
   ymdf(j) = y(j)
   wmdf(j) = w(j)
}
#call dblepr("imax:",-1,imax,1)
#call dblepr("y:",-1,y,6)
#call dblepr("w:",-1,w,6)
#call dblepr("ymdf:",-1,ymdf,6)
#call dblepr("wmdf:",-1,wmdf,6)

call unimode(ymdf,wmdf,y1,w1,y2,w2,ind,kt,imax,n)
#call labpepr("Got past first call to unimode.",-1)

if(imode < 0) {
	mse = ssemin/dble(n)
	if(ymdf(k1)>=ymdf(k2)) imode=dble(k1)
	else imode=dble(k2)
}
else {
	sse = 0.d0
	do j = 1,n {
		sse = sse + (ymdf(j)-y(j))**2
	}
	mse = sse/dble(n)
}
return
end
