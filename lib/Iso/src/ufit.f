C Output from Public domain Ratfor, version 1.03
      subroutine ufit(y,w,imode,ymdf,wmdf,mse,y1,w1,y2,w2,ind,kt,n)
      implicit double precision(a-h,o-z)
      double precision imode, mse, imax
      dimension y(n), w(n), ymdf(n), wmdf(n),y1(n), w1(n), y2(n), w2(n),
     * ind(n), kt(n)
      if(imode .lt. 0)then
      m = n-1
      tau = 1.5d0
      imax = -1.d0
      ssemin = 1.d200
      do23002 i = 1,m 
      do23004 j = 1,n 
      ymdf(j) = y(j)
      wmdf(j) = w(j)
23004 continue
23005 continue
      call unimode(ymdf,wmdf,y1,w1,y2,w2,ind,kt,tau,n)
      sse = 0.d0
      do23006 j = 1,n 
      sse = sse + (ymdf(j)-y(j))**2
23006 continue
23007 continue
      if(sse .lt. ssemin)then
      ssemin = sse
      imax = tau
      endif
      tau = tau+1.d0
23002 continue
23003 continue
      k1 = int(imax-0.5d0)
      k2 = int(imax+0.5d0)
      else
      imax = imode
      endif
      do23010 j = 1,n 
      ymdf(j) = y(j)
      wmdf(j) = w(j)
23010 continue
23011 continue
      call unimode(ymdf,wmdf,y1,w1,y2,w2,ind,kt,imax,n)
      if(imode .lt. 0)then
      mse = ssemin/dble(n)
      if(ymdf(k1).ge.ymdf(k2))then
      imode=dble(k1)
      else
      imode=dble(k2)
      endif
      else
      sse = 0.d0
      do23016 j = 1,n 
      sse = sse + (ymdf(j)-y(j))**2
23016 continue
23017 continue
      mse = sse/dble(n)
      endif
      return
      end
