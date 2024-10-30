subroutine pava(y,w,kt,n)
implicit double precision(a-h,o-z)
logical same
dimension y(n), w(n), kt(n)

# Note: `kt' <--> `keep track' (of the level sets).

do i = 1,n {
	kt(i) = i
}

if(n==1) return

repeat{
	same = .true.
	do i = 2,n {
		if(y(i-1) > y(i)) {
			k1 = kt(i)
			k2 = kt(i-1)
			do j = 1,n {
				if(kt(j)==k1) kt(j) = k2
			}
			wnew = w(i-1) + w(i)
			ynew = (w(i-1)*y(i-1)+w(i)*y(i))/wnew
			do j = 1,n {
				if(kt(j)==k2) {
					y(j) = ynew
					w(j) = wnew
				}
			}
			same = .false.
		}
	}
	if(same) break
}

return
end
