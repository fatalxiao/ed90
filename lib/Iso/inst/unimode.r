subroutine unimode(y,w,y1,w1,y2,w2,ind,kt,tau,n)
implicit double precision(a-h,o-z)
dimension y(n), w(n), y1(n), w1(n), y2(n), w2(n), ind(n), kt(n)

# Handle the linear ordering cases:
if(tau >= dble(n)) {
	call pava(y,w,kt,n)
	return
}

if(tau <= 1.d0) {
	do i = 1,n {
		j = n+1-i
		y2(i) = y(j)
                w2(i) = w(j)
        }
	call pava(y2,w2,kt,n)
	do i = 1,n {
		j = n+1-i
		y(i) = y2(j)
		w(i) = w2(j)
	}
	return
}

k1 = 0
k2 = 0
do i = 1,n {
	if(i < tau) {
		y1(i) = y(i)
		w1(i) = w(i)
		k1 = k1+1
	}
	if(i > tau) {
		j = n+1-i
		y2(j) = y(i)
		w2(j) = w(i)
		k2 = k2+1
	}
}

if(k1==0) {
    call rexit("The index of the mode is 0.\n")
}
if(k2==0) {
    call rexit("The index of the mode is one more than the number of indices.\n")
}

if(k1+k2 == n) {
	call pava(y1,w1,kt,k1)
	do i = 1,k1 {
		y(i) = y1(i)
		w(i) = w1(i)
	}
	call pava(y2,w2,kt,k2)
	do i = 1,k2 {
		j = n+1-i
		y(j) = y2(i)
		w(j) = w2(i)
	}
	return
}

if(k1+k2 == n-1) {
	yk = y(k1+1)
	call pava(y1,w1,kt,k1)
	call pava(y2,w2,kt,k2)
	i1 = 1
	i2 = 1
	i  = 1
	repeat{
		if(i1 <= k1) t1 = y1(i1)
		else t1 = y2(k2)+1.d10
		if(i2 <= k2) t2 = y2(i2)
		else t2 = y1(k1)+1.d10

		if(t1 < t2) {
			y(i) = y1(i1)
			ind(i) = i1
			i1 = i1+1
		}
		else {
			y(i) = y2(i2)
			ind(i) = n-i2+1
			i2 = i2+1
		}
		i = i + 1
		if(i == n) break
	}
	y(n) = yk
	ind(n) = k1+1
	do i = 1,n {
		w1(ind(i)) = w(i)
	}
	do i = 1,n {
		w(i) = w1(i)
	}
	call pava(y,w,kt,n)
	do i = 1,n {
		y1(ind(i)) = y(i)
		w1(ind(i)) = w(i)
	}
	do i = 1,n {
		y(i) = y1(i)
		w(i) = w1(i)
	}
} else {
    call rexit("The total length of the monotone segments is neither n nor n-1.")
}
return
end
