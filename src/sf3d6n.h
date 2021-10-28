c   Wedge  shape function values
	
	tmp1 = sq3(1,1,ig)
	tmp2 = sq3(1,2,ig)
	tmp3 = sq3(1,3,ig)
	tmp4 = sq3(1,4,ig)
	tmp5 = sq3(1,5,ig)
	tmp6 = sq3(1,6,ig)


	xr11 =
     &	+ tmp1 * x1(1, 1) + tmp2 * x1(1, 2)
     &	+ tmp3 * x1(1, 3) + tmp4 * x1(1, 4)
     &	+ tmp5 * x1(1, 5) + tmp6 * x1(1, 6)

	xr12 =
     &	+ tmp1 * x1(2, 1) + tmp2 * x1(2, 2)
     &	+ tmp3 * x1(2, 3) + tmp4 * x1(2, 4)
     &	+ tmp5 * x1(2, 5) + tmp6 * x1(2, 6)

	xr13 =
     &	+ tmp1 * x1(3, 1) + tmp2 * x1(3, 2)
     &	+ tmp3 * x1(3, 3) + tmp4 * x1(3, 4)
     &	+ tmp5 * x1(3, 5) + tmp6 * x1(3, 6)


	tmp1 = sq3(2,1,ig)
	tmp2 = sq3(2,2,ig)
	tmp3 = sq3(2,3,ig)
	tmp4 = sq3(2,4,ig)
	tmp5 = sq3(2,5,ig)
	tmp6 = sq3(2,6,ig)


	xr21 =
     &	+ tmp1 * x1(1, 1) + tmp2 * x1(1, 2)
     &	+ tmp3 * x1(1, 3) + tmp4 * x1(1, 4)
     &	+ tmp5 * x1(1, 5) + tmp6 * x1(1, 6)

	xr22 =
     &	+ tmp1 * x1(2, 1) + tmp2 * x1(2, 2)
     &	+ tmp3 * x1(2, 3) + tmp4 * x1(2, 4)
     &	+ tmp5 * x1(2, 5) + tmp6 * x1(2, 6)

	xr23 =
     &	+ tmp1 * x1(3, 1) + tmp2 * x1(3, 2)
     &	+ tmp3 * x1(3, 3) + tmp4 * x1(3, 4)
     &	+ tmp5 * x1(3, 5) + tmp6 * x1(3, 6)


	tmp1 = sq3(3,1,ig)
	tmp2 = sq3(3,2,ig)
	tmp3 = sq3(3,3,ig)
	tmp4 = sq3(3,4,ig)
	tmp5 = sq3(3,5,ig)
	tmp6 = sq3(3,6,ig)


	xr31 =
     &	+ tmp1 * x1(1, 1) + tmp2 * x1(1, 2)
     &	+ tmp3 * x1(1, 3) + tmp4 * x1(1, 4)
     &	+ tmp5 * x1(1, 5) + tmp6 * x1(1, 6)

	xr32 =
     &	+ tmp1 * x1(2, 1) + tmp2 * x1(2, 2)
     &	+ tmp3 * x1(2, 3) + tmp4 * x1(2, 4)
     &	+ tmp5 * x1(2, 5) + tmp6 * x1(2, 6)

	xr33 =
     &	+ tmp1 * x1(3, 1) + tmp2 * x1(3, 2)
     &	+ tmp3 * x1(3, 3) + tmp4 * x1(3, 4)
     &	+ tmp5 * x1(3, 5) + tmp6 * x1(3, 6)



c	jacobian

	cf11 = + (xr22*xr33 - xr32*xr23)
	cf12 = - (xr12*xr33 - xr32*xr13)
	cf13 = + (xr12*xr23 - xr22*xr13)
	cf21 = - (xr21*xr33 - xr31*xr23)
	cf22 = + (xr11*xr33 - xr31*xr13)
	cf23 = - (xr11*xr23 - xr21*xr13)
	cf31 = + (xr21*xr32 - xr31*xr22)
	cf32 = - (xr11*xr32 - xr31*xr12)
	cf33 = + (xr11*xr22 - xr21*xr12)

c	write(6,180)xr21,xr33,cf21
c	write(6,180)xr11,xr33, cf22
c	stop
c	write(65,180)cf22,xr11,xr33
	det = ( xr11 * cf11
     & 	      + xr21 * cf12
     &	      + xr31 * cf13 )
c	write(6,180)cf21,cf22,det
c	stop
        cf11 = cf11/det
	cf12 = cf12/det
	cf13 = cf13/det
	cf21 = cf21/det
	cf22 = cf22/det
	cf23 = cf23/det
	cf31 = cf31/det
	cf32 = cf32/det
	cf33 = cf33/det
cc	write(6,180)cf21,cf22
c	stop

c	global first derivatives
      	sh(0,1)=sq3(0,1,ig)
	sh(1,1)=sq3(1,1,ig)*cf11+sq3(2,1,ig)*cf12+sq3(3,1,ig)*cf13
	sh(2,1)=sq3(1,1,ig)*cf21+sq3(2,1,ig)*cf22+sq3(3,1,ig)*cf23
	sh(3,1)=sq3(1,1,ig)*cf31+sq3(2,1,ig)*cf32+sq3(3,1,ig)*cf33
     	sh(0,2)=sq3(0,2,ig)
	sh(1,2)=sq3(1,2,ig)*cf11+sq3(2,2,ig)*cf12+sq3(3,2,ig)*cf13
	sh(2,2)=sq3(1,2,ig)*cf21+sq3(2,2,ig)*cf22+sq3(3,2,ig)*cf23
	sh(3,2)=sq3(1,2,ig)*cf31+sq3(2,2,ig)*cf32+sq3(3,2,ig)*cf33
      	sh(0,3)=sq3(0,3,ig)
	sh(1,3)=sq3(1,3,ig)*cf11+sq3(2,3,ig)*cf12+sq3(3,3,ig)*cf13
	sh(2,3)=sq3(1,3,ig)*cf21+sq3(2,3,ig)*cf22+sq3(3,3,ig)*cf23
	sh(3,3)=sq3(1,3,ig)*cf31+sq3(2,3,ig)*cf32+sq3(3,3,ig)*cf33
      	sh(0,4)=sq3(0,4,ig)
	sh(1,4)=sq3(1,4,ig)*cf11+sq3(2,4,ig)*cf12+sq3(3,4,ig)*cf13
	sh(2,4)=sq3(1,4,ig)*cf21+sq3(2,4,ig)*cf22+sq3(3,4,ig)*cf23
	sh(3,4)=sq3(1,4,ig)*cf31+sq3(2,4,ig)*cf32+sq3(3,4,ig)*cf33
     	sh(0,5)=sq3(0,5,ig)
	sh(1,5)=sq3(1,5,ig)*cf11+sq3(2,5,ig)*cf12+sq3(3,5,ig)*cf13
	sh(2,5)=sq3(1,5,ig)*cf21+sq3(2,5,ig)*cf22+sq3(3,5,ig)*cf23
	sh(3,5)=sq3(1,5,ig)*cf31+sq3(2,5,ig)*cf32+sq3(3,5,ig)*cf33
        sh(0,6)=sq3(0,6,ig)
	sh(1,6)=sq3(1,6,ig)*cf11+sq3(2,6,ig)*cf12+sq3(3,6,ig)*cf13
	sh(2,6)=sq3(1,6,ig)*cf21+sq3(2,6,ig)*cf22+sq3(3,6,ig)*cf23
	sh(3,6)=sq3(1,6,ig)*cf31+sq3(2,6,ig)*cf32+sq3(3,6,ig)*cf33

