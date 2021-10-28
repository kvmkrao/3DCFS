c	shape function values

	tmp1 = sql(1,1)
	tmp2 = sql(1,2)
	tmp3 = sql(1,3)
	tmp4 = sql(1,4)
	tmp5 = sql(1,5)
	tmp6 = sql(1,6)
	tmp7 = sql(1,7)
	tmp8 = sql(1,8)

	xr11 =
     &	+ tmp1 * x(1, 1) + tmp2 * x(1, 2)
     &	+ tmp3 * x(1, 3) + tmp4 * x(1, 4)
     &	+ tmp5 * x(1, 5) + tmp6 * x(1, 6)
     &	+ tmp7 * x(1, 7) + tmp8 * x(1, 8)
	xr12 =
     &	+ tmp1 * x(2, 1) + tmp2 * x(2, 2)
     &	+ tmp3 * x(2, 3) + tmp4 * x(2, 4)
     &	+ tmp5 * x(2, 5) + tmp6 * x(2, 6)
     &	+ tmp7 * x(2, 7) + tmp8 * x(2, 8)
	xr13 =
     &	+ tmp1 * x(3, 1) + tmp2 * x(3, 2)
     &	+ tmp3 * x(3, 3) + tmp4 * x(3, 4)
     &	+ tmp5 * x(3, 5) + tmp6 * x(3, 6)
     &	+ tmp7 * x(3, 7) + tmp8 * x(3, 8)

	tmp1 = sql(2,1)
	tmp2 = sql(2,2)
	tmp3 = sql(2,3)
	tmp4 = sql(2,4)
	tmp5 = sql(2,5)
	tmp6 = sql(2,6)
	tmp7 = sql(2,7)
	tmp8 = sql(2,8)

	xr21 =
     &	+ tmp1 * x(1, 1) + tmp2 * x(1, 2)
     &	+ tmp3 * x(1, 3) + tmp4 * x(1, 4)
     &	+ tmp5 * x(1, 5) + tmp6 * x(1, 6)
     &	+ tmp7 * x(1, 7) + tmp8 * x(1, 8)
	xr22 =
     &	+ tmp1 * x(2, 1) + tmp2 * x(2, 2)
     &	+ tmp3 * x(2, 3) + tmp4 * x(2, 4)
     &	+ tmp5 * x(2, 5) + tmp6 * x(2, 6)
     &	+ tmp7 * x(2, 7) + tmp8 * x(2, 8)
	xr23 =
     &	+ tmp1 * x(3, 1) + tmp2 * x(3, 2)
     &	+ tmp3 * x(3, 3) + tmp4 * x(3, 4)
     &	+ tmp5 * x(3, 5) + tmp6 * x(3, 6)
     &	+ tmp7 * x(3, 7) + tmp8 * x(3, 8)

	tmp1 = sql(3,1)
	tmp2 = sql(3,2)
	tmp3 = sql(3,3)
	tmp4 = sql(3,4)
	tmp5 = sql(3,5)
	tmp6 = sql(3,6)
	tmp7 = sql(3,7)
	tmp8 = sql(3,8)

	xr31 =
     &	+ tmp1 * x(1, 1) + tmp2 * x(1, 2)
     &	+ tmp3 * x(1, 3) + tmp4 * x(1, 4)
     &	+ tmp5 * x(1, 5) + tmp6 * x(1, 6)
     &	+ tmp7 * x(1, 7) + tmp8 * x(1, 8)
	xr32 =
     &	+ tmp1 * x(2, 1) + tmp2 * x(2, 2)
     &	+ tmp3 * x(2, 3) + tmp4 * x(2, 4)
     &	+ tmp5 * x(2, 5) + tmp6 * x(2, 6)
     &	+ tmp7 * x(2, 7) + tmp8 * x(2, 8)
	xr33 =
     &	+ tmp1 * x(3, 1) + tmp2 * x(3, 2)
     &	+ tmp3 * x(3, 3) + tmp4 * x(3, 4)
     &	+ tmp5 * x(3, 5) + tmp6 * x(3, 6)
     &	+ tmp7 * x(3, 7) + tmp8 * x(3, 8)


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

	det = ( xr11 * cf11
     & 	      + xr21 * cf12
     &	      + xr31 * cf13 )

        cf11 = cf11/det
	cf12 = cf12/det
	cf13 = cf13/det
	cf21 = cf21/det
	cf22 = cf22/det
	cf23 = cf23/det
	cf31 = cf31/det
	cf32 = cf32/det
	cf33 = cf33/det

c	global first derivatives
      	sh(0,1)=sql(0,1)
	sh(1,1)=sql(1,1)*cf11+sql(2,1)*cf12+sql(3,1)*cf13
	sh(2,1)=sql(1,1)*cf21+sql(2,1)*cf22+sql(3,1)*cf23
	sh(3,1)=sql(1,1)*cf31+sql(2,1)*cf32+sql(3,1)*cf33
     	sh(0,2)=sql(0,2)
	sh(1,2)=sql(1,2)*cf11+sql(2,2)*cf12+sql(3,2)*cf13
	sh(2,2)=sql(1,2)*cf21+sql(2,2)*cf22+sql(3,2)*cf23
	sh(3,2)=sql(1,2)*cf31+sql(2,2)*cf32+sql(3,2)*cf33
      	sh(0,3)=sql(0,3)
	sh(1,3)=sql(1,3)*cf11+sql(2,3)*cf12+sql(3,3)*cf13
	sh(2,3)=sql(1,3)*cf21+sql(2,3)*cf22+sql(3,3)*cf23
	sh(3,3)=sql(1,3)*cf31+sql(2,3)*cf32+sql(3,3)*cf33
      	sh(0,4)=sql(0,4)
	sh(1,4)=sql(1,4)*cf11+sql(2,4)*cf12+sql(3,4)*cf13
	sh(2,4)=sql(1,4)*cf21+sql(2,4)*cf22+sql(3,4)*cf23
	sh(3,4)=sql(1,4)*cf31+sql(2,4)*cf32+sql(3,4)*cf33
     	sh(0,5)=sql(0,5)
	sh(1,5)=sql(1,5)*cf11+sql(2,5)*cf12+sql(3,5)*cf13
	sh(2,5)=sql(1,5)*cf21+sql(2,5)*cf22+sql(3,5)*cf23
	sh(3,5)=sql(1,5)*cf31+sql(2,5)*cf32+sql(3,5)*cf33
        sh(0,6)=sql(0,6)
	sh(1,6)=sql(1,6)*cf11+sql(2,6)*cf12+sql(3,6)*cf13
	sh(2,6)=sql(1,6)*cf21+sql(2,6)*cf22+sql(3,6)*cf23
	sh(3,6)=sql(1,6)*cf31+sql(2,6)*cf32+sql(3,6)*cf33
        sh(0,7)=sql(0,7)
	sh(1,7)=sql(1,7)*cf11+sql(2,7)*cf12+sql(3,7)*cf13
	sh(2,7)=sql(1,7)*cf21+sql(2,7)*cf22+sql(3,7)*cf23
	sh(3,7)=sql(1,7)*cf31+sql(2,7)*cf32+sql(3,7)*cf33
        sh(0,8)=sql(0,8)
	sh(1,8)=sql(1,8)*cf11+sql(2,8)*cf12+sql(3,8)*cf13
	sh(2,8)=sql(1,8)*cf21+sql(2,8)*cf22+sql(3,8)*cf23
	sh(3,8)=sql(1,8)*cf31+sql(2,8)*cf32+sql(3,8)*cf33

