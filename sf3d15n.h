c   Wedge  shape function values
	
	tmp1 = sq3(1,1,ig)
	tmp2 = sq3(1,2,ig)
	tmp3 = sq3(1,3,ig)
	tmp4 = sq3(1,4,ig)
	tmp5 = sq3(1,5,ig)
	tmp6 = sq3(1,6,ig)
	tmp7 = sq3(1,7,ig)
	tmp8 = sq3(1,8,ig)
	tmp9 = sq3(1,9,ig)
	tmp10 = sq3(1,10,ig)
	tmp11 = sq3(1,11,ig)
	tmp12 = sq3(1,12,ig)
	tmp13 = sq3(1,13,ig)
	tmp14 = sq3(1,14,ig)
        tmp15 = sq3(1,15,ig)
        tmp16 = sq3(1,16,ig)
        tmp17 = sq3(1,17,ig)
        tmp18 = sq3(1,18,ig)

	xr11 =
     &	+ tmp1 * x1(1, 1)  + tmp2 * x1(1, 2)
     &	+ tmp3 * x1(1, 3)  + tmp4 * x1(1, 4)
     &	+ tmp5 * x1(1, 5)  + tmp6 * x1(1, 6)
     &	+ tmp7 * x1(1, 7)  + tmp8 * x1(1, 8)
     &	+ tmp9 * x1(1, 9)  + tmp10 * x1(1, 10)
     &	+ tmp11 * x1(1, 11)+ tmp12 * x1(1, 12)
     &	+ tmp13 * x1(1, 13)+ tmp14 * x1(1, 14)
     &	+ tmp15 * x1(1, 15)+ tmp16 * x1(1, 16)
     &  + tmp17 * x1(1, 17)+ tmp18 * x1(1, 18) 

	xr12 =
     &	+ tmp1 * x1(2, 1) + tmp2 * x1(2, 2)
     &	+ tmp3 * x1(2, 3) + tmp4 * x1(2, 4)
     &	+ tmp5 * x1(2, 5) + tmp6 * x1(2, 6)
     &	+ tmp9 * x1(2, 9)  + tmp10 * x1(2, 10)
     &	+ tmp11 * x1(2, 11)+ tmp12 * x1(2, 12)
     &	+ tmp13 * x1(2, 13)+ tmp14 * x1(2, 14)
     &	+ tmp15 * x1(2, 15)+ tmp16 * x1(2, 16)
     &  + tmp17 * x1(2, 17)+ tmp18 * x1(2, 18)

	xr13 =
     &	+ tmp1 * x1(3, 1) + tmp2 * x1(3, 2)
     &	+ tmp3 * x1(3, 3) + tmp4 * x1(3, 4)
     &	+ tmp5 * x1(3, 5) + tmp6 * x1(3, 6)
     &	+ tmp9 * x1(3, 9)  + tmp10 * x1(3, 10)
     &	+ tmp11 * x1(3, 11)+ tmp12 * x1(3, 12)
     &	+ tmp13 * x1(3, 13)+ tmp14 * x1(3, 14)
     &	+ tmp15 * x1(3, 15)+ tmp16 * x1(3, 16)
     &  + tmp17 * x1(3, 17)+ tmp18 * x1(3, 18)


	tmp1 = sq3(2,1,ig)
	tmp2 = sq3(2,2,ig)
	tmp3 = sq3(2,3,ig)
	tmp4 = sq3(2,4,ig)
	tmp5 = sq3(2,5,ig)
	tmp6 = sq3(2,6,ig)
	tmp7 = sq3(2,7,ig)
	tmp8 = sq3(2,8,ig)
	tmp9 = sq3(2,9,ig)
	tmp10 = sq3(2,10,ig)
	tmp11 = sq3(2,11,ig)
	tmp12 = sq3(2,12,ig)
	tmp13 = sq3(2,13,ig)
	tmp14 = sq3(2,14,ig)
	tmp15 = sq3(2,15,ig)
        tmp16 = sq3(2,16,ig)
        tmp17 = sq3(2,17,ig)
        tmp18 = sq3(2,18,ig)

	xr21 =
     &	+ tmp1  * x1(1, 1) + tmp2 * x1(1, 2)
     &	+ tmp3  * x1(1, 3) + tmp4 * x1(1, 4)
     &	+ tmp5  * x1(1, 5) + tmp6 * x1(1, 6)
     &	+ tmp7  * x1(1, 7) + tmp8 * x1(1, 8)
     &	+ tmp9  * x1(1, 9) + tmp10 * x1(1, 10)
     &	+ tmp11 * x1(1, 11) + tmp12 * x1(1, 12)
     &	+ tmp13 * x1(1, 13) + tmp14 * x1(1, 14)
     &	+ tmp15 * x1(1, 15) + tmp16 * x1(1, 16)
     &  + tmp17 * x1(1, 17) + tmp18 * x1(1, 18) 

	xr22 =
     &	+ tmp1 * x1(2, 1) + tmp2 * x1(2, 2)
     &	+ tmp3 * x1(2, 3) + tmp4 * x1(2, 4)
     &	+ tmp5 * x1(2, 5) + tmp6 * x1(2, 6)
     &	+ tmp7  * x1(2, 7) + tmp8 * x1(2, 8)
     &	+ tmp9  * x1(2, 9) + tmp10 * x1(2, 10)
     &	+ tmp11 * x1(2, 11) + tmp12 * x1(2, 12)
     &	+ tmp13 * x1(2, 13) + tmp14 * x1(2, 14)
     &	+ tmp15 * x1(2, 15) + tmp16 * x1(2, 16)
     &  + tmp17 * x1(2, 17) + tmp18 * x1(2, 18) 


	xr23 =
     &	+ tmp1 * x1(3, 1) + tmp2 * x1(3, 2)
     &	+ tmp3 * x1(3, 3) + tmp4 * x1(3, 4)
     &	+ tmp5 * x1(3, 5) + tmp6 * x1(3, 6)
     &	+ tmp7  * x1(3, 7) + tmp8 * x1(3, 8)
     &	+ tmp9  * x1(3, 9) + tmp10 * x1(3, 10)
     &	+ tmp11 * x1(3, 11) + tmp12 * x1(3, 12)
     &	+ tmp13 * x1(3, 13) + tmp14 * x1(3, 14)
     &	+ tmp15 * x1(3, 15) + tmp16 * x1(3, 16)
     &  + tmp17 * x1(3, 17) + tmp18 * x1(3, 18)

	tmp1 = sq3(3,1,ig)
	tmp2 = sq3(3,2,ig)
	tmp3 = sq3(3,3,ig)
	tmp4 = sq3(3,4,ig)
	tmp5 = sq3(3,5,ig)
	tmp6 = sq3(3,6,ig)
	tmp7 = sq3(3,7,ig)
	tmp8 = sq3(3,8,ig)
	tmp9 = sq3(3,9,ig)
	tmp10 = sq3(3,10,ig)
	tmp11 = sq3(3,11,ig)
	tmp12 = sq3(3,12,ig)
	tmp13 = sq3(3,13,ig)
	tmp14 = sq3(3,14,ig)
	tmp15 = sq3(3,15,ig)
        tmp16 = sq3(3,16,ig)
        tmp17 = sq3(3,17,ig)
        tmp18 = sq3(3,18,ig)


	xr31 =
     &	+ tmp1 * x1(1, 1) + tmp2 * x1(1, 2)
     &	+ tmp3 * x1(1, 3) + tmp4 * x1(1, 4)
     &	+ tmp5 * x1(1, 5) + tmp6 * x1(1, 6)
     &	+ tmp7 * x1(1, 7) + tmp8 * x1(1, 8)
     &	+ tmp9 * x1(1, 9) + tmp10 * x1(1, 10)
     &	+ tmp11 * x1(1, 11) + tmp12 * x1(1, 12)
     &	+ tmp13 * x1(1, 13) + tmp14 * x1(1, 14)
     &	+ tmp15 * x1(1, 15) + tmp16 * x1(1, 16)
     &  + tmp17 * x1(1, 17) + tmp18 * x1(1, 18)



	xr32 =
     &	+ tmp1 * x1(2, 1) + tmp2 * x1(2, 2)
     &	+ tmp3 * x1(2, 3) + tmp4 * x1(2, 4)
     &	+ tmp5 * x1(2, 5) + tmp6 * x1(2, 6)
     &	+ tmp7 * x1(2, 7) + tmp8 * x1(2, 8)
     &	+ tmp9 * x1(2, 9) + tmp10 * x1(2, 10)
     &	+ tmp11 * x1(2, 11) + tmp12 * x1(2, 12)
     &	+ tmp13 * x1(2, 13) + tmp14 * x1(2, 14)
     &	+ tmp15 * x1(2, 15) + tmp16 * x1(2, 16)
     &  + tmp17 * x1(2, 17) + tmp18 * x1(2, 18)


	xr33 =
     &	+ tmp1 * x1(3, 1) + tmp2 * x1(3, 2)
     &	+ tmp3 * x1(3, 3) + tmp4 * x1(3, 4)
     &	+ tmp5 * x1(3, 5) + tmp6 * x1(3, 6)
     &	+ tmp7 * x1(3, 7) + tmp8 * x1(3, 8)
     &	+ tmp9 * x1(3, 9) + tmp10 * x1(3, 10)
     &	+ tmp11 * x1(3, 11) + tmp12 * x1(3, 12)
     &	+ tmp13 * x1(3, 13) + tmp14 * x1(3, 14)
     &	+ tmp15 * x1(3, 15) + tmp16 * x1(3, 16)
     &  + tmp17 * x1(3, 17) + tmp18 * x1(3, 18)



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
!            sh(0,i) = sh0(0,i)
!            sh(1,i) = (cof11*sh0(1,i)+cof12*sh0(2,i)+cof13*sh0(3,i))/det
!            sh(2,i) = (cof21*sh0(1,i)+cof22*sh0(2,i)+cof23*sh0(3,i))/det
!            sh(3,i) = (cof31*sh0(1,i)+cof32*sh0(2,i)+cof33*sh0(3,i))/det

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

        sh(0,7)=sq3(0,7,ig)
	sh(1,7)=sq3(1,7,ig)*cf11+sq3(2,7,ig)*cf12+sq3(3,7,ig)*cf13
	sh(2,7)=sq3(1,7,ig)*cf21+sq3(2,7,ig)*cf22+sq3(3,7,ig)*cf23
	sh(3,7)=sq3(1,7,ig)*cf31+sq3(2,7,ig)*cf32+sq3(3,7,ig)*cf33

        sh(0,8)=sq3(0,8,ig)
	sh(1,8)=sq3(1,8,ig)*cf11+sq3(2,8,ig)*cf12+sq3(3,8,ig)*cf13
	sh(2,8)=sq3(1,8,ig)*cf21+sq3(2,8,ig)*cf22+sq3(3,8,ig)*cf23
	sh(3,8)=sq3(1,8,ig)*cf31+sq3(2,8,ig)*cf32+sq3(3,8,ig)*cf33

        sh(0,9)=sq3(0,9,ig)
	sh(1,9)=sq3(1,9,ig)*cf11+sq3(2,9,ig)*cf12+sq3(3,9,ig)*cf13
	sh(2,9)=sq3(1,9,ig)*cf21+sq3(2,9,ig)*cf22+sq3(3,9,ig)*cf23
	sh(3,9)=sq3(1,9,ig)*cf31+sq3(2,9,ig)*cf32+sq3(3,9,ig)*cf33

        sh(0,10)=sq3(0,10,ig)
	sh(1,10)=sq3(1,10,ig)*cf11+sq3(2,10,ig)*cf12+sq3(3,10,ig)*cf13
	sh(2,10)=sq3(1,10,ig)*cf21+sq3(2,10,ig)*cf22+sq3(3,10,ig)*cf23
	sh(3,10)=sq3(1,10,ig)*cf31+sq3(2,10,ig)*cf32+sq3(3,10,ig)*cf33

        sh(0,11)=sq3(0,11,ig)
	sh(1,11)=sq3(1,11,ig)*cf11+sq3(2,11,ig)*cf12+sq3(3,11,ig)*cf13
	sh(2,11)=sq3(1,11,ig)*cf21+sq3(2,11,ig)*cf22+sq3(3,11,ig)*cf23
	sh(3,11)=sq3(1,11,ig)*cf31+sq3(2,11,ig)*cf32+sq3(3,11,ig)*cf33

        sh(0,12)=sq3(0,12,ig)
	sh(1,12)=sq3(1,12,ig)*cf11+sq3(2,12,ig)*cf12+sq3(3,12,ig)*cf13
	sh(2,12)=sq3(1,12,ig)*cf21+sq3(2,12,ig)*cf22+sq3(3,12,ig)*cf23
	sh(3,12)=sq3(1,12,ig)*cf31+sq3(2,12,ig)*cf32+sq3(3,12,ig)*cf33

        sh(0,13)=sq3(0,13,ig)
	sh(1,13)=sq3(1,13,ig)*cf11+sq3(2,13,ig)*cf12+sq3(3,13,ig)*cf13
	sh(2,13)=sq3(1,13,ig)*cf21+sq3(2,13,ig)*cf22+sq3(3,13,ig)*cf23
	sh(3,13)=sq3(1,13,ig)*cf31+sq3(2,13,ig)*cf32+sq3(3,13,ig)*cf33

        sh(0,14)=sq3(0,14,ig)
	sh(1,14)=sq3(1,14,ig)*cf11+sq3(2,14,ig)*cf12+sq3(3,14,ig)*cf13
	sh(2,14)=sq3(1,14,ig)*cf21+sq3(2,14,ig)*cf22+sq3(3,14,ig)*cf23
	sh(3,14)=sq3(1,14,ig)*cf31+sq3(2,14,ig)*cf32+sq3(3,14,ig)*cf33

        sh(0,15)=sq3(0,15,ig)
	sh(1,15)=sq3(1,15,ig)*cf11+sq3(2,15,ig)*cf12+sq3(3,15,ig)*cf13
	sh(2,15)=sq3(1,15,ig)*cf21+sq3(2,15,ig)*cf22+sq3(3,15,ig)*cf23
	sh(3,15)=sq3(1,15,ig)*cf31+sq3(2,15,ig)*cf32+sq3(3,15,ig)*cf33

        sh(0,16)=sq3(0,16,ig)
        sh(1,16)=sq3(1,16,ig)*cf11+sq3(2,16,ig)*cf12+sq3(3,16,ig)*cf13
        sh(2,16)=sq3(1,16,ig)*cf21+sq3(2,16,ig)*cf22+sq3(3,16,ig)*cf23
        sh(3,16)=sq3(1,16,ig)*cf31+sq3(2,16,ig)*cf32+sq3(3,16,ig)*cf33

        sh(0,17)=sq3(0,17,ig)
        sh(1,17)=sq3(1,17,ig)*cf11+sq3(2,17,ig)*cf12+sq3(3,17,ig)*cf13
        sh(2,17)=sq3(1,17,ig)*cf21+sq3(2,17,ig)*cf22+sq3(3,17,ig)*cf23
        sh(3,17)=sq3(1,17,ig)*cf31+sq3(2,17,ig)*cf32+sq3(3,17,ig)*cf33

        sh(0,18)=sq3(0,18,ig)
        sh(1,18)=sq3(1,18,ig)*cf11+sq3(2,18,ig)*cf12+sq3(3,18,ig)*cf13
        sh(2,18)=sq3(1,18,ig)*cf21+sq3(2,18,ig)*cf22+sq3(3,18,ig)*cf23
        sh(3,18)=sq3(1,18,ig)*cf31+sq3(2,18,ig)*cf32+sq3(3,18,ig)*cf33

