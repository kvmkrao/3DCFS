        subroutine sf315o(ss, x, det, sh, level, xr)
        implicit none
!       include 'global.h'
	integer nsd,nen
	parameter (nen=18)
        parameter(nsd= 3)
        real*8 sh0(9,nen)
        real*8 sh(0:9,nen)
        real*8 xr(nsd,nsd), x(nsd,nen), ss(nsd)
!        real*8 cof11,cof12,cof13
!        real*8 cof21,cof22,cof23
!        real*8 cof31,cof32,cof33
        integer level, i
        real*8 t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12
	real*8 t13,t14,t15

	real*8 X1x,X1e, X1z,Y1x, Y1e, Y1z, Z1x, Z1e, Z1z
	real*8 X1xx, X1ee, X1zz, X1xe, X1xz, X1ez
	real*8 Y1xx, Y1ee, Y1zz, Y1xe, Y1xz, Y1ez
	real*8 Z1xx, Z1ee, Z1zz, Z1xe, Z1xz, Z1ez 

	real*8 det,detx, dete, detz
	real*8 x1,x2,x3,e1,e2,e3,z1,z2,z3

	real*8 x2Xx,x2Xe,x2Xz,x2Yx,x2Ye,x2Yz,x2Zx,x2Ze,x2Zz
	real*8 x3XX,x3YY,x3ZZ,x3XY,x3YZ,x3ZX 
	real*8 e2Xx,e2Xe,e2Xz,e2Yx,e2Ye,e2Yz,e2Zx,e2Ze,e2Zz
	real*8 e3XX,e3YY,e3ZZ,e3XY,e3YZ,e3ZX 
	real*8 z2Xx,z2Xe,z2Xz,z2Yx,z2Ye,z2Yz,z2Zx,z2Ze,z2Zz
	real*8 z3XX,z3YY,z3ZZ,z3XY,z3YZ,z3ZX 
	real*8 sum0, sum1,sum2,sum3
	real*8 sum11,sum22,sum33,sum12,sum13,sum23

	sum0 = 0.0d0
	sum1 = 0.0d0
	sum2 = 0.0d0
	sum3 = 0.0d0
	sum11 = 0.0d0
	sum22 = 0.0d0
	sum33 = 0.0d0
	sum12 = 0.0d0
	sum13 = 0.0d0
	sum23 = 0.0d0

!.......Temporary variables.
!        real*8  tmp1,tmp2,tmp3

!!        write(105,*) ss(1),ss(2),ss(3)
!!        write(105,*) xx(1,1),xx(2,1),xx(3,1)
!!        write(105,*) xx(1,2),xx(2,2),xx(3,2)
!!        write(105,*) xx(1,3),xx(2,3),xx(3,3)
!!        write(105,*) xx(1,4),xx(2,4),xx(3,4)
!!        write(105,*) xx(1,5),xx(2,5),xx(3,5)
!!        write(105,*) xx(1,6),xx(2,6),xx(3,6)
!!        write(105,*)
!!        return

!.......Using the temp variables to reduce redundant calculations.
        t1 = 2.0*ss(1)-1.0
        t2 = 2.0*ss(2)-1.0
        t3 = 1.0-ss(1)-ss(2)
        t4 = ss(3)-1.0
        t5 = 1.0-ss(3)*ss(3)
        t6 = ss(3)+1.0
	t7 = 4.0*ss(1)-1.0
	t8 = -3.0+4.0*ss(1)+4.0*ss(2)
	t9 = 4.0*ss(2)-1.0
	t10 = 1.0-ss(3)*ss(3)
	t11 = ss(3)*(ss(3)-1.0)/2.0
	t12 = ss(3)*(ss(3)+1.0)/2.0
	t13 = (2.0*ss(3)-1.0)/2.0
	t14 = (2.0*ss(3)+1.0)/2.0
	t15 = -2.0*ss(3)

         sh(0,1) = ss(1)*t1*t11
	 sh(0,2) = ss(2)*t2*t11
	 sh(0,3) = (1.0-2.0*ss(1)-2.0*ss(2))*t3*t11
         sh(0,4) = 4.0*ss(1)*ss(2)*t11
	 sh(0,5) = 4.0*ss(1)*t3*t11
	 sh(0,6) = 4.0*ss(2)*t3*t11

         sh(0,7)  = ss(1)*t1* t10
         sh(0,8)  = ss(2)*t2* t10
         sh(0,9)  = (1.0-2.0*ss(1)-2.0*ss(2))*t3* t10
         sh(0,10) = 4.0*ss(1)*ss(2)* t10
         sh(0,11) = 4.0*ss(1)*t3*t10
         sh(0,12) = 4.0*ss(2)*t3*t10

         sh(0,13) = ss(1)*t1*t12
         sh(0,14) = ss(2)*t2*t12
         sh(0,15) = (1.0-2.0*ss(1)-2.0*ss(2))*t3*t12
         sh(0,16) = 4.0*ss(1)*ss(2)*t12
         sh(0,17) = 4.0*ss(1)*t3*t12
         sh(0,18) = 4.0*ss(2)*t3*t12

!	,xi 
         sh0(1,1) = t7*t11
         sh0(1,2) = 0.0d0
         sh0(1,3) = t8*t11
         sh0(1,4) =  4.0*ss(2)*t11
         sh0(1,5) = (4.0-8.0*ss(1)-4.0*ss(2))*t11
         sh0(1,6) = -4.0*ss(2)*t11
         sh0(1,7) = t7*t10
         sh0(1,8) = 0.0d0
         sh0(1,9) = t8*t10
         sh0(1,10) =  4.0*ss(2)*t10
         sh0(1,11) = (4.0-8.0*ss(1)-4.0*ss(2))*t10
         sh0(1,12) = -4.0*ss(2)*t10
         sh0(1,13) = t7*t12
         sh0(1,14) = 0.0d0
         sh0(1,15) = t8*t12
         sh0(1,16) =  4.0*ss(2)*t12
         sh0(1,17) = (4.0-8.0*ss(1)-4.0*ss(2))*t12
         sh0(1,18) = -4.0*ss(2)*t12

!	,eta	
         sh0(2,1) = 0.0d0
         sh0(2,2) = t9*t11
         sh0(2,3) = t8*t11
         sh0(2,4) =  4.0*ss(1)*t11
         sh0(2,5) = -4.0*ss(1)*t11
         sh0(2,6) = (4.0-8.0*ss(2)-4.0*ss(1))*t11
         sh0(2,7) = 0.0d0
         sh0(2,8) = t9*t10
         sh0(2,9) = t8*t10
         sh0(2,10) =  4.0*ss(1)*t10
         sh0(2,11) = -4.0*ss(1)*t10
         sh0(2,12) = (4.0-8.0*ss(2)-4.0*ss(1))*t10
         sh0(2,13) = 0.0d0
         sh0(2,14) = t9*t12
         sh0(2,15) = t8*t12
         sh0(2,16) =  4.0*ss(1)*t12
         sh0(2,17) = -4.0*ss(1)*t12
         sh0(2,18) = (4.0-8.0*ss(2)-4.0*ss(1))*t12

! 	,zeta
         sh0(3,1) = ss(1)*t1*t13
         sh0(3,2) = ss(2)*t2*t13
         sh0(3,3) = (1.0-2.0*ss(1)-2.0*ss(2))*t3*t13
         sh0(3,4) = 4.0*ss(1)*ss(2)*t13
         sh0(3,5) = 4.0*ss(1)*t3*t13
         sh0(3,6) = 4.0*ss(2)*t3*t13
         sh0(3,7) = ss(1)*t1*t15
         sh0(3,8) = ss(2)*t2*t15
         sh0(3,9) = (1.0-2.0*ss(1)-2.0*ss(2))*t3*t15
         sh0(3,10) = 4.0*ss(1)*ss(2)*t15
         sh0(3,11) = 4.0*ss(1)*t3*t15
         sh0(3,12) = 4.0*ss(2)*t3*t15
         sh0(3,13) = ss(1)*t1*t14
         sh0(3,14) = ss(2)*t2*t14
         sh0(3,15) = (1.0-2.0*ss(1)-2.0*ss(2))*t3*t14
         sh0(3,16) = 4.0*ss(1)*ss(2)*t14
         sh0(3,17) = 4.0*ss(1)*t3*t14
         sh0(3,18) = 4.0*ss(2)*t3*t14

        do i=1,nen
        sum0 = sum0 + sh(0,i)
	sum1 = sum1 + sh0(1,i)
	sum2 = sum2 + sh0(2,i)
	sum3 = sum3 + sh0(3,i)
        end do
!	write(6,*)"SUM0 SUM1 SUM2 SUM3"
!	write(6,*)sum0, sum1, sum2, sum3


! 	,xi xi
        sh0(4,1) = 4.0*t11
        sh0(4,2) = 0.0d0
        sh0(4,3) = 4.0*t11
        sh0(4,4) = 0.0d0
        sh0(4,5) = 0.0d0
        sh0(4,6) = -8.0*t11
        sh0(4,7) = 4.0*t10
        sh0(4,8) = 0.0d0
        sh0(4,9) = 4.0*t10
        sh0(4,10) = 0.0d0
        sh0(4,11) = 0.0d0
        sh0(4,12) = -8.0*t10
        sh0(4,13) = 4.0*t12
        sh0(4,14) = 0.0d0
        sh0(4,15) = 4.0*t12
        sh0(4,16) = 0.0d0
        sh0(4,17) = 0.0d0
        sh0(4,18) = -8.0*t12

!	,eta eta
 
        sh0(5,1) = 0.0d0
        sh0(5,2) = 4.0*t11
        sh0(5,3) = 4.0*t11
        sh0(5,4) = 0.0d0
        sh0(5,5) = 0.0d0
        sh0(5,6) = -8.0*t11
        sh0(5,7) = 0.0d0
        sh0(5,8) = 4.0*t10
        sh0(5,9) = 4.0*t10
        sh0(5,10) = 0.0d0
        sh0(5,11) = 0.0d0
        sh0(5,12) = -8.0*t10
        sh0(5,13) = 0.0d0
        sh0(5,14) = 4.0*t12
        sh0(5,15) = 4.0*t12
        sh0(5,16) = 0.0d0
        sh0(5,17) = 0.0d0
        sh0(5,18) = -8.0*t12

!	,zeta zeta
         sh0(6,1) = ss(1)*t1
         sh0(6,2) = ss(2)*t2
         sh0(6,3) = (1.0-2.0*ss(1)-2.0*ss(2))*t3
         sh0(6,4) = 4.0*ss(1)*ss(2)
         sh0(6,5) = 4.0*ss(1)*t3
         sh0(6,6) = 4.0*ss(2)*t3

         sh0(6,7)  = -ss(1)*t1*2.0d0
         sh0(6,8)  = -ss(2)*t2*2.0d0
         sh0(6,9)  = -(1.0-2.0*ss(1)-2.0*ss(2))*t3*2.0d0
         sh0(6,10) = -4.0*ss(1)*ss(2)*2.0d0
         sh0(6,11) = -4.0*ss(1)*t3*2.0d0
         sh0(6,12) = -4.0*ss(2)*t3*2.0d0

         sh0(6,13) = ss(1)*t1
         sh0(6,14) = ss(2)*t2
         sh0(6,15) = (1.0-2.0*ss(1)-2.0*ss(2))*t3
         sh0(6,16) = 4.0*ss(1)*ss(2)
         sh0(6,17) = 4.0*ss(1)*t3
         sh0(6,18) = 4.0*ss(2)*t3

!	,xi eta
         sh0(7,1) = 0.0d0
         sh0(7,2) = 0.0d0
         sh0(7,3) = 4.0*t11
         sh0(7,4) = 4.0*t11
         sh0(7,5) = -4.0*t11
         sh0(7,6) = -4.0*t11
         sh0(7,7) = 0.0d0
         sh0(7,8) = 0.0d0
         sh0(7,9) = 4.0*t10
         sh0(7,10) = 4.0*t10
         sh0(7,11) = -4.0*t10
         sh0(7,12) = -4.0*t10
         sh0(7,13) = 0.0d0
         sh0(7,14) = 0.0d0
         sh0(7,15) = 4.0*t12
         sh0(7,16) = 4.0*t12
         sh0(7,17) = -4.0*t12
         sh0(7,18) = -4.0*t12

!	,xi zeta
         sh0(8,1) = t7*t13
         sh0(8,2) = 0.0d0
         sh0(8,3) = t8*t13
         sh0(8,4) =  4.0*ss(2)*t13
         sh0(8,5) = (4.0-8.0*ss(1)-4.0*ss(2))*t13
         sh0(8,6) = -4.0*ss(2)*t13
         sh0(8,7) = t7*t15
         sh0(8,8) = 0.0d0
         sh0(8,9) = t8*t15
         sh0(8,10) =  4.0*ss(2)*t15
         sh0(8,11) = (4.0-8.0*ss(1)-4.0*ss(2))*t15
         sh0(8,12) = -4.0*ss(2)*t15
         sh0(8,13) = t7*t14
         sh0(8,14) = 0.0d0
         sh0(8,15) = t8*t14
         sh0(8,16) =  4.0*ss(2)*t14
         sh0(8,17) = (4.0-8.0*ss(1)-4.0*ss(2))*t14
         sh0(8,18) = -4.0*ss(2)*t14

!	,eta zeta
         sh0(9,1) = 0.0d0
         sh0(9,2) = t9*t13
         sh0(9,3) = (-3.0+4.0*ss(2)+4.0*ss(1))*t13
         sh0(9,4) =  4.0*ss(1)*t13
         sh0(9,5) = -4.0*ss(1)*t13
         sh0(9,6) = (4.0-8.0*ss(2)-4.0*ss(1))*t13
         sh0(9,7) = 0.0d0
         sh0(9,8) = t9*t15
         sh0(9,9) = (-3.0+4.0*ss(2)+4.0*ss(1))*t15
         sh0(9,10) =  4.0*ss(1)*t15
         sh0(9,11) = -4.0*ss(1)*t15
         sh0(9,12) = (4.0-8.0*ss(2)-4.0*ss(1))*t15
         sh0(9,13) = 0.0d0
         sh0(9,14) = t9*t14
         sh0(9,15) = (-3.0+4.0*ss(2)+4.0*ss(1))*t14
         sh0(9,16) =  4.0*ss(1)*t14
         sh0(9,17) = -4.0*ss(1)*t14
         sh0(9,18) = (4.0-8.0*ss(2)-4.0*ss(1))*t14

        do i=1,nen
	        sum11 = sum11 + sh0(4,i)
	        sum22 = sum22 + sh0(5,i)
	        sum33 = sum33 + sh0(6,i)
	        sum12 = sum12 + sh0(7,i)
	        sum13 = sum13 + sh0(8,i)
	        sum23 = sum23 + sh0(9,i)
        end do

!        write(6,*)"SUM11 SUM22 SUM33 "
!        write(6,*)sum11, sum22, sum33 
!	write(6,*)"SUM12 SUM13 SUM23"
!	write(6,*)sum12, sum13, sum23

        X1x = 0.0d0
        X1e = 0.0d0
	X1z = 0.0d0
        Y1x = 0.0d0
        Y1e = 0.0d0
	Y1z = 0.0d0
	Z1x = 0.0d0
	Z1e = 0.0d0
	Z1z = 0.0d0

        X1xx = 0.0d0
        X1ee = 0.0d0
	X1zz = 0.0d0
        X1xe = 0.0d0
	X1xz = 0.0d0
	X1ez = 0.0d0

        Y1xx = 0.0d0
        Y1ee = 0.0d0
	Y1zz = 0.0d0
        Y1xe = 0.0d0
	Y1xz = 0.0d0
	Y1ez = 0.0d0

	Z1xx = 0.0d0
        Z1ee = 0.0d0
	Z1zz = 0.0d0
        Z1xe = 0.0d0
	Z1xz = 0.0d0
	Z1ez = 0.0d0

      do 210 i = 1,nen

        X1x = X1x + sh0(1,i) * x(1,i)
        X1e = X1e + sh0(2,i) * x(1,i)
	X1z = X1z + sh0(3,i) * x(1,i)

        Y1x = Y1x + sh0(1,i) * x(2,i)
        Y1e = Y1e + sh0(2,i) * x(2,i)
	Y1z = Y1z + sh0(3,i) * x(2,i)

        Z1x = Z1x + sh0(1,i) * x(3,i)
        Z1e = Z1e + sh0(2,i) * x(3,i)
	Z1z = Z1z + sh0(3,i) * x(3,i)

	X1xx = X1xx + sh0(4,i) *x(1,i)
	X1ee = X1ee + sh0(5,i) *x(1,i)
	X1zz = X1zz + sh0(6,i) *x(1,i)
	X1xe = X1xe + sh0(7,i) *x(1,i)
	X1xz = X1xz + sh0(8,i) *x(1,i)
	X1ez = X1ez + sh0(9,i) *x(1,i)

	Y1xx = Y1xx + sh0(4,i) *x(2,i)
	Y1ee = Y1ee + sh0(5,i) *x(2,i)
	Y1zz = Y1zz + sh0(6,i) *x(2,i)
	Y1xe = Y1xe + sh0(7,i) *x(2,i)
	Y1xz = Y1xz + sh0(8,i) *x(2,i)
	Y1ez = Y1ez + sh0(9,i) *x(2,i)

	Z1xx = Z1xx + sh0(4,i) *x(3,i)
	Z1ee = Z1ee + sh0(5,i) *x(3,i)
	Z1zz = Z1zz + sh0(6,i) *x(3,i)
	Z1xe = Z1xe + sh0(7,i) *x(3,i)
	Z1xz = Z1xz + sh0(8,i) *x(3,i)
	Z1ez = Z1ez + sh0(9,i) *x(3,i)

 210   continue

	det = X1x*(Y1e*Z1z-Z1e*Y1z) - X1e*(Y1x*Z1z-Z1x*Y1z) 
     $				    + X1z*(Y1x*Z1e-Z1x*Y1e)

!	write(6,*) "DET of J"
!	write(6,*) det

        detx = X1xx * (Y1e*Z1z-Z1e*Y1z) + X1x*(Y1xe*Z1z+Y1e*Z1xz)
     $                                  - X1x*(Z1xe*Y1z+Z1e*Y1xz)
     $        -X1xe * (Y1x*Z1z-Z1x*Y1z) - X1e*(Y1xx*Z1z+Y1x*Z1xz)
     $                                  + X1e*(Z1xx*Y1z+Z1x*Y1xz)
     $        +X1xz * (Y1x*Z1e-Z1x*Y1e) + X1z*(Y1xx*Z1e+Y1x*Z1xe)
     $                                  - X1z*(Z1xx*Y1e+Z1x*Y1xe)

	dete = X1xe * (Y1e*Z1z-Z1e*Y1z) + X1x*(Y1ee*Z1z+Y1e*Z1ez)
     $                                  - X1x*(Z1ee*Y1z+Z1e*Y1ez)
     $        -X1ee * (Y1x*Z1z-Z1x*Y1z) - X1e*(Y1xe*Z1z+Y1x*Z1ez)
     $                                  + X1e*(Z1xe*Y1z+Z1x*Y1ez)
     $        +X1ez * (Y1x*Z1e-Z1x*Y1e) + X1z*(Y1xe*Z1e+Y1x*Z1ee)
     $                                  - X1z*(Z1xe*Y1e+Z1x*Y1ee)

	detz = X1xz * (Y1e*Z1z-Z1e*Y1z) + X1x*(Y1ez*Z1z+Y1e*Z1zz)
     $                                  - X1x*(Z1ez*Y1z+Z1e*Y1zz)
     $        -X1ez * (Y1x*Z1z-Z1x*Y1z) - X1e*(Y1xz*Z1z+Y1x*Z1zz)
     $                                  + X1e*(Z1xz*Y1z+Z1x*Y1zz)
     $        +X1zz * (Y1x*Z1e-Z1x*Y1e) + X1z*(Y1xz*Z1e+Y1x*Z1ez)
     $                                  - X1z*(Z1xz*Y1e+Z1x*Y1ez)

	x1 =  (Y1e*Z1z-Z1e*Y1z)/det 
	x2 = -(X1e*Z1z-Z1e*X1z)/det  
	x3 =  (X1e*Y1z-Y1e*X1z)/det
	e1 = -(Y1x*Z1z-Z1x*Y1z)/det
	e2 =  (X1x*Z1z-Z1x*X1z)/det
	e3 = -(X1x*Y1z-Y1x*X1z)/det
	z1 =  (Y1x*Z1e-Z1x*Y1e)/det
	z2 = -(X1x*Z1e-Z1x*X1e)/det
	z3 =  (X1x*Y1e-Y1x*X1e)/det

!	write(6,*)x1,x2,x3,e1,e2,e3,z1,z2,z3

	xr(1,1) =  x1 !x1x
	xr(1,2) =  x2 !x1y
	xr(1,3) =  x3 !x1z
	xr(2,1) =  e1 !e1x
	xr(2,2) =  e2 !e1y
	xr(2,3) =  e3 !e1z
	xr(3,1) =  z1 !z1x
	xr(3,2) =  z2 !z1y
	xr(3,3) =  z3 !z1z

        do i = 1,nen
         sh(1,i) =  sh0(1,i) * x1 + sh0(2,i) * e1 + sh0(3,i)*z1
	 sh(2,i) =  sh0(1,i) * x2 + sh0(2,i) * e2 + sh0(3,i)*z2
	 sh(3,i) =  sh0(1,i) * x3 + sh0(2,i) * e3 + sh0(3,i)*z3
        end do

!	if(level.eq.1) return 

!	check
	x2Xx = (Y1xe* Z1z + Y1e * Z1xz
     $	      - Y1xz* Z1e - Y1z * Z1xe) - x1*detx
	x2Xe = (Y1ee* Z1z + Y1e * Z1ez
     $	      - Y1ez* Z1e - Y1z * Z1ee) - x1*dete 
	x2Xz = (Y1ez* Z1z + Y1e * Z1zz
     $	      - Y1zz* Z1e - Y1z * Z1ez) - x1*detz
!	check

	x2Yx = -(X1xe* Z1z + X1e * Z1xz
     $	        -X1xz* Z1e - X1z * Z1xe) + x2*detx
	x2Ye = -(X1ee* Z1z + X1e * Z1ez
     $	        -X1ez* Z1e - X1z * Z1ee) + x2*dete 
	x2Yz = -(X1ez* Z1z + X1e * Z1zz
     $	        -X1zz* Z1e - X1z * Z1ez) + x2*detz 

	x2Zx = (X1xe* Y1z + X1e * Y1xz
     $	      - X1xz* Y1e - X1z * Y1xe) - x3*detx
	x2Ze = (X1ee* Y1z + X1e * Y1ez
     $	       -X1ez* Y1e - X1z * Y1ee) - x3*dete
	x2Zz = (X1ez* Y1z + X1e * Y1zz 
     $	       -X1zz* Y1e - X1z * Y1ez) - x3*detz
!	check
        x3XX = (x2Xx * x1 + x2Xe * e1 + x2Xz * z1)/det
!	check 
        x3YY = (x2Yx * x2 + x2Ye * e2 + x2Yz * z2)/det
        x3ZZ = (x2Zx * x3 + x2Ze * e3 + x2Zz * z3)/det
        x3XY = (x2Xx * x2 + x2Xe * e2 + x2Xz * z2)/det
        x3YZ = (x2Yx * x3 + x2Ye * e3 + x2Yz * z3)/det
        x3ZX = (x2Zx * x1 + x2Ze * e1 + x2Zz * z1)/det
!	write(6,*)"x3XX, x3YY, x3ZZ, x3XY, x3YZ, x3ZX"
!	write(6,*) x3XX, x3YY, x3ZZ, x3XY, x3YZ, x3ZX

!	check
	e2Xx = -(Y1xx* Z1z + Y1x * Z1xz
     $	       - Y1xz* Z1x - Y1z * Z1xx) + e1*detx
	e2Xe = -(Y1xe* Z1z + Y1x * Z1ez
     $	       - Y1ez* Z1x - Y1z * Z1xe) + e1*dete
	e2Xz = -(Y1xz* Z1z + Y1x * Z1zz
     $	       - Y1zz* Z1x - Y1z * Z1xz) + e1*detz
!	check 

	e2Yx = (X1xx* Z1z + X1x * Z1xz
     $	      - X1xz* Z1x - X1z * Z1xx) - e2*detx
	e2Ye = (X1xe* Z1z + X1x * Z1ez
     $	      - X1ez* Z1x - X1z * Z1xe) - e2*dete
	e2Yz = (X1xz* Z1z + X1x * Z1zz
     $	      - X1zz* Z1x - X1z * Z1xz) - e2*detz

	e2Zx = -(X1xx* Y1z + X1x * Y1xz
     $	       - X1xz* Y1x - X1z * Y1xx) + e3*detx
	e2Ze = -(X1xe* Y1z + X1x * Y1ez
     $	       - X1ez* Y1x - X1z * Y1xe) + e3*dete
	e2Zz = -(X1xz* Y1z + X1x * Y1zz
     $	       - X1zz* Y1x - X1z * Y1xz) + e3*detz

        e3XX = (e2Xx * x1 + e2Xe * e1 + e2Xz * z1)/det
        e3YY = (e2Yx * x2 + e2Ye * e2 + e2Yz * z2)/det
        e3ZZ = (e2Zx * x3 + e2Ze * e3 + e2Zz * z3)/det
        e3XY = (e2Xx * x2 + e2Xe * e2 + e2Xz * z2)/det
        e3YZ = (e2Yx * x3 + e2Ye * e3 + e2Yz * z3)/det
        e3ZX = (e2Zx * x1 + e2Ze * e1 + e2Zz * z1)/det
!	write(6,*)"e3XX, e3YY, e3ZZ, e3XY, e3YZ, e3ZX"
!	write(6,*)e3XX, e3YY, e3ZZ, e3XY, e3YZ, e3ZX


!	check
	z2Xx = (Y1xx* Z1e + Y1x * Z1xe
     $	      - Y1xe* Z1x - Y1e * Z1xx) - z1*detx
	z2Xe = (Y1xe* Z1e + Y1x * Z1ee
     $	      - Y1ee* Z1x - Y1e * Z1xe) - z1*dete
	z2Xz = (Y1xz* Z1e + Y1x * Z1ez
     $	      - Y1ez* Z1x - Y1e * Z1xz) - z1*detz
!	check 
	z2Yx = -(X1xx* Z1e + X1x * Z1xe
     $	       - X1xe* Z1x - X1e * Z1xx) + z2*detx
	z2Ye = -(X1xe* Z1e + X1x * Z1ee
     $	       - X1ee* Z1x - X1e * Z1xe) + z2*dete
	z2Yz = -(X1xz* Z1e + X1x * Z1ez
     $	       - X1ez* Z1x - X1e * Z1xz) + z2*detz

	z2Zx = (X1xx* Y1e + X1x * Y1xe
     $	      - X1xe* Y1x - X1e * Y1xx) - z3*detx
	z2Ze = (X1xe* Y1e + X1x * Y1ee
     $	      - X1ee* Y1x - X1e * Y1xe) - z3*dete
	z2Zz = (X1xz* Y1e + X1x * Y1ez
     $	      - X1ez* Y1x - X1e * Y1xz) - z3*detz

        z3XX = (z2Xx * x1 + z2Xe * e1 + z2Xz * z1)/det
        z3YY = (z2Yx * x2 + z2Ye * e2 + z2Yz * z2)/det
        z3ZZ = (z2Zx * x3 + z2Ze * e3 + z2Zz * z3)/det
        z3XY = (z2Xx * x2 + z2Xe * e2 + z2Xz * z2)/det
        z3YZ = (z2Yx * x3 + z2Ye * e3 + z2Yz * z3)/det
        z3ZX = (z2Zx * x1 + z2Ze * e1 + z2Zz * z1)/det
!	write(6,*)"z3XX, z3YY, z3ZZ, z3XY, z3YZ, z3ZX"
!	write(6,*)z3XX, z3YY, z3ZZ, z3XY, z3YZ, z3ZX

      do 110 i = 1,nen

        sh(4,i) = sh0(4,i) *x1*x1 + sh0(5,i)*e1*e1 
     $		 +sh0(6,i) *z1*z1 
     $           + 2.0* sh0(7,i)*x1*e1
     $           + 2.0* sh0(8,i)*x1*z1
     $           + 2.0* sh0(9,i)*e1*z1
     $           + sh0(1,i)*x3XX
     $           + sh0(2,i)*e3XX
     $           + sh0(3,i)*z3XX
!	write(6,*)"sh0(4,i),sh0(5,i),sh0(6,i),sh0(7,i), sh0(8,i),sh0(9,i)"
!	write(6,*) sh0(4,i),sh0(5,i),sh0(6,i),sh0(7,i), sh0(8,i),sh0(9,i)
!	write(6,*)"x1,e1,z1,x3XX,e3XX,z3XX"
!	write(6,*)x1,e1,z1,x3XX,e3XX,z3XX


       sh(5,i) = sh0(4,i) *x2*x2 + sh0(5,i)*e2*e2
     $		 +sh0(6,i) *z2*z2 
     $           + 2.0* sh0(7,i)*x2*e2
     $           + 2.0* sh0(8,i)*x2*z2
     $           + 2.0* sh0(9,i)*e2*z2
     $           + sh0(1,i)*x3YY
     $           + sh0(2,i)*e3YY
     $           + sh0(3,i)*z3YY

       sh(6,i) = sh0(4,i) *x3*x3 + sh0(5,i)*e3*e3
     $		 +sh0(6,i) *z3*z3 
     $           + 2.0* sh0(7,i)*x3*e3
     $           + 2.0* sh0(8,i)*x3*z3
     $           + 2.0* sh0(9,i)*e3*z3
     $           + sh0(1,i)*x3ZZ
     $           + sh0(2,i)*e3ZZ
     $           + sh0(3,i)*z3ZZ

       sh(7,i) = sh0(4,i) *x1*x2 + sh0(5,i)*e1*e2
     $		 +sh0(6,i) *z1*z2 
     $           + sh0(7,i)*(x1*e2+x2*e1)
     $           + sh0(8,i)*(x1*z2+x2*z1)
     $           + sh0(9,i)*(e1*z2+e2*z1)
     $           + sh0(1,i)*x3XY
     $           + sh0(2,i)*e3XY
     $           + sh0(3,i)*z3XY

      sh(8,i) = sh0(4,i) *x1*x3 + sh0(5,i)*e1*e3
     $		 +sh0(6,i) *z1*z3
     $           + sh0(7,i)*(x1*e3+x3*e1)
     $           + sh0(8,i)*(x1*z3+x3*z1)
     $           + sh0(9,i)*(e1*z3+e3*z1)
     $           + sh0(1,i)*x3ZX
     $           + sh0(2,i)*e2ZX
     $           + sh0(3,i)*z3ZX

      sh(9,i) = sh0(4,i) *x2*x3 + sh0(5,i)*e2*e3
     $		 +sh0(6,i) *z2*z3
     $           + sh0(7,i)*(x2*e3+x3*e2)
     $           + sh0(8,i)*(x2*z3+x3*z2)
     $           + sh0(9,i)*(e2*z3+e3*z2)
     $           + sh0(1,i)*x3YZ
     $           + sh0(2,i)*e2YZ
     $           + sh0(3,i)*z3YZ
110    continue

        if (level .eq. 2 ) return
        return
        end

!        if (level.eq.1) return
!
!        endsubroutine

!***************************************************************************
        subroutine lgau3tq(igaus,ngaus,gaus,weis)
!
        implicit none
        real*8  gaus(3,9),weis(9)
        integer igaus,ngaus
!
        goto (1,2,3,4),igaus
!
 1      continue
        ngaus = 1
!
        gaus(1,1) = + 1.0d0/3.0d0
        gaus(2,1) = + 1.0d0/3.0d0
        gaus(3,1) = + 1.0d0/3.0d0
!
        weis(1) = 1.0d0
        return
!
 2      continue
	ngaus = 3

        gaus(1,1) = + 0.5d0
        gaus(2,1) = + 0.5d0
        gaus(3,1) = + 0.0d0

        gaus(1,2) = + 0.0d0
        gaus(2,2) = + 0.5d0
        gaus(3,2) = + 0.5d0

        gaus(1,3) = + 0.5d0
        gaus(2,3) = + 0.0d0
        gaus(3,3) = + 0.5d0

	weis(1)   = + 1.0d0/3.0d0
	weis(2)   = + 1.0d0/3.0d0
	weis(3)   = + 1.0d0/3.0d0

	return

3      continue
        ngaus = 4
!
        gaus(1,1) = + 1.0d0/3.0d0
        gaus(2,1) = + 1.0d0/3.0d0 
        gaus(3,1) = + 1.0d0/3.0d0

        gaus(1,2) = + 0.60d0
        gaus(2,2) = + 0.20d0
        gaus(3,2) = + 0.20d0

        gaus(1,3) = + 0.20d0
        gaus(2,3) = + 0.60d0
        gaus(3,3) = + 0.20d0

        gaus(1,4) = + 0.20d0
        gaus(2,4) = + 0.20d0
        gaus(3,4) = + 0.60d0

	weis(1)   = -27.0d0/48.0d0
	weis(2)   =  25.0d0/48.0d0
	weis(3)   =  25.0d0/48.0d0
	weis(4)   =  25.0d0/48.0d0
        return
!
!.....  lumped at the nodes
!
 4      continue
        ngaus = 7

        gaus(1,1) = 1.0d0/3.0d0
        gaus(2,1) = 1.0d0/3.0d0 
        gaus(3,1) = 1.0d0/3.0d0

        gaus(1,2) = 0.0597158717d0
        gaus(2,2) = 0.4701420641d0
        gaus(3,2) = 0.4701420641d0

        gaus(1,3) = 0.4701420641d0
        gaus(2,3) = 0.0597158717d0
        gaus(3,3) = 0.4701420641d0

        gaus(1,4) = 0.4701420641d0
        gaus(2,4) = 0.4701420641d0
        gaus(3,4) = 0.0597158717d0

        gaus(1,5) = 0.7974269853d0
        gaus(2,5) = 0.1012865073d0
        gaus(3,5) = 0.1012865073d0

        gaus(1,6) = 0.1012865073d0
        gaus(2,6) = 0.7974269853d0
        gaus(3,6) = 0.1012865073d0

        gaus(1,7) = 0.1012865073d0
        gaus(2,7) = 0.1012865073d0
        gaus(3,7) = 0.7974269853d0

        weis(1) = 0.2250000000d0
        weis(2) = 0.1323941527d0
        weis(3) = 0.1323941527d0
        weis(4) = 0.1323941527d0
        weis(5) = 0.1259391805d0
        weis(6) = 0.1259391805d0
	weis(7) = 0.1259391805d0
!
        return
        end

***************************************************************************
        subroutine sf36o(ss, xx, det, sh, level, xr)
        implicit none
!        include 'global.h'
	integer nsd,nen
	parameter (nen=6)
        parameter(nsd= 3)
        real*8 sh0(0:9,nen), sh(0:3,nen)
        real*8 xr(nsd,nsd), xx(nsd,nen), ss(nsd)
        real*8 cof11,cof12,cof13
        real*8 cof21,cof22,cof23
        real*8 cof31,cof32,cof33
        real*8 det
        integer level, i

!.......Temporary variables.
        real*8  tmp1,tmp2,tmp3

!!        write(105,*) ss(1),ss(2),ss(3)
!!        write(105,*) xx(1,1),xx(2,1),xx(3,1)
!!        write(105,*) xx(1,2),xx(2,2),xx(3,2)
!!        write(105,*) xx(1,3),xx(2,3),xx(3,3)
!!        write(105,*) xx(1,4),xx(2,4),xx(3,4)
!!        write(105,*) xx(1,5),xx(2,5),xx(3,5)
!!        write(105,*) xx(1,6),xx(2,6),xx(3,6)
!!        write(105,*)
!!        return

!.......Using the temp variables to reduce redundant calculations.
        tmp1 = 0.5d0*(1-ss(3))
        tmp2 = 0.5d0*(1+ss(3))
        tmp3 = (1-ss(1)-ss(2))
!.......shape functions Na.

        sh0(0,1) = tmp1*(ss(1))
        sh0(0,2) = tmp1*(ss(2))
        sh0(0,3) = tmp1*tmp3
        sh0(0,4) = tmp2*(ss(1))
        sh0(0,5) = tmp2*(ss(2))
        sh0(0,6) = tmp2*tmp3
!.......shape functions derivatives Na,xi
        sh0(1,1) =  tmp1
        sh0(1,2) =  0.0
        sh0(1,3) = -tmp1
        sh0(1,4) =  tmp2
        sh0(1,5) =  0.0
        sh0(1,6) = -tmp2

!.......shape functions derivatives Na,eta
        sh0(2,1) =  0.0
        sh0(2,2) =  tmp1
        sh0(2,3) = -tmp1
        sh0(2,4) =  0.0
        sh0(2,5) =  tmp2
        sh0(2,6) = -tmp2

!.......shape functions derivatives Na,chi
        sh0(3,1) = -0.5*ss(1)
        sh0(3,2) = -0.5*ss(2)
        sh0(3,3) = -0.5*tmp3
        sh0(3,4) = -sh0(3,1)
        sh0(3,5) = -sh0(3,2)
        sh0(3,6) = -sh0(3,3)

        xr(1,1) = xx(1,1) - xx(1,3)     ! x,xi
        xr(1,2) = xx(1,2) - xx(1,3)     ! x,eta
        xr(1,3) = 0.0                   ! x,chi
        xr(2,1) = xx(2,1) - xx(2,3)     ! y,xi
        xr(2,2) = xx(2,2) - xx(2,3)     ! y,eta
        xr(2,3) = 0.0                   ! y,chi
        xr(3,1) = 0.0                   ! z,xi
        xr(3,2) = 0.0                   ! z,eta
        xr(3,3) = 0.5*(xx(3,4)-xx(3,1)) ! z,chi

!.......determinant of the jacobian matrix
        det = xr(3,3)*(xr(1,1)*xr(2,2)-xr(1,2)*xr(2,1))

        cof11 =  xr(2,2)*xr(3,3)
        cof12 = -xr(2,1)*xr(3,3)
        cof13 =  0.0
        cof21 = -xr(1,2)*xr(3,3)
        cof22 =  xr(1,1)*xr(3,3)
        cof23 =  0.0
        cof31 =  0.0
        cof32 =  0.0
        cof33 =  xr(1,1)*xr(2,2)-xr(1,2)*xr(2,1)
!!        write(105,*) 'nen',nen
!!        write(105,*) cof11,cof12
!!        write(105,*) cof21,cof22
!!        write(105,*) cof33, det
!!        write(105,*)
!!        write(105,*) ss(1),ss(2),ss(3)
!!        write(105,*)
!!        write(105,*) sh0(0,1),sh0(0,2),sh0(0,3)
!!        write(105,*) sh0(0,4),sh0(0,5),sh0(0,6)
!!        write(105,*)
!!        write(105,*) sh0(1,1),sh0(1,2),sh0(1,3)
!!        write(105,*) sh0(1,4),sh0(1,5),sh0(1,6)
!!        write(105,*)
!!        write(105,*) sh0(2,1),sh0(2,2),sh0(2,3)
!!        write(105,*) sh0(2,4),sh0(2,5),sh0(2,6)
!!        write(105,*)
!!        write(105,*) sh0(3,1),sh0(3,2),sh0(3,3)
!!        write(105,*) sh0(3,4),sh0(3,5),sh0(3,6)

!.......derivatives Na,x; Na,y; Na,z
        do i=1,nen
            sh(0,i) = sh0(0,i)
            sh(1,i) = (cof11*sh0(1,i)+cof12*sh0(2,i)+cof13*sh0(3,i))/det
            sh(2,i) = (cof21*sh0(1,i)+cof22*sh0(2,i)+cof23*sh0(3,i))/det
            sh(3,i) = (cof31*sh0(1,i)+cof32*sh0(2,i)+cof33*sh0(3,i))/det
        enddo

        xr(1,1) = cof11/det           
        xr(2,1) = cof12/det           
        xr(3,1) = cof13/det           
        xr(1,2) = cof21/det           
        xr(2,2) = cof22/det           
        xr(3,2) = cof23/det           
        xr(1,3) = cof31/det           
        xr(2,3) = cof32/det           
        xr(3,3) = cof33/det           

!FIXME..debugging
!!        write(105,*)
!!        write(105,*) sh(0,1),sh(0,2),sh(0,3)
!!        write(105,*) sh(0,4),sh(0,5),sh(0,6)
!!        write(105,*)
!!        write(105,*) sh(1,1),sh(1,2),sh(1,3)
!!        write(105,*) sh(1,4),sh(1,5),sh(1,6)
!!        write(105,*)
!!        write(105,*) sh(2,1),sh(2,2),sh(2,3)
!!        write(105,*) sh(2,4),sh(2,5),sh(2,6)
!!        write(105,*)
!!        write(105,*) sh(3,1),sh(3,2),sh(3,3)
!!        write(105,*) sh(3,4),sh(3,5),sh(3,6)
!!        write(105,*)
!!        write(105,*) sh(1,1)+sh(1,4),sh(1,2)+sh(1,5),sh(1,3)+sh(1,6)
!!        write(105,*) sh(2,1)+sh(2,4),sh(2,2)+sh(2,5),sh(2,3)+sh(2,6)
!!        write(105,*) sh(3,1)+sh(3,4),sh(3,2)+sh(3,5),sh(3,3)+sh(3,6)
!!        write(105,*) xx(1,1),xx(2,1),xx(3,1)
!!        write(105,*) xx(1,2),xx(2,2),xx(3,2)
!!        write(105,*) xx(1,3),xx(2,3),xx(3,3)
!!        stop
        
        if (level.eq.1) return

        endsubroutine

***************************************************************************

        subroutine lgau3t(igaus,ngaus,gaus,weis)
!
        implicit none

        real*8  gaus(3,9),weis(9)
        integer igaus,ngaus
!
        goto (1,2,3,4),igaus
!
 1      continue
        ngaus = 1
!
        gaus(1,1) = + 0.333333333333333
        gaus(2,1) = + 0.333333333333333
        gaus(3,1) = + 0.000000000000000
!
        weis(1) = 1.000000000000000
        return
!
 2      continue
 3      continue
        ngaus = 9
!
        gaus(1,1) = + 0.166666666666667d0
        gaus(2,1) = + 0.166666666666667d0
        gaus(3,1) = - sqrt(3.0d0/5)

        gaus(1,2) = + 0.166666666666667d0
        gaus(2,2) = + 0.666666666666667d0
        gaus(3,2) =   gaus(3,1)

        gaus(1,3) = + 0.666666666666667d0
        gaus(2,3) = + 0.166666666666667d0
        gaus(3,3) =   gaus(3,1)

        gaus(1,4) = + 0.166666666666667d0
        gaus(2,4) = + 0.166666666666667d0
        gaus(3,4) = - gaus(3,1)

        gaus(1,5) = + 0.166666666666667d0
        gaus(2,5) = + 0.666666666666667d0
        gaus(3,5) = - gaus(3,1)

        gaus(1,6) = + 0.666666666666667d0
        gaus(2,6) = + 0.166666666666667d0
        gaus(3,6) = - gaus(3,1)

        gaus(1,7) = + 0.166666666666667d0
        gaus(2,7) = + 0.166666666666667d0
        gaus(3,7) =   0.0d0

        gaus(1,8) = + 0.166666666666667d0
        gaus(2,8) = + 0.666666666666667d0
        gaus(3,8) =   0.0d0

        gaus(1,9) = + 0.666666666666667d0
        gaus(2,9) = + 0.166666666666667d0
        gaus(3,9) =   0.0d0

        weis(1) = 5.0d0/54.0
        weis(2) = weis(1)
        weis(3) = weis(1)
        weis(4) = weis(1)
        weis(5) = weis(1)
        weis(6) = weis(1)
        weis(7) = 4.0d0/27.0
        weis(8) = weis(7)
        weis(9) = weis(7)

!
        return
!
!.....  lumped at the nodes
!
 4      continue
        ngaus = 6
        gaus(1,1) = + 1.0d0
        gaus(1,2) = + 0.0d0
        gaus(1,3) = + 0.0d0
        gaus(1,4) = + 1.0d0
        gaus(1,5) = + 0.0d0
        gaus(1,6) = + 0.0d0
        gaus(2,1) = + 0.0d0
        gaus(2,2) = + 1.0d0
        gaus(2,3) = + 0.0d0
        gaus(2,4) = + 0.0d0
        gaus(2,5) = + 1.0d0
        gaus(2,6) = + 0.0d0
        gaus(3,1) = - 1.0d0
        gaus(3,2) = - 1.0d0
        gaus(3,3) = - 1.0d0
        gaus(3,4) = + 1.0d0
        gaus(3,5) = + 1.0d0
        gaus(3,6) = + 1.0d0
        weis(1) = 0.166666666666667d0
        weis(2) = 0.166666666666667d0
        weis(3) = 0.166666666666667d0
        weis(4) = 0.166666666666667d0
        weis(5) = 0.166666666666667d0
        weis(6) = 0.166666666666667d0
!
        return
        end


