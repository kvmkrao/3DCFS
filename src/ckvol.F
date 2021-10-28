	subroutine ckvol(ien,xl)

c**************************************************************************
c	BLOCK FOR 2D INCOMPRESSIBLE SEMI-DISCRETE GAL/SUPG/PSPG METHOD
c       TO COMPUTE RESIDUALS FOR MATRIX FREE GMRES TECHNIQUE. 
c	LOOPS HAVE BEEN UNROLLED TO ACHIEVE SPEED UP
c	WRITTEN BY VINAY KALRO 05/04/93
c****************************************************************************
	implicit none

	real*8 mk
	parameter (mk=0.333333333)
	include "global.h"

C...EXTERNAL VARIABLES
	real*8  xl(nsd,nnl)
	integer ien(nen,nec)

C...INTERNAL VARIABLES
	real*8 x1(3,nen) !,sh(0:3,nen),xr(3,3),cf(3,3)

	real*8 eff0,hk,evol
C**************************************************************************
c......ARRAYS FOR SHAPE ROUTINES
	real*8  tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp18
	real*8  tmp9,tmp10,tmp11,tmp12,tmp13,tmp14,tmp15,tmp16,tmp17
	real*8  xr11,xr12,xr13,xr21,xr22,xr23,xr31,xr32,xr33
	real*8  cf11,cf12,cf13,cf21,cf22,cf23,cf31,cf32,cf33
	integer na,ig,inl,k,isd,jsd,imat,iee,inode,level,ie,idf,i,j
	real*8  h1,h2,h3,h4,h5,h6
C**************************************************************************
	integer l1

        real*8  sh(0:9,nen),xs(nsd,nsd),det
        real*8  gaus(3,9)
        real*8  weis(9)



	do ie = 1,nec
		evol = 0.0

C...LOCALIZE  ARRAYS
	        do inl = 1,nen
		k = ien(inl,ie)
        		do isd = 1,nsd
                       x1(isd,inl)  = xl(isd,k)
                	enddo
        	enddo

	do ig=1,nquad3

C...SHAPE FUNCTION SUBROUTINE
!1		       if(nen.eq.8) then
!1		       		include "sf3d8n.h"
!!		       elseif(nen.eq.15.or.nen.eq.6)then
				include "sf3d15n.h"
!!			else
!!     				include "sf3d4n.h"
!!		       endif
	
C...WEIGHTS
		det=abs(det)
		eff0 = wq3(ig)*det
		evol = evol +  eff0

	enddo ! quadrature loop
		hk = abs(evol)**(1.0/3.0)
!		if(nen.eq.18.or.nen.eq.6)then
		   tmp1 = x1(1,2)-x1(1,1)
                   tmp2 = x1(2,2)-x1(2,1)
                   tmp3 = x1(3,2)-x1(3,1)
                   h1 = sqrt( tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3 )
                   tmp1 = x1(1,3)-x1(1,1)
                   tmp2 = x1(2,3)-x1(2,1)
                   tmp3 = x1(3,3)-x1(3,1)
                   h2 = sqrt( tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3 )
                   tmp1 = x1(1,2)-x1(1,3)
                   tmp2 = x1(2,2)-x1(2,3)
                   tmp3 = x1(3,2)-x1(3,3)
                   h3 = sqrt( tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3 )
	 	   hk=0.5*min(h1,h2,h3)
!		end if
			
		he(ie) = hk
!	murali
		write(6,*) he(ie)
!	murali
		volume(ie) = evol	


	enddo ! element loop 

	return
	end

