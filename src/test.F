!********************************************************************
!		modified trans_data2 subroutine 
!********************************************************************
	subroutine trans_data(ien,x1,umag_trans,vort_trans,
     $		   he_trans,ntrans)

	implicit none
	include "global.h"
	integer ien(nen,nec),rng(neface,nec) !,npart(20)
	real*8  x1(nsd,nnl)
	integer mtrans
	parameter (mtrans = 51)
	integer ntran,ntrans,ia

	real*8  umag_trans(nnl),he_trans(nnl)
	real*8  vort_trans(nnl)
	integer itrans,n,inen,inode,inode1,isd,idof,ig
	real*8 xt(2,mtrans)
	real*8 dist,MinDistance
	integer nn2d,ne2d,ie,m,itmp,i,j,ieln
	real*8 tmp1,tmp2,tmp3

	ntran = 1
	xt(1,1) = 0.2d0
	xt(2,1) = 0.0d0


	nn2d = nn/nslices
	ne2d = ne/nslices

!	ELEMENT LOOP STARTS HERE 
	do 200 ie=1,nec

	do j=1,nen
	  inode = ien(j,ie)
!	    write(6,*) inode,x1(1,inode),x1(2,inode)
	  MinDistance=1e10
	  do i=1,ntran
	  tmp1 = x1(1,inode)-xt(1,i)
	  tmp2 = x1(2,inode)-xt(2,i)
	  dist = dsqrt(tmp1*tmp1+tmp2*tmp2)
	  if(dist.lt.MinDistance) then
	  MinDistance = dist
	  itrans        = i
	  end if
	  end do
!	FINDING THE SLICE NUMBER ON WHICH THE NODE LIES  
	  m=(inode/nn2d) + 1
          if(mod(inode,nn2d).eq.0) m=m-1
!	FINDING THE ELEMENT NUMBER IN SPANWISE DIRECTION  
	  ieln = (ie/ne2d) + 1
	  if(mod(ie,ne2d).eq.0) ieln=ieln-1
	  ia = m+itrans-1

	  he_trans(inode)   = he_t(ia)
	  vort_trans(inode) = vort_t(ia)
	  umag_trans(inode) = 0.0 !umag_t(ia)

	end do

200	continue
!	ELEMENT LOOP END HERE 

	return
	end
