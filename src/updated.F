c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	updated.fcm                                                          c
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine updated(d,b,ind,in,nl)
	implicit none
	integer ind,in,nl
	real*8 d(ind,in),b(ind,in)
	integer inl,idf
	include "fmalloc.h"
	include "global.h"
	
	do inl = 1,in
		do idf = 1,ind
			d(idf,inl)=d(idf,inl)+b(idf,inl)
		enddo
	enddo

c	if (my_pe.eq.0) print * , 'nnc  nnl:',nnc,nnl
c	if (my_pe.eq.1) print * , 'nnc nnl:',in,nl
c	if (my_pe.eq.2) print * , 'nnc nnl:',in,nl
c	if (my_pe.eq.3) print * , 'nnc nnl:',in,nl


	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine updated_vis(d,b,ind,in,nl)
        implicit none
        integer ind,in,nl
        real*8 d(ind,in),b(ind,in)
        integer inl,idf
        include "fmalloc.h"
        include "global.h"

	do inl = 1,in
		do idf = 1,ind
			d(idf,inl)=b(idf,inl)
		enddo
	enddo

	return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine updated_sum(d_sum,dn1,ind,in)
	implicit none
	integer inl,idf,ind,in
	real *8 d_sum(ind,in),dn1(ind,in)
	include "fmalloc.h"
	include "global.h"

	do inl = 1,in
		do idf = 1,ind
			d_sum(idf,inl)=d_sum(idf,inl)+dn1(idf,inl)
		end do
	end do
	
	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine updated_pro(d_pro,dn1,ind,in)
	implicit none
	integer ind,in,inl
	real *8 d_pro(6,in),dn1(ind,in)
c	include "mpi.h"
	include "global.h"
c        include "ncomm.h"
        include "fmalloc.h"

	do inl = 1,in
		d_pro(1,inl)=d_pro(1,inl)+dn1(1,inl)*dn1(1,inl)
		d_pro(2,inl)=d_pro(2,inl)+dn1(1,inl)*dn1(2,inl)
		d_pro(3,inl)=d_pro(3,inl)+dn1(1,inl)*dn1(3,inl)
		d_pro(4,inl)=d_pro(4,inl)+dn1(2,inl)*dn1(2,inl)
		d_pro(5,inl)=d_pro(5,inl)+dn1(2,inl)*dn1(3,inl)
		d_pro(6,inl)=d_pro(6,inl)+dn1(3,inl)*dn1(3,inl)
	end do

	return
	end
