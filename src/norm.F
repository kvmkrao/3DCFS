	subroutine getnorm(v1,v2,norm,neq)
	implicit none
	include "global.h"
	include "fmalloc.h"
	integer neq
	real*8 v1(neq),v2(neq),norm,sdot,tmp


	norm = sdot (neq,v1,1,v2,1)

#if defined(UNI)
#else
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

C	if(numprocs.ge.16) then
C...USER IMPLEMENTED TREES FASTER ON T3D
C		call vk_rduce(norm,1)
C		call vk_bcast(norm,1)
C	else
	tmp = norm
#if defined(UNI)
#else
	call MPI_REDUCE(tmp,norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(norm,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif
C	endif
		
	

	return
	end
