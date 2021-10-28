	subroutine getkey(key)
	implicit none
	include "global.h"
	character*255 key

C...BCAST DONE IN C
	call rdchar(key)
	return
	end

	subroutine getint(key)
	implicit none
	include "global.h"
	integer  key

	if(my_pe.eq.0)call rdint(key)
#if defined(UNI)
#else
	call MPI_BCAST(key,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

	return
	end

	subroutine getreal(key)
	implicit none
	include "global.h"
	real*8 key

	if(my_pe.eq.0) call rdreal(key)
#if defined(UNI)
#else
	call MPI_BCAST(key,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif

	return
	end
