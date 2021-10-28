	subroutine tstamp(time, itime)
	include "mpf.h"
	real* 8 time(*)
	integer itime
	real*8 tbegin, tend
	save tbegin, tend

#if defined(UNI)
	tbegin=0.0; tend=0.0;
#else

	call MPI_BARRIER (MPI_COMM_WORLD,ierr)
	tend = MPI_WTIME()
#endif

	if (itime.gt.0) time(itime) = time(itime) + tend - tbegin
	tbegin = tend

	return
	end

c***************************************************************************
