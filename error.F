	subroutine error(string, param, fatal)
	include "mpf.h"

c	print a message and optional parameter on stderr, stop if error is fatal
	character string*(*)
	integer param,my_pe,ierr
	logical fatal

#if defined(UNI)
#else
	call MPI_COMM_RANK(MPI_COMM_WORLD,my_pe,ierr)
#endif


	if(my_pe.eq.0) then
		if (param.ne.-999) then
			write(unit=0, fmt='(a$)') string(1:mylen(string))
			write(unit=0, fmt=*) param
		else
			write(unit=0, fmt='(a)') string(1:mylen(string))
		end if
	endif

	if (fatal) then
		!call exit(1)
		stop
	endif

	return
	end
