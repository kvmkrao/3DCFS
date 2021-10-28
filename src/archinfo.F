	subroutine archinfo(iunit)
	implicit none
	include "global.h"
	integer iunit

	if (my_pe.eq.0) write(iunit,1002) numprocs

 1002	format(/"HP-Intel Xeon CLUSTER with ",i5," processing nodes")

	return
	end
