C...HEADER FILE FOR INTERFACE BETWEEN VOLUME AND SURFACE MESH 
C	NODEMAP: MAPS FLUID NODES  TO SURFACE NODES
C	NORMAL : SPECIFIES WHICH SIDE OF MEMBRANE THE FLUID NODE LIES ON
C	IA_F2S : POINTS TO LOCATIONS WHERE PN_0 receives data FROM OTHERS
C	PERM   : REORDERING OF SURFACE NODES ON PN_0
C	ja_f2s : CONTAINS FLUID NODES FOR WHICH INFO IS EXCHANGED  

        include "mpf.h"

        integer request,ierr,numthreads,istat(MPI_STATUS_SIZE)
        common/ccmpi/request,ierr,numthreads,istat

	real*8 normal(1)	
	integer  nodemap(1),ia_f2s(1),ja_f2s(1),perm(1)
	integer nfs,nn_f2s
	pointer (permptr,perm),(ja_f2sptr,ja_f2s),(normalptr,normal)
	pointer (nodemapptr,nodemap),(ia_f2sptr,ia_f2s)
	common/cf2s/nodemapptr,ia_f2sptr,permptr,ja_f2sptr,normalptr,nfs,nn_f2s

