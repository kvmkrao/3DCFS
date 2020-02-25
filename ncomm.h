C...COMM
C	include "/usr/include/mpp/mpi.h"
	integer ja_rcv(1),ja_snd(1),ia_rcv(1),ia_snd(1),snd_home(1),rcv_home(1)
	integer nnbs,nnbr,nni,nne,nnrtot,snd_nbrs,rcv_nbrs
	integer jstat

	integer SSTAT( MPI_STATUS_SIZE,512),RSTAT(MPI_STATUS_SIZE,512)
	integer IRQ(512),ISQ(512)

	pointer (ja_rcvptr,ja_rcv),(ja_sndptr,ja_snd)
	pointer (ia_rcvptr,ia_rcv),(ia_sndptr,ia_snd)
	pointer(snd_homeptr,snd_home),(rcv_homeptr,rcv_home)
	common/ccomm/ja_rcvptr,ja_sndptr,snd_homeptr,rcv_homeptr,
     $		     ia_rcvptr,ia_sndptr,nnbs,nnbr,nni,nne,
     $		     nnrtot,snd_nbrs,rcv_nbrs
  


C...COMMUNICATION SETUP FOR CEBE
	integer colors(1),o2n(1),o2nloc(1),colflag
	pointer (colorsptr,colors),(o2nptr,o2n),(o2nlocptr,o2nloc)
	common/cebecomm/colorsptr,o2nptr,colflag
