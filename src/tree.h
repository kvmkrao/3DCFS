	integer tn,tl,tsend(1),trecv(1)
	integer ides,isrc,ierr,request,istat(MPI_STATUS_SIZE),my_pe
	pointer (tsendptr,tsend),(trecvptr,trecv)
	common/ctree/tn,tl,ides,isrc,ierr,request,istat,my_pe,tsendptr,trecvptr
	
