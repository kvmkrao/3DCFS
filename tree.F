	subroutine treesetup(numprocs,my_id)
	implicit none
	include "mpf.h"
	include "tree.h"
	include "fmalloc.h"
	integer numprocs,my_id
	integer j,ip,np,il,tkc,tkp,k

	k = 1
	do while(k.le.numprocs)
		k = k*2
	enddo
	k = k/2

C...COMPUTE SIZE OF TREE
	j = ((numprocs-1)/k + 1)*k

	tn = j

	my_pe = my_id

C...COMPUTE LEVELS IN TREE 
	tl = 0
	do while(j.ne.1)
		j = j/2
		tl = tl + 1
	enddo

	if(tl.gt.0) then
	tsendptr = malloc(tl*isize)
	trecvptr = malloc(tl*isize)
	endif

	do j = 1,tl
		tsend(j) =-1
		trecv(j) =-1
	enddo

	np = 1
	tkp = tn
	tkc = tkp/2

	do il = 1,tl
		do ip = 1,np
			isrc = (ip-1)*(tkp)
			ides =  isrc + tkc
			if(ides.lt.numprocs) then
				if(isrc.eq.my_pe) tsend(il) = 1
				if(ides.eq.my_pe) trecv(il) = 1
			endif
		enddo
		np  = np *2
		tkp = tkp/2
		tkc = tkc/2
	enddo


c	if(my_pe.eq.0) write(6,*)'TREE LEVELS',tl

	return

	end


	subroutine vk_bcast(value,num)
C...USE TREE STRUCTURE TO BROADCAST A VALUE
	implicit none
	include "mpf.h"
	include "tree.h"
	integer num
	real*8 value(num)
	integer il,tk

	tk = tn/2

C...LOOP OVER LEVELS
#if defined(UNI)
#else
	do il = 1,tl
C...LEFT CHILD (SAME AS PARENT)
C...RIGHT CHILD 
		isrc = my_pe - tk
		ides = my_pe + tk
		if(tsend(il).ne.-1)
     $		        call MPI_ISEND(value,num,MPI_DOUBLE_PRECISION,ides,
     $				      100,MPI_COMM_WORLD,request,ierr)
		if(trecv(il).ne.-1)
     $		        call MPI_RECV(value,num,MPI_DOUBLE_PRECISION,isrc,
     $				      100,MPI_COMM_WORLD,istat,ierr)
		tk = tk/2
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	enddo
#endif


	return
	end

		
C...USE TREE STRUCTURE TO REDUCE A VALUE
	subroutine vk_rduce(value,num)
	implicit none
	include "mpf.h"
	include "tree.h"
	include "fmalloc.h"
	integer num
	real*8 value(num),tmp(100)
C	pointer (tmpptr,tmp)
	integer il,tk,j

C	tmpptr = malloc(fsize*num)
	tk = 1

C...LOOP OVER LEVELS
#if defined(UNI)
#else
	do il = tl,1,-1
C...LEFT CHILD (SAME AS PARENT)
		do j = 1,num
			tmp(j) = value(j)
		enddo
C...RIGHT CHILD 
		ides = my_pe - tk
		isrc = my_pe + tk
		if(trecv(il).ne.-1) then
      		       call MPI_ISEND(value,num,MPI_DOUBLE_PRECISION,ides,
     $				      100,MPI_COMM_WORLD,request,ierr)
		endif
		if(tsend(il).ne.-1) then
      		        call MPI_RECV(value,num,MPI_DOUBLE_PRECISION,isrc,
     $				      100,MPI_COMM_WORLD,istat,ierr)
			do j = 1,num
			value(j) = value(j) + tmp(j)
			enddo
		endif

		tk = tk*2

		call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	enddo
#endif

C	call free(tmpptr)

	return
	end

		
			
C...USE TREE STRUCTURE TO FIND MAX VALUE IN AN ARRAY 
	subroutine vk_max(fmax)
	implicit none
	include "mpf.h"
	include "tree.h"
	include "fmalloc.h"
	real*8  fmax,tmp
	integer il,tk,j

	tk = 1

C...LOOP OVER LEVELS
#if defined(UNI)
#else
	do il = tl,1,-1
C...LEFT CHILD (SAME AS PARENT)
		tmp = fmax
C...RIGHT CHILD 
		ides = my_pe - tk
		isrc = my_pe + tk
		if(trecv(il).ne.-1) then
      		       call MPI_ISEND(fmax,1,MPI_DOUBLE_PRECISION,ides,
     $				      100,MPI_COMM_WORLD,request,ierr)
		endif
		if(tsend(il).ne.-1) then
      		        call MPI_RECV(fmax,1,MPI_DOUBLE_PRECISION,isrc,
     $				      100,MPI_COMM_WORLD,istat,ierr)
			if(tmp.gt.fmax) fmax = tmp 
		endif

		tk = tk*2

		call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	enddo

        call MPI_BCAST(fmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif

	return
	end

	
C...USE TREE STRUCTURE TO FIND MIN VALUE IN AN ARRAY 
	subroutine vk_min(fmin)
	implicit none
	include "mpf.h"
	include "tree.h"
	include "fmalloc.h"
	integer n
	real*8  fmin,tmp
	integer il,tk,j

	tk = 1

C...LOOP OVER LEVELS
#if defined(UNI)
#else
	do il = tl,1,-1
C...LEFT CHILD (SAME AS PARENT)
		tmp = fmin
C...RIGHT CHILD 
		ides = my_pe - tk
		isrc = my_pe + tk
		if(trecv(il).ne.-1) then
      		       call MPI_ISEND(fmin,1,MPI_DOUBLE_PRECISION,ides,
     $				      100,MPI_COMM_WORLD,request,ierr)
		endif
		if(tsend(il).ne.-1) then
      		        call MPI_RECV(fmin,1,MPI_DOUBLE_PRECISION,isrc,
     $				      100,MPI_COMM_WORLD,istat,ierr)
			if(tmp.lt.fmin)fmin = tmp 
		endif

		tk = tk*2

		call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	enddo

        call MPI_BCAST(fmin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif


	return
	end
