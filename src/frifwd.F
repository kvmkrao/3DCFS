	subroutine printf(ien)
	real*8 ien(nsd,nnc)
	include "global.h"
	 write(61+my_pe,111)((ien(i,j),i=1,nsd),j=1,nnc)
111	format (3f25.9) !(6i8)
	return
	end

        subroutine printi(ien)
        integer ien(5,nec)
        include "global.h"
         write(81+my_pe,111)((ien(i,j),i=1,5),j=1,nec)
111     format (5i8)
        return
        end

	subroutine frif(filesb1,ienb,dim1,dim2)
	implicit none
        include "global.h"
	include "fmalloc.h"
	character*255 filesb1
	integer dim1, dim2
	integer ienb(dim1,nec),iens(dim1,dim2)
	integer offset, i, j1, j2, nnn, mmm

	pointer (iensptr, iens)
	iensptr = malloc(dim1*dim2*isize)

c	if(my_pe.eq.0) write(6,*) 'dim1,dim2',dim1,dim2

          nnn=(dim2-1)/numprocs + 1
	 mmm=nnn
	if(my_pe.eq.numprocs-1) nnn=dim2-(numprocs-1)*mmm
	offset = mmm*my_pe + 1
        open(unit=2,file=filesb1,status='unknown',
     .     access='DIRECT',RECL=isize*dim1*dim2)
	read (2,rec=1) iens
	close(2)

	if(my_pe.eq.0) write(6,*) 'read mien',offset

	j1 = 0
	do j2 = offset, offset+nec-1
	  j1 = j1+1
	  do i = 1,dim1
	   ienb(i,j1) = iens(i,j2)
	  end do
	end do

	call MPI_BARRIER (MPI_COMM_WORLD,ierr)

	call free(iensptr)

	return
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine frdf(filesb2,ienb,dim1,dim2)
        implicit none
        include "global.h"
        include "fmalloc.h"
	integer dim1, dim2, nnn, mmm
        real*8 ienb(dim1,nnc),iens(dim1,dim2)
        integer offset, i, j1, j2
	character*255 filesb2
        pointer (iensptr, iens)
        iensptr = malloc(dim1*dim2*fsize)
	
	  nnn=(dim2-1)/numprocs + 1
	  mmm=nnn
	  if(my_pe.eq.numprocs-1) nnn=dim2-(numprocs-1)*mmm
        offset = mmm*my_pe + 1

        open(unit=2,file=filesb2,status='unknown',
     .     access='DIRECT',RECL=fsize*dim1*dim2)
        read (2,rec=1) iens
        close(2)

        j1 = 0
        do j2 = offset, offset+nnc-1
          j1 = j1+1
          do i = 1,dim1
           ienb(i,j1) = iens(i,j2)
          end do
        end do

	if(my_pe.eq.0) write(6,*) 'read mxyz',offset

        call MPI_BARRIER (MPI_COMM_WORLD,ierr)

        call free(iensptr)

        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine fwdf(filesb,ienb,dim1,dim2)
	implicit none
        include "global.h"
        include "fmalloc.h"
        integer dim1, dim2, nnn, mmm
        real*8 ienb(dim1,1),iens(dim1,dim2)
        integer offset, i, j, j1, j2, nsent,nrcv,ireq,irequest
        character*255 filesb
        pointer (iensptr, iens)
        iensptr = malloc(dim1*dim2*fsize)
        
          nnn=(dim2-1)/numprocs + 1
          mmm=nnn
          if(my_pe.eq.numprocs-1) nnn=dim2-(numprocs-1)*mmm
	nsent = dim1*nnc


	do i = 1,numprocs-1
	  if(my_pe.eq.i)
     $		call MPI_ISEND (ienb(1,1),nsent,MPI_DOUBLE_PRECISION,0,
     $		100,MPI_COMM_WORLD,istat,ierr)

	 nrcv=mmm*dim1
	 if(i.eq.numprocs-1) nrcv=(dim2-(numprocs-1)*mmm)*dim1
	 offset = mmm*i+1
	 if(my_pe.eq.0) then
     	  call MPI_IRECV(iens(1,offset),nrcv,MPI_DOUBLE_PRECISION,i,
     $	  100,MPI_COMM_WORLD,ireq,ierr)
	call MPI_WAIT(ireq,irequest,ierr)
	end if

	end do

	if(my_pe.eq.0) then
	 do j = 1,nnc
	 do i = 1,dim1
	  iens(i,j)=ienb(i,j)
	 end do
	 end do

	 open(unit=2,file=filesb,status='unknown',
     .     access='DIRECT',RECL=fsize*dim1*dim2)
         write (2,rec=1) iens
         close(2)
	end if
	call free(iensptr)
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	return
	end
