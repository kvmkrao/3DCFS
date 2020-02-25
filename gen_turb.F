CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine genic_turb(d_t)
	implicit none
        include "global.h"
        real* 8 d_t(ndf_t,nnc)
	integer idf,in
!        real* 8 xi,yi,rad,theta,cost,sint,velr,velt,romega
	real*8 iniv

	iniv = 1.0e-10
	do in = 1,nnc
		 do idf = 1,ndf_t
	         d_t(idf,in)= iniv !ic_t(1,nrng+1)
		 end do
	end do

        return
        end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine genf_turb(ien,f,rng)
!	SET BOUNDARY CONDITION ARRAY 
!	MAKES USE OF A GLOBAL TEMP ARRAY
	implicit none
	include "global.h"
	include "fmalloc.h"

	integer rng(neface,nec),ien(nen,nec)
	integer btype,inface
	real*8  fl(ndf_t,nnl),f(ndf_t,nnc),rt,lt
	pointer (flptr,fl)
	integer inode,iface,ie,irng,jnode,iq,knode,i
	integer inl,jnl,idf,isrc,idst,iwho,itag,ilen,iproc
	integer nqt
	
	flptr = malloc(ndf_t*nnl*fsize)


C...INITIALIZE
	do inl = 1,nnl
		do idf = 1,ndf_t
			fl(idf,inl) = 0.0
		enddo
	enddo

        rt = 0.0
        lt = 0.0

        do ie=1,nec
            do iface=1,neface
                irng = rng(iface,ie)
                do inface=1,nnface
                    knode=facemap(iface,inface,elem_type)
                    if(knode.ne.0) then
                        inode=ien(knode,ie)
!	write(6,*)"knode inode ie irng",knode,inode,ie,irng
                        if(irng.eq.1) then
			  fl(1,inode) = 0.0d0
			else if(irng.eq.3) then
                          fl(1,inode) = 1.0e-10
                        endif
                    endif
                enddo
            enddo
        enddo


	call scatter(fl,f,ndf_t,2)
	call free(flptr)

	return
	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	subroutine genid_sd_turb(ien,id,rng)
C	SET ID ARRAY FOR EDDY VISCOSITY
C	MAKES USE OF A GLOBAL TEMP ARRAY

	implicit none
	include "global.h"
	include "fmalloc.h"

	integer rng(neface,nec),ien(nen,nec)
	integer btype,inface
	real*8  idl(ndf_t,nnl),idm(nsd,nnl),id(ndf_t,nnc),mid(nsd,nnc)
	pointer (idlptr,idl),(idmptr,idm)
	integer inode,iface,ie,irng,jnode,knode,iq
	integer inl,jnl,idf,isrc,idst,iwho,itag,ilen,iproc
	integer nqt,nqt_t
	integer p_pe,n_pe

	idlptr = malloc(ndf_t*nnl*fsize)

	nq_t = 0
	nqt_t = 0

	

C...INITIALIZE
	call fclear(idl,ndf_t*nnl)
	call fclear(id,ndf_t*nnc)

        do ie=1,nec
            do iface=1,neface
                irng=rng(iface,ie)
                if((irng.ge.1).and.(irng.le.nrng)) then
                    do inode=1,nnface
                        knode=facemap(iface,inode,elem_type)
                        if(knode.ne.0) then
                            jnode=ien(knode,ie)
!test murali
!                           if(irng.eq.5) then
!                       if(my_pe.eq.0) write(6,*)knode,jnode,ie,irng
!                           end if
!test murali
                            do idf=1,ndf_t
                                if(bc_t(idf,irng).eq.1) then
                                    idl(idf,jnode)=1.0
                                endif
                            enddo
                        endif
                    enddo
                endif
            enddo
        enddo


	call scatter(idl,id,ndf_t,1)
	call free(idlptr)

!	write(6,*) "ndf_t", ndf_t
	do inl = 1,nnc
		do idf = 1,ndf_t
			if(id(idf,inl).eq.0.0)  then
			id(idf,inl)=1.0
			nq_t     = nq_t + 1
!			write(6,*) nq_t
			else
			id(idf,inl) = 0.0
			endif
		enddo
	enddo
!	write(6,*) "nq_t", nq_t


#if defined(UNI)
#else
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

C...FIND TOTAL NUMBER OF EQUATIONS ON FIRST PROCESSOR AND BROADCAST TO REST
#if defined(UNI)
	nqt_t=nq_t; nq_t = nqt_t;
#else

	call MPI_REDUCE(nq_t,nqt_t,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
	if(my_pe.eq.0) nq_t = nqt_t
	call MPI_BCAST(nq_t,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
!	write(6,*)nq_t
	return
	end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine setd_turb(dn1_t,fn_t,id_t)
	implicit none
	include "global.h"
	integer inc,idf

	real* 8  dn1_t(ndf_t,nnc), fn_t(ndf_t,nnc), id_t(ndf_t,nnc)

	do inc=1,nnc
		do idf=1,ndf_t
		if (id_t(idf,inc).eq.0.0) dn1_t(idf,inc) = fn_t(idf,inc)         
		end do
	end do

	return
	end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

