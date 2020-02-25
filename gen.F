!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine genid_sd(ien,id,rng,mid)
!       SET ID ARRAY 
!       MAKES USE OF A GLOBAL TEMP ARRAY

        implicit none
        include "global.h"
        include "fmalloc.h"

        integer rng(neface,nec),ien(nen,nec)
        integer btype,inen
        real*8  idl(ndf,nnl),idm(nsd,nnl),id(ndf,nnc),mid(nsd,nnc)
        pointer (idlptr,idl),(idmptr,idm)
        integer inode,iface,ie,irng,jnode,knode,iq
        integer inl,jnl,idf,isrc,idst,iwho,itag,ilen,iproc
        integer nqt
        integer p_pe,n_pe

        idlptr = malloc(ndf*nnl*fsize)

        nq = 0
        nqt = 0

!...INITIALIZE
        call fclear(idl,ndf*nnl)
        call fclear(id,ndf*nnc)

!debug
!       write(6,*)'elem_type: ',elem_type, 'nnface: ',nnface,' neface: ',neface
!debug

!.......ie = element no
!.......iface = face no on that element
!.......irng = bnd no on which the face lies
!.......inode = node no on that face(1 to 3 for triangle and 1 to 4 for
!.......rectangle)
!.......knode = local node no for the wedge (1 to 6)
!.......jnode = global node no
!.......bc(idf,irng) tells whether the boundary condition is specified
!.......or not
        do ie=1,nec
            do iface=1,neface
                irng=rng(iface,ie)
                if((irng.ge.1).and.(irng.le.nrng)) then
                    do inode=1,nnface
                        knode=facemap(iface,inode,elem_type)
                        if(knode.ne.0) then
                            jnode=ien(knode,ie)
!test murali
!			    if(irng.eq.5) then
!     			if(my_pe.eq.0) write(6,*)knode,jnode,ie,irng	
!			    end if
!test murali
                            do idf=1,ndf
                                if(bc(idf,irng).eq.1) then
                                    idl(idf,jnode)=1.0
                                endif
                            enddo
                        endif
                    enddo
                endif
            enddo
        enddo


        call tstamp(time,3)
        call scatter(idl,id,ndf,1)
        call tstamp(time,5)
!FIXME..decide if it needs to be freed or not

C...FIX PRESSURE BOUNDARY CONDITION
        if(pnode.gt.0) then
            p_pe = (pnode-1)/maxnnc
            n_pe = pnode - p_pe*maxnnc
            if(my_pe.eq.p_pe) then
                    id(4,n_pe) = 1.0
!                       write(6,*) 'FIXING PRESSURE ON ',my_pe
            endif
        endif
        

!.......following is to get lm array.
        do inl = 1,nnc
            do idf = 1,ndf
                if(id(idf,inl).eq.0.0) then
                    nq     = nq + 1
                    id(idf,inl)=1.0
                else
                id(idf,inl) = 0.0
                endif
            enddo
        enddo

        call gather(idl,id,ndf)


!!        do ie=1,nec
!!            do inen=1,nen
!!                write(140+my_pe,*) (lm(idf,inen,ie),idf=1,ndf)
!!            enddo
!!        enddo

        call free(idlptr)

#if defined(UNI)
#else
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

C...FIND TOTAL NUMBER OF EQUATIONS ON FIRST PROCESSOR AND BROADCAST TO REST
#if defined(UNI)
        nqt=nq; nq = nqt;
#else
        call MPI_REDUCE(nq,nqt,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(my_pe.eq.0) nq = nqt

        call MPI_BCAST(nq,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

#endif

        if(meshmv) then

            mnq = 0
            nqt = 0

            call fclear(mid,mnqm)
            idmptr = malloc(nsd*nnl*fsize)
            call fclear(idm,nsd*nnl)

            do ie = 1,nec
                do iface =1,neface
                    do irng = 1,nrng
                        if(rng(iface,ie).eq.irng) then

                            do inode = 1,nnface

                                knode=facemap(iface,inode,elem_type)
                                if(knode.ne.0)then
                                    jnode = ien(knode,ie)

                                    do idf = 1,nsd
                                        if(bcm(idf,irng).eq.1) idm(idf,jnode)=1.0
                                    enddo
                                end if
                            enddo

                            goto 30

                        endif

                    enddo
30                  continue
                enddo
            enddo

            call tstamp(time,3)
            call scatter(idm,mid,nsd,1)
            call tstamp(time,5)
            call free(idmptr)

C...ACCOUNT FOR BOUNDARY LAYER NODES
            if(turb_model.eq.1) then
                do inl = 1,nnc
                    if(hwall(inl).lt.blyr) then
                        do idf = 1,nsd
                            mid(idf,inl) = 1.0
                        enddo
                    endif
                enddo
            endif


            do inl = 1,nnc
                do idf = 1,nsd
                    if(mid(idf,inl).eq.0.0)  then
                        mid(idf,inl)=1.0
                        mnq     = mnq + 1
                    else
                        mid(idf,inl) = 0.0
                    endif
                enddo
            enddo

#if defined(UNI)
#else
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C...FIND TOTAL NUMBER OF EQUATIONS ON FIRST PROCESSOR AND BROADCAST TO REST
            call MPI_REDUCE(mnq,nqt,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            if(my_pe.eq.0) mnq = nqt
            call MPI_BCAST(mnq,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

        endif

        return
        end





CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine genid(ien,id,rng,mid)
C.......SPACE-TIME
C       SET ID ARRAY 
C       MAKES USE OF A GLOBAL TEMP ARRAY

        implicit none
        include "global.h"
        include "fmalloc.h"
       
        integer rng(neface,nec),ien(nen,nec)
        integer btype
        real*8  idl(ndf,nnl),idm(nsd,nnl),id(ndf,nnc),mid(nsd,nnc)
        pointer (idlptr,idl),(idmptr,idm)
        integer inode,iface,ie,irng,jnode,knode,iq
        integer inl,jnl,idf,isrc,idst,iwho,itag,ilen,iproc
        integer nqt
        integer p_pe,n_pe

        idlptr = malloc(ndf*nnl*fsize)

        nq = 0
        nqt = 0

C...INITIALIZE
        call fclear(idl,ndf*nnl)
        call fclear(id,ndf*nnc)

        do ie = 1,nec
            do iface = 1,neface
                do irng = 1,nrng
                    if(rng(iface,ie).eq.irng) then
                        do inode = 1,nnface
                            knode=facemap(iface,inode,elem_type)
                            if(knode.ne.0)then 
                                jnode = ien(knode,ie)
                                do idf = 1,ndf
                                    if(bc(idf,irng).eq.1) idl(idf,jnode)=1.0
                                enddo
                            end if
                        enddo

                        goto 20

                    endif

                enddo
 20             continue
            enddo
        enddo

        call tstamp(time,3)
        call scatter(idl,id,ndf,1)
        call tstamp(time,5)
        call free(idlptr)

C...FIX PRESSURE BOUNDARY CONDITION
        if(pnode.gt.0) then
            p_pe = (pnode-1)/maxnnc
            n_pe = pnode - p_pe*maxnnc
            if(my_pe.eq.p_pe) then
                id(4,n_pe) = 1.0
                write(6,*) 'FIXING PRESSURE ON ',my_pe
            endif
        endif
        

        do inl = 1,nnc
            do idf = 1,ndf
                if(id(idf,inl).eq.0.0)  then
                    id(idf,inl)=1.0
                    nq     = nq + 1
                else
                    id(idf,inl) = 0.0
                endif
            enddo
        enddo

#if defined(UNI)
#else
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

C...FIND TOTAL NUMBER OF EQUATIONS ON FIRST PROCESSOR AND BROADCAST TO REST
#if defined(UNI)
        nqt=nq; nq = 2*nqt;
#else
        call MPI_REDUCE(nq,nqt,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(my_pe.eq.0) nq = 2*nqt
        call MPI_BCAST(nq,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

        if(meshmv) then

        mnq = 0
        nqt = 0

        call fclear(mid,mnqm)
        idmptr = malloc(nsd*nnl*fsize)
        call fclear(idm,nsd*nnl)

        do ie = 1,nec
            do iface =1,neface
                do irng = 1,nrng
                    if(rng(iface,ie).eq.irng) then

                        do inode = 1,nnface

                            knode=facemap(iface,inode,elem_type)
                            if(knode.ne.0)then
                                jnode = ien(knode,ie)
                                do idf = 1,nsd
                                    if(bcm(idf,irng).eq.1) idm(idf,jnode)=1.0
                                enddo
                            end if

                        enddo

                        goto 30

                    endif

                enddo
 30             continue
            enddo
        enddo

        call tstamp(time,3)
        call scatter(idm,mid,nsd,1)
        call tstamp(time,5)
        call free(idmptr)

C...ACCOUNT FOR BOUNDARY LAYER NODES
        if(turb_model.eq.1) then
            do inl = 1,nnc
                if(hwall(inl).lt.blyr) then
                    do idf = 1,nsd
                        mid(idf,inl) = 1.0
                    enddo
                endif
            enddo
        endif


        do inl = 1,nnc
            do idf = 1,nsd
                if(mid(idf,inl).eq.0.0)  then
                    mid(idf,inl)=1.0
                    mnq     = mnq + 1
                else
                    mid(idf,inl) = 0.0
                endif
            enddo
        enddo

#if defined(UNI)
#else
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C...FIND TOTAL NUMBER OF EQUATIONS ON FIRST PROCESSOR AND BROADCAST TO REST
        call MPI_REDUCE(mnq,nqt,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(my_pe.eq.0) mnq = nqt
        call MPI_BCAST(mnq,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

        endif


        return
        end



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine genf(ien,f,rng,ang)
!       SET BOUNDARY CONDITION ARRAY 
!       MAKES USE OF A GLOBAL TEMP ARRAY
        implicit none
        include "global.h"
        include "fmalloc.h"

        real*8  fl(ndf,nnl),f(ndf,nnc),rt,lt
!.......Just a place holder, dimensions are garbage.
        real*8  ang(nnc)          
!.......get rid of the following, won't be needed.
        integer ibnd
        real*8  cv,gas_c,ub,ubc
        integer rng(neface,nec),ien(nen,nec)
        integer btype
        pointer (flptr,fl)
        integer inode,iface,ie,irng,jnode,iq,knode,i,k,inface
        integer inl,jnl,idf,isrc,idst,iwho,itag,ilen,iproc
        integer nqt
        flptr = malloc(ndf*nnl*fsize)

        gas_c    = 1.0d0/gamma
        cv       = gas_c/(gamma-1.0d0)


C...INITIALIZE
        do inl = 1,nnl
            do idf = 1,ndf
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
			  fl(2,inode) = 0.0d0
                          fl(3,inode) = 0.0d0
                          fl(4,inode) = 0.0d0
			if(irng.eq.2) then
                          fl(3,inode) = 0.0d0
                          fl(4,inode) = 0.0d0
                        else if (irng.eq.3) then
		          fl(1,inode) = 1.0d0
		          fl(2,inode) = rmach
		          fl(3,inode) = 0.0d0
			  fl(4,inode) = 0.0d0
		          fl(5,inode) = 1.0d0*(cv*1.0d0+0.5d0*(rmach**2)) != 6.28571436 for M=3
                        else if (irng.eq.4) then
                          fl(4,inode) = 0.0d0
			else if(irng.eq.5) then
			  fl(3,inode) = 0.0d0
                          fl(4,inode) = 0.0d0
			else
			     fl(4,inode) = 0.0d0 
			end if
                        endif
                    endif
                enddo
            enddo
        enddo

!.......application of bleed
!.......somehow calculate ub based on percentage mass removal
!	9% bleed is applying only on the upstream of the throat

!.......application of bleed ends.

!!!        do ie = 1,nec
!!!             do iface =1,neface
!!!                  do irng = 1,nrng
!!!                      if(rng(iface,ie).eq.irng) then
!!!
!!!                          do inode = 1,nnface
!!!
!!!                              knode=facemap(iface,inode,elem_type)
!!!                              if(knode.ne.0)then
!!!                                  jnode = ien(knode,ie)
!!!
!!!                                  do idf = 1,ndf
!!!                                  if(bc(idf,irng).eq.1) then
!!!                                      fl(idf,jnode)=ic(idf,irng)
!!!                                  endif
!!!                                  enddo
!!!                              end if
!!!                          enddo
!!!
!!!                          goto 20
!!!
!!!                      endif
!!!                  enddo
!!! 20               continue
!!!             enddo
!!!        enddo

        call tstamp(time,3)
        call scatter(fl,f,ndf,2)
        call tstamp(time,5)
        call free(flptr)


        return
        end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine bndryinfo(ien,f,rng)
C       SET BOUNDARY CONDITION ARRAY
C       MAKES USE OF A GLOBAL TEMP ARRAY
        implicit none
        include "global.h"
        include "fmalloc.h"
 
        integer rng(neface,nec),ien(nen,nec)
        integer btype
        real*8  fl(ndf,nnl),f(ndf,nnc),dminmax(2,5)
        pointer (flptr,fl)
        integer inode,iface,ie,irng,jnode,iq,knode,i
        integer inl,jnl,idf,isrc,idst,iwho,itag,ilen,iproc
        integer nqt
 
        flptr = malloc(ndf*nnl*fsize)
 
        call gather(fl,f,ndf)
        do i=1,ndf
        dminmax(1,i)=+1000000.0
        dminmax(2,i)=-1000000.0
        enddo
 
        do ie = 1,nec
                 do iface =1,neface
                        do irng = 1,nrng
                            if(rng(iface,ie).eq.1) then
 
                                do inode = 1,nnface 
 
                                        knode=facemap(iface,inode,elem_type)
                                        if(knode.ne.0)then
                                        jnode = ien(knode,ie)
 
                                        do i=1,ndf
                                        dminmax(1,i)=min(dminmax(1,i),fl(i,jnode))
                                        dminmax(2,i)=max(dminmax(2,i),fl(i,jnode))
                                        enddo
                                        end if
                                enddo
 
                                goto 20
 
                            endif
                        enddo
 20     continue
                 enddo
        enddo
        if(my_pe.eq.0) then
                  write(6,'("dmin(PARA)=",4e16.8)') dminmax(1,1),dminmax(1,2),
     .                  dminmax(1,3),dminmax(1,4)
                  write(6,'("dmax(PARA)=",4e16.8)') dminmax(2,1),dminmax(2,2),
     .                  dminmax(2,3),dminmax(2,4)
        endif
 
        call free(flptr)
 
 
        return
        end    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine genx(ien,f,rng)
C       SET BOUNDARY CONDITION ARRAY
C       MAKES USE OF A GLOBAL TEMP ARRAY
C       converted from genf by vk@rice.edu.Apr08.2002@2111 Holly Hall
        implicit none
        include "global.h"
        include "fmalloc.h"

        integer rng(neface,nec),ien(nen,nec)
        integer btype
        real*8  fl(nsd,nnl),f(nsd,nnc),rt,lt
        pointer (flptr,fl)
        integer inode,iface,ie,irng,jnode,iq,knode,i
        integer inl,jnl,idf,isrc,idst,iwho,itag,ilen,iproc
        integer nqt

        flptr = malloc(nsd*nnl*fsize)


C...INITIALIZE
        do inl = 1,nnl
                do idf = 1,nsd
                        fl(idf,inl) = 0.0
                enddo
        enddo

        rt = 0.0
        lt = 0.0
        do ie = 1,nec
                 do iface =1,neface
                        do irng = 1,nrng
                            if(rng(iface,ie).eq.sprng) then

                                do inode = 1,nnface

                                        knode=facemap(iface,inode,elem_type)
                                        jnode = ien(knode,ie)

                                        do idf = 1,nsd
                                        if(bc(idf,irng).eq.1) then
                                        fl(idf,jnode)=dt*ic(idf,irng)
                                        endif
                                        enddo

                                enddo

                                goto 20

                            endif
                        enddo
 20     continue
                 enddo
        enddo

        call tstamp(time,3)
        call scatter(fl,f,nsd,2)
        call tstamp(time,5)
        call free(flptr)


        return
        end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine genic(x,d)
        implicit none
        include "global.h"
        real* 8 d(ndf,nnc),x(nsd,nnc)
        integer idf,in,i,jnode,k
        real* 8 xi,yi,rad,theta,cost,sint,velr,velt,romega
        real*8 ufree, vfree, wfree
        real*8 gas_c, cv, rhofree, tempfree
        real*8 init_d(ndf)     
        integer ie,iface,irng,inface,knode,inode
        integer rng(neface,nec),ien(nen,nec)


        ufree    = rmach
        vfree    = 0.0d0
        wfree    = 0.0d0
        gas_c    = 1.0d0/gamma
        cv       = gas_c/(gamma-1.0d0)
        rhofree  = 1.0d0
        tempfree = 1.0d0

!.......initial values of various dofs.
        init_d(1) = rhofree
        init_d(2) = rhofree*ufree
        init_d(3) = rhofree*vfree
        init_d(4) = rhofree*wfree
        init_d(5) = rhofree*(cv*tempfree + 
     $      0.5*(ufree*ufree + vfree*vfree + wfree*wfree))

        do in = 1,nnc
            do idf = 1,ndf
!!!                d(idf,in)  = ic(idf,nrng+1)
                d(idf,in) = init_d(idf)
            enddo
        enddo

        go to 200
        do ie=1,nec
            do iface=1,neface
                irng = rng(iface,ie)
                do inface=1,nnface
                    knode=facemap(iface,inface,elem_type)
                    if(knode.ne.0) then
                        inode=ien(knode,ie)
                        if(irng.eq.3) then
                        d(1,inode) = init_d(1)
                        d(2,inode) = rmach !init_d(2)
                        d(3,inode) = init_d(3)
                        d(4,inode) = init_d(4)
                        d(5,inode) = init_d(5)
                        end if
                     end if
                end do
            end do
        end do
        write(6,*)"HI... I am in GENIC subroutine"
200     continue


c       do k = 1,numprocs
c c         write(91+my_pe,11) ((d(i,jnode),i=1,2),jnode=1,nnc)
c           if(my_pe.eq.0.and.jnode.eq.1)write(92,*) tt
c 11        format(2e16.9)
c       enddo


        return
        end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine setd(d,f,id)
        implicit none
        include "global.h"
        integer inc,idf

        real* 8  d(ndf,nnc), f(ndf,nnc), id(ndf,nnc)

        do inc=1,nnc
                do idf=1,ndf
                        if (id(idf,inc).eq.0.0) d(idf,inc) = f(idf,inc)         
                end do
        end do

        return
        end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine updatebc(ien,rng,dn1,dn2,mvsrf)
        implicit none
        include "global.h"
        include "fmalloc.h"

        integer rng(neface,nec),ien(nen,nec)
        real*8  dn1(ndf,nnc),dn2(ndf,nnc),mvsrf(nnc)
        real*8  dl(nnl)
        pointer (dlptr,dl)
        integer inode,iface,ie,jnode,knode,idf,i,k


        if(its.eq.1) then

C...FORM mvsrf
                dlptr = malloc(fsize*nnl)
                call fclear(dl,nnl)

                do ie = 1,nec
                         do iface =1,neface
                                    if(rng(iface,ie).eq.force_surface) then
                                        do inode = 1,nnface
                                        knode=facemap(iface,inode,elem_type)
                                        jnode = ien(knode,ie)
!					write(6,*)"jnode knode ie",jnode,knode,ie
                                        dl(jnode)=1.0
                                        enddo
                                     endif
                        enddo
                enddo

                call scatter(dl,mvsrf,1,1)
                call free(dlptr)

        endif

        

C...UPDATE VELOCITY ON THE DESIRED BOUNDARY
        do inode = 1,nnc
                if(mvsrf(inode).gt.0.0) then
                        do idf = 1,nsd
                                dn1(idf,inode) = dn2(idf,inode)
                        enddo
                endif
        enddo

        return
        end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine testmv(ien,rng,xn1,xn2,mvsrf)
        implicit none
        include "global.h"
        include "fmalloc.h"

        integer rng(neface,nec),ien(nen,nec)
        real*8  xn1(nsd,nnc),xn2(nsd,nnc),mvsrf(nnc)
        integer inode,iface,ie,jnode,knode,idf

        

c       do ie = 1,nec           !*********************************************
c        do inode = 1,6
c        jnode = ien(inode,ie) 
c        if(my_pe.eq.4.and.jnode.eq.1234)write(6,*)'INODE,IE In P-5',inode,ie,xn1(1,j),xn1(2,j),xn1(3,j) 
c        end do 
c        end do!**********************************

C...MOVE MESH
        do inode = 1,nnc
                if(mvsrf(inode).gt.0.0) then
                do idf = 1,nsd
                xn2(idf,inode) = xn2(idf,inode)+ ic(idf,force_surface)
                enddo
                endif
        enddo


        return
        end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine rotatebc (ien,rng,d1,x,nt1)
c       subroutine rotatebc (ien)
        implicit none
        include "global.h"
        include "fmalloc.h"
        include "ncomm.h"
        integer rng(neface,nec),ien(nen,nec)
        real *8 d1(ndf,nnl),x(nsd,nnl),ut,theta,theta_t,nt1
        integer ie,iface,inode,jnode,knode,i,k,iunit

        ut=perc*vel
        theta=perc*nt1

        do  ie = 1,nec
        do  iface =1,neface
                if(rng(iface,ie).eq.force_surface) then
                do  inode = 1,nnface
                         knode=facemap(iface,inode,elem_type)
                        jnode = ien(knode,ie)


c                        write(91+my_pe,11) jnode,(d1(i,jnode),i=1,2)
c                       if(my_pe.eq.1) write(92,11) jnode
c 11                    format(i6,2x,4e16.9)



c                       if((x(1,jnode).ge.0).and.(x(2,jnode).ge.0))
c     $                         theta_t=theta+asin(x(2,jnode))
c                       if((x(1,jnode).le.0).and.(x(2,jnode).ge.0))
c     $                         theta_t=theta-asin(x(2,jnode))+2.0*asin(1.0)
c                       if((x(1,jnode).le.0).and.(x(2,jnode).le.0))
c     $                         theta_t=theta+asin(dabs(x(2,jnode)))+2.0*asin(1.0)
c                       if((x(1,jnode).ge.0).and.(x(2,jnode).le.0))     
c     $                         theta_t=theta-asin(dabs(x(2,jnode)))+4.0*asin(1.0)

                        d1(1,jnode) = -ut*sin(theta_t)
                        d1(2,jnode) =  ut*cos(theta_t)


c                       do k = 1,numprocs
c                       write(91+my_pe,11) jnode,(d1(i,jnode),i=1,2),
c     $                         (x(i,jnode),i=1,2)     
c                      if(my_pe.eq.0.and.jnode.eq.1)write(92,*) 'nt1,theta',nt1,theta
c 11                    format(i6,2x,2e16.9,2x,2e16.9)
c                       enddo


                enddo
                endif
        enddo
        enddo

c               write(6,*)'Sine inverse',my_pe,nec
        return
        end
