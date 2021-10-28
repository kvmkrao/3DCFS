!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine commsetup(ien)
        implicit none
        include "global.h"
        include "ncomm.h"
        include "fmalloc.h"
        integer ien(nen,nec)
        integer hsh(1),nl2h(1),tmph(1),bnd(1),nbr(1),nbs(1)
        pointer (hshptr,hsh),(nl2hptr,nl2h)
        pointer (tmphptr,tmph),(bndptr,bnd)
        pointer (nbsptr,nbs),(nbrptr,nbr)
        integer ie,inl,inode,jnode,knode
        integer iloc,ihome,itmp1,itmp2,iptr,jptr,icol
        integer isrc,ides,iproc,i,j,k,l
        real*8 pbnds
        integer imin,imax,nn_p


!...ALLOCATE A HASH TABLE
!..FIND MINIMUM AND MAXIMUM NODE NUMBER ON PROCESSOR
        imin = nn
        imax = 1

        do ie = 1,nec
            do inl = 1,nen
                inode = ien(inl,ie)
                if(inode.lt.imin)imin=inode     
                if(inode.gt.imax)imax=inode     
            enddo
        enddo


!       do i = 1,numprocs               !************************************
!       if(my_pe.eq.i-1)
!       write(81+my_pe,23)((ien(k,l),k=1,nen),l=1,nec)
! 23    format(6i9)
!       end do
!       write(41+my_pe,*)nnc
        
!       if (my_pe.eq.4)
!     $ write(6,*)'I M IN NCOMM.F::::::::::::',ien(1,1200)

        nn_p = imax-imin + 1
!       do iproc = 1,numprocs
!       if((my_pe+1).eq.iproc.and.nn_p.gt.nn)write(6,607) my_pe,nn_p,imin,imax
!       enddo
! 607   format(4i8)
!       write(6,*)'my_pe,nn_p,nnc',my_pe,nn_p,nnc
        hshptr = malloc(isize*nn_p)
        do inode =1,nn_p
            hsh(inode) = 0
        enddo


        if (debug) call error("COMMSETUP: STAGE 0",-999,.false.)

!...COUNT TOTAL NUMBER OF P NODES(nnrtot = nnl)
        do ie = 1,nec
            do inl = 1,nen
                inode = ien(inl,ie)
                jnode = inode -imin +1
                hsh(jnode) =1
            enddo
        enddo

        nnrtot = 0
        do inode = 1,nn_p
            if(hsh(inode).eq.1) nnrtot = nnrtot+1
        enddo

!...NNL:  NUMBER OF NODES AT THE CLUSTER LEVEL
        nnl = nnrtot
!                       write(41+my_pe,*)nnc,nnl!***********************************
!                       write(6,*)'my_pe,nnl,nn,maxnnc',my_pe,nnl,nn,maxnnc


!...MAKE MAPPING FROM P TO G NODES
        ja_rcvptr = malloc(isize*nnrtot)
        nl2hptr   = malloc(isize*nnrtot)
        nbrptr    = malloc(isize*numprocs)
        nbsptr    = malloc(isize*numprocs)
        do iproc = 1,numprocs
            nbr(iproc) = 0
        enddo
        do iproc = 1,numprocs
            nbs(iproc) = 0
        enddo

!...COMPUTE NL2H ARRAY (CONTAINS  HOMES OF EACH P NODE)
!...RCVLOC ARRAY( WILL EVENTUALLY CONTAIN THE CORRESPONDING G NODE)
!...RCVLOC AND NL2H TOGETHER GIVES ME G<->P CORRESPONDENCE
!...ALSO FIND NEIGHBOURS FOR EACH PN
        nnrtot = 0
        do inode = 1,nn_p
            if(hsh(inode).eq.1) then
                jnode = inode + imin -1
                nnrtot = nnrtot+1
                ja_rcv(nnrtot) = jnode
                iproc  = (jnode-1)/maxnnc + 1
                if(nbr(iproc).eq.0.and.iproc-1.ne.my_pe) 
     $              rcv_nbrs = rcv_nbrs + 1
                nbr(iproc) = nbr(iproc) + 1
                ihome = iproc -1
                nl2h(nnrtot) = ihome
            endif
                                
        enddo


!********************************************************************
!       print * , 'TEST INFORMATION'
!       print * , 'my_pe',my_pe,'nn_p',nn_p,'jnode(1)',ja_rcv(1),'jnode(nnp)',ja_rcv(nn_p)
!*********************************************************************

!...ALLOCATE  AND FORM IA_RCV ARRAY
        ia_rcvptr = malloc(isize*(rcv_nbrs+1))
        iproc = my_pe + 1
        ia_rcv(1)  = 1 + nbr(iproc)
        nbr(iproc) = 0
        if(rcv_nbrs.gt.0) rcv_homeptr = malloc(isize*rcv_nbrs)
        
        rcv_nbrs = 0
        do iproc = 1,numprocs
            if(nbr(iproc).gt.0) then
                rcv_nbrs = rcv_nbrs + 1
                ihome = iproc-1
                rcv_home(rcv_nbrs) = ihome
                ia_rcv(rcv_nbrs+1) = ia_rcv(rcv_nbrs) + nbr(iproc)
            endif
        enddo

!                       write(41+my_pe,*)nnc,nnl!**************************O.K.****
!...FIND SIZE OF RECEIVE BUFFER AND NUMBER OF ON_PN NODES
        nnbr = ia_rcv(rcv_nbrs+1)-ia_rcv(1)
        nne = nnbr
        nni = nnrtot - nnbr
                
!...ORDER ON-PROCESSOR NODES FIRST
        jnode = 0
        do inode = 1,nnrtot
            if(nl2h(inode).eq.my_pe) then
                jnode = jnode + 1
                hsh(jnode) = inode
            endif
        enddo
!                       write(41+my_pe,*)nnc,nnl!**************************O.K. 

!...ORDER OFF-PROCESSER NODES
        do iproc = 1,numprocs
            ihome = iproc -1
            if(my_pe.ne.ihome) then
                do inode = 1,nnrtot
                    if(nl2h(inode).eq.ihome) then
                        jnode = jnode + 1
                        hsh(jnode) = inode
                    endif
                enddo
            endif
        enddo

!               write(41+my_pe,*)nnc,nnl!**************************O.K.

        call free(nl2hptr)
!                                       if(my_pe.eq.0)write(6,*)'I AM OK HERE '
!...REORDER ja_rcv ARRAYS BASED ON HASH  TABLE
        tmphptr = malloc(isize*nnrtot)
        do inode = 1,nnrtot
            tmph(inode) = ja_rcv(hsh(inode))
        enddo
        do inode = 1,nnrtot
            ja_rcv(inode)  = tmph(inode)
        enddo

        call free(tmphptr)

!...MODIFY  JA_RCV      
        do inode = 1,nnrtot
            jnode = ja_rcv(inode)
            knode = jnode -imin + 1
            hsh(knode) = inode
            ja_rcv(inode)= mod(jnode - 1,maxnnc) + 1
        enddo


!...REORDER IEN SO IT FORMS A MAP OF E<->P NODES
        do ie = 1,nec
            do inl = 1,nen
                inode = ien(inl,ie)
                jnode = inode -imin + 1
                ien(inl,ie) = hsh(jnode)
            enddo
        enddo   
        call free(hshptr)
                        
!                       write(41+my_pe,*)nnc,nnl!**************************O.K.
        

!...FIND OUT AMOUNTS OF DATA TO BE TRANSFERRED BETWEEN PROCESSORS
        if (debug) call error("COMMSETUP: STAGE 1",-999,.false.)

!...BLOCKING SENDS AND RECEIVES
#if defined(UNI)
#else
        do iproc = 1,numprocs-1
            isrc = my_pe + iproc
            if(isrc.gt.numprocs-1) isrc = isrc - numprocs
            ides = my_pe - iproc 
            if(ides.lt.0) ides = ides + numprocs
            i = ides + 1
            j = isrc + 1
            call MPI_ISEND(nbr(i),1,MPI_INTEGER,ides,100,
     $              MPI_COMM_WORLD,request,ierr)
            call MPI_RECV(nbs(j),1,MPI_INTEGER,isrc,100,
     $              MPI_COMM_WORLD,istat,ierr)

            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        enddo
#endif

!                       write(41+my_pe,*)nnc,nnl!**************************OK

!...NBR ARRAY NO LONGER REQUIRED
        call free(nbrptr)

!...REDUCE DATA FROM NBS TO FIND ACTUAL SND_NBRS 
        snd_nbrs = 0
        do iproc = 1,numprocs
            if(nbs(iproc).gt.0) then
                snd_nbrs = snd_nbrs + 1
            endif
        enddo


!...FORM A CSR STRUCTURE FOR THE  SEND BUFFERS  (DISCOUNT ON-PN NODES)
        if(snd_nbrs.ne.0)  snd_homeptr = malloc(isize*snd_nbrs)
        ia_sndptr   = malloc(isize*(snd_nbrs+1))
        ihome = 0
        do iproc = 1,numprocs   
            if(nbs(iproc).gt.0) then
                ihome = ihome + 1
                snd_home(ihome) = iproc -1
                ia_snd(ihome) = nbs(iproc)
            endif
        enddo
        
        i = ia_snd(1)
        ia_snd(1) = 1
        do iproc = 2,snd_nbrs+1
            j = ia_snd(iproc)
            ia_snd(iproc) = ia_snd(iproc-1) + i
            i = j
        enddo    

!...COMPUTE SIZE OF SEND BUFFER
        nnbs = ia_snd(snd_nbrs+1)-1

        call free(nbsptr)

!               write(41+my_pe,*)nnc,nnl!**************************OK
        
!...CONSTRUCT ja_snd ARRAY
        if (debug) call error("COMMSETUP: STAGE 2",-999,.false.)
        if(nnbs.gt.0) ja_sndptr = malloc(isize*nnbs)
!...SEND INFO FROM ja_rcv
#if defined(UNI)
#else
        do iproc = 1,rcv_nbrs
            ides = rcv_home(iproc)
            iptr = ia_rcv(iproc)
            iloc = ia_rcv(iproc+1) - ia_rcv(iproc)
            call MPI_ISEND(ja_rcv(iptr),iloc,MPI_INTEGER,ides,100,
     $          MPI_COMM_WORLD,IRQ(iproc),ierr)
        enddo
!               write(41+my_pe,*)nnc,nnl!**************************OK

!...COLLECT INFO TO ja_snd
        do iproc = 1,snd_nbrs
            isrc = snd_home(iproc)
            iptr = ia_snd(iproc)
            iloc = ia_snd(iproc+1) -ia_snd(iproc)
            call MPI_IRECV(ja_snd(iptr),iloc,MPI_INTEGER,isrc,100,
     $          MPI_COMM_WORLD,ISQ(iproc),ierr)
        enddo
!                       write(41+my_pe,*)nnc,nnl!**************************OK

        call MPI_WAITALL(rcv_nbrs,IRQ,RSTAT,ierr)       
        if(snd_nbrs.gt.0) call MPI_WAITALL(snd_nbrs,ISQ,SSTAT,ierr)

!                       write(41+my_pe,*)nnc,nnl!****************************

!       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!                       write(41+my_pe,*)nnc,nnl!***********************not ok*
#endif
        if (debug) call error("COMMSETUP:  COMPLETE",-999,.false.)
!                       write(41+my_pe,*)nnc,nnl!***********************************
        return
        end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine gather(xl,xn,n)
        implicit none
        include "global.h"
        include "ncomm.h"
        include "fmalloc.h"
        integer iproc,iptr,ides,isrc,i,j,k,iloc,jloc,inl,ie,n
        integer nsize,ist,ind
        real*8  xn(n,nnc),xl(n,nnrtot)
        real*8  bufr(n,1),bufs(n,1)
        pointer (bufrptr,bufr),(bufsptr,bufs)

!...CLEAR P ARRAYS
        do j = 1,nnrtot
            do i = 1,n
                xl(i,j) = 0.0
            enddo
        enddo
#if defined(UNI)
        do j = 1,nnl; do i=1,n; xl(i,j) = xn(i,j); enddo; enddo;
        return;
#else

!...ALLOCATE MEMORY FOR SEND  BUFFER 
        if(nnbs.gt.0) then
            bufsptr = malloc(fsize*nnbs*n)
        endif


!...OFF PN GATHER 

!...NON BLOCKING RECEIVES
        do iproc = 1,rcv_nbrs
            isrc = rcv_home(iproc)
            iptr = ia_rcv(iproc)
            nsize = n*(ia_rcv(iproc+1) - ia_rcv(iproc))
            call MPI_IRECV(xl(1,iptr),nsize,MPI_DOUBLE_PRECISION,isrc,100,
     $          MPI_COMM_WORLD,IRQ(iproc),ierr)
        enddo

!...FILL IN THE SEND BUFFERS
        ist = ia_snd(1)
        ind = ia_snd(snd_nbrs+1) -1
        jloc = 0
        do  iloc =  ist,ind
            jloc = jloc+1
            k = ja_snd(iloc)
            do j = 1,n
                bufs(j,jloc) = xn(j,k)
            enddo
        enddo

!...NON BLOCKING SENDS
        do iproc = 1,snd_nbrs
            ides = snd_home(iproc)
            nsize = n*(ia_snd(iproc+1) - ia_snd(iproc))
            iptr   = ia_snd(iproc)          
            call MPI_ISEND(bufs(1,iptr),nsize,MPI_DOUBLE_PRECISION,ides,
     $          100,MPI_COMM_WORLD,ISQ(iproc),ierr)
        enddo



!...ON PN GATHER 
        do j = 1,nni
            iloc = ja_rcv(j)
            do i = 1,n
                xl(i,j) = xn(i,iloc)
            enddo
        enddo

        call MPI_WAITALL(rcv_nbrs,IRQ,RSTAT,ierr)       

!...DEALLOCATE MEMORY FOR BUFFER 

!       if(colflag.eq.1) call permute(xl,ndf,nnl,o2n,1)

        if(snd_nbrs.gt.0)call MPI_WAITALL(snd_nbrs,ISQ,SSTAT,ierr)      
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(nnbs.gt.0) call free(bufsptr)
        return
#endif
        end
        

                

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine scatter(xl,xn,n,iopt)

!...IOPT = 1 SCATTER
!...IOPT = 2 GLOBALIZE/BNC
        implicit none
        include "global.h"
        include "ncomm.h"
        include "fmalloc.h"
        integer iproc,iptr,ides,isrc,i,j,k,n,iloc,jloc,inl,ie
        integer nsize,ist,ind,iopt
        real*8 xn(n,nnc),xl(n,nnrtot)
        real*8 bufr(n,1),bufs(n,1)
        pointer (bufrptr,bufr),(bufsptr,bufs)


!       if(colflag.eq.1) call permute(xl,ndf,nnl,o2n,2)

!...ALLOCATE MEMORY FOR BUFFER 
        if(nnbs.gt.0) then
            bufsptr = malloc(fsize*nnbs*n)
        endif


        do j = 1,nnc
            do i = 1,n
                xn(i,j) = 0.0
            enddo
        enddo
#if defined(UNI)
        do j=1,nnc; do i=1,n; xn(i,j)=xl(i,j); enddo; enddo;
        return;
#endif


                

!...OFF PN SCATTER

!...NON BLOCKING SENDS
        do iproc = 1,rcv_nbrs
            ides = rcv_home(iproc)
            iptr = ia_rcv(iproc)
            nsize = n*(ia_rcv(iproc+1) - ia_rcv(iproc))
            call MPI_ISEND(xl(1,iptr),nsize,MPI_DOUBLE_PRECISION,ides,
     $          100,MPI_COMM_WORLD,IRQ(iproc),ierr)
        enddo

!...NON BLOCKING  RECEIVES
        do iproc = 1,snd_nbrs
            isrc = snd_home(iproc)
            iptr   = ia_snd(iproc)          
            nsize = n*(ia_snd(iproc+1) - ia_snd(iproc))
            call MPI_IRECV(bufs(1,iptr),nsize,MPI_DOUBLE_PRECISION,isrc,
     $          100,MPI_COMM_WORLD,ISQ(iproc),ierr)
        enddo

!...ON PN SCATTER
        do j = 1,nni
            jloc = ja_rcv(j)
            do i = 1,n
                xn(i,jloc) = xl(i,j)
            enddo
        enddo

        if(snd_nbrs.gt.0)call MPI_WAITALL(snd_nbrs,ISQ,SSTAT,ierr)      

!...FILL OUT THE SEND BUFFERS
        ist = ia_snd(1)
        ind = ia_snd(snd_nbrs+1) -1
        jloc = 0
        goto(10,20),iopt
 10     continue
        do  iloc =  ist,ind
            jloc = jloc+1
            k = ja_snd(iloc)
            do j = 1,n
                xn(j,k) = xn(j,k) + bufs(j,jloc) 
            enddo
        enddo
 
        goto 30

 20     continue
        if (debug) call error("SCATTER:  GLOBALIZE",-999,.false.)
        do  iloc =  ist,ind
            jloc = jloc+1
            k = ja_snd(iloc)
            do j = 1,n
                if(xn(j,k).eq.0.0) xn(j,k)=bufs(j,jloc) 
            enddo
        enddo

 30     continue
!...DEALLOCATE MEMORY FOR BUFFER 
        call MPI_WAITALL(rcv_nbrs,IRQ,RSTAT,ierr)       
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(nnbs.gt.0) call free(bufsptr)

        return
        end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine commdealloc
        implicit none
        include "global.h"
        include "ncomm.h"
        integer iother

        if(snd_homeptr.ne.0) call free(snd_homeptr)
        if(rcv_homeptr.ne.0) call free(rcv_homeptr)
        if(ia_rcvptr.ne.0) call free(ia_rcvptr)
        if(ia_sndptr.ne.0) call free(ia_sndptr)
        if(ia_rcvptr.ne.0) call free(ia_rcvptr)
        if(ja_sndptr.ne.0) call free(ja_sndptr)
        if(ja_rcvptr.ne.0) call free(ja_rcvptr)

        if(iprecond.eq.3) then
            call free(colorsptr)
            call free(o2nptr)
        endif
        
        return
        end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
