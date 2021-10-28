!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     * Main driver for the compressible fluid flow code with 
!       stabilized semi-discrete formaulation (f90)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine hypo

        implicit none

        include "global.h"
        include "fmalloc.h"
        integer ien(nen,nec),rng(neface,nec)
!!        real*8  selm(nee,nee,nec)
        real*8  idn(ndf,nnc),mvsrf(nnc),nt1
        real*8  xn1(nsd,nnc),x1(nsd,nnl)
        real*8  xn2(nsd,nnc),x2(nsd,nnl)
        real*8  xn(nsd,nnc),x_surf(5,5)
        real*8  dno(ndf,nnc),do(ndf,nnl),fn(ndf,nnc)
        real*8  dn1(ndf,nnc),d1(ndf,nnl),bg1(ndf,nnc),p1(ndf,nnl)
        real*8  dn1_glob(ndf,nnc),dno_glob(ndf,nnc)
        real*8  dn1_trans(ndf,nnc)
        real*8  d1_glob(ndf,nnl),do_glob(ndf,nnl)
        real*8  dn2(ndf,nnc),d2(ndf,nnl),bg2(ndf,nnc),p2(ndf,nnl)
        real*8  cbl(nnl)
        real*8  dn1_sum(ndf,nnc),dn1_pro(6,nnc),eddy_vis(1,nnc)
        real*8  d1_sum(ndf,nnl),d1_pro(6,nnl),eddy_vis1(1,nnl),zlen
        real*8  ang(nnc),slipang(nnc)
        
        real*8 vloc1(ndf,nnl), wg1(nqc)
        real*8 vloc2(ndf,nnl), wg2(nqc)
        real*8 id_t(ndf_t,nnc)
        real*8 dno_t(ndf_t,nnc),do_t(ndf_t,nnl),fn_t(ndf_t,nnc)
        real*8 dn1_t(ndf_t,nnc),d1_t(ndf_t,nnl)
        real*8  bg_t(ndf_t,nnc), p_t(ndf_t,nnl)
        real*8 vloc_t(ndf_t,nnl),wg_t(nqc_t)

        integer  ntrans 
        parameter(ntrans=1)
	real*8 ynode(nnc)!u have to check values @transition point
        integer inode,jnode,ielm,ie,inen
        integer rank,inface
	integer ios
          
	real*8  ynoden(nnl)
        real*8  prs_exit_ini,prs_exit_fin,time_start,time_exit
        real*8  p_down

	real*8 umag_trans(nnl), umagn_trans(nnc)
	real*8 he_trans(nnl), hen_trans(nnc)
	real*8 dist_trans(nnl), distn_trans(nnc)
	real*8 vort_trans(nnl),vortn_trans(nnc) 
	real*8 disfwall(nnl),disfwalln(nnc)

        pointer (ienptr,ien)
!!        pointer (selmptr,selm)
        pointer (rngptr,rng),(idnptr,idn)
        pointer (angptr,ang)
        pointer (slipangptr,slipang)
        pointer (mvsrfptr,mvsrf)
        pointer (xnptr,xn)
        pointer (xn1ptr,xn1),(x1ptr,x1)
        pointer (dn1ptr,dn1),(d1ptr,d1)
        pointer (d1_globptr,d1_glob)
	pointer (dn1_globptr,dn1_glob)
	pointer (dn1_transptr,dn1_trans)
        pointer (do_globptr,do_glob)
        pointer (bg1ptr,bg1),(p1ptr,p1)
        pointer (xn2ptr,xn2),(x2ptr,x2)
        pointer (dn2ptr,dn2),(d2ptr,d2)
        pointer (bg2ptr,bg2),(p2ptr,p2)
        pointer (fnptr,fn),(doptr,do),(dnoptr,dno)
        pointer (dno_globptr,dno_glob)
        pointer (cblptr,cbl)
        pointer (eddy_visptr,eddy_vis),(eddy_vis1ptr,eddy_vis1)
        pointer (dn1_sumptr,dn1_sum),(dn1_proptr,dn1_pro) 
        pointer (d1_sumptr,d1_sum),(d1_proptr,d1_pro)

        pointer (vloc1ptr,vloc1), (wg1ptr,wg1)
        pointer (vloc2ptr,vloc2), (wg2ptr,wg2)

        pointer  (id_tptr,id_t)
        pointer  (dn1_tptr,dn1_t) ,(d1_tptr,d1_t)
        pointer  (dno_tptr,dno_t) ,(do_tptr,do_t)
        pointer  (bg_tptr,bg_t)   ,(p_tptr,p_t)
        pointer  (fn_tptr,fn_t)  
        pointer  (vloc_tptr,vloc_t),(wg_tptr,wg_t)

	pointer  (ynode_ptr,ynode)
	pointer  (ynoden_ptr,ynoden)

	pointer  (umag_tptr,umag_trans), (umagn_tptr,umagn_trans)
	pointer  (he_tptr,he_trans), (hen_tptr,hen_trans) 
	pointer  (dist_tptr,dist_trans), (distn_tptr,distn_trans)
	pointer  (vort_tptr,vort_trans), (vortn_tptr,vortn_trans)
	pointer  (disfwall_ptr,disfwall), (disfwalln_ptr,disfwalln)
  

        real*8 mx(nsd,nnl),md(nsd,nnl),mid(nsd,nnc),mdn(nsd,nnc)
        pointer (mxptr,mx),(mdptr,md),(midptr,mid),(mdnptr,mdn)


!...TEMP ARRAYS FOR CHECKPOINTING
        real*8  dn(ndf,nnc),dn_t(ndf_t,nnc)
        pointer (dnptr,dn) ,(dn_tptr,dn_t)

!...ARRAYS FOR SURFACE MESH
        real*8 x_s(nsd,nn_s),d_s(nsd,nn_s),v_s(nsd,nn_s),p_s(nn_s)
        real*8 ds_old(nsd,nn_s)
        integer ien_s(nen_s,ne_s)
        pointer (x_sptr,x_s),(d_sptr,d_s),(v_sptr,v_s),(p_sptr,p_s)
        pointer (ds_optr,ds_old)
        pointer (ien_sptr,ien_s)
        integer mflag

        data time /30*0.0/

        integer i,j,idf,jdf,inl,nsize
        real*8 norm,res1,res2
        real*8 dminmax(2,5),dxminmax(2,5)
        real*8 norm_t,res1_t,res2_t,del_t,dummy
        integer temp(nrng), counter, irng, iface
!	integer eqn(ndf,nnc)
!       To generate probe directory

!       write(6,*)'I M IN HYPO,my_pe',my_pe,'nec',nec

        if (debug) call error("hypo: initialization",-999,.false.)
        call tstamp(time, 0)

        i = 1
        res1  = 0.0
        res2  = 0.0
        res1_t= 0.0
        res2_t= 0.0
        del_t = 0.0

!...READ IN MESH FILES
        ienptr     = malloc(nec*nen*isize)
!!        selmptr    = malloc(nee*nee*nec*fsize)
        rngptr     = malloc(nec*neface *isize)
        angptr     = malloc(nnc*fsize)
        slipangptr = malloc(nnc*fsize)
        xn1ptr     = malloc(nnc*nsd*fsize)
        xnptr      = malloc(nnc*nsd*fsize)
!.......skip for the time being
        call frif(mienfile,ien,nen,ne)
        call frif(mrngfile,rng,neface,ne)
        call frdf(mxyzfile,xn,nsd,nn)
        call frdf(mangfile,ang,1,nn)
        call frdf(msangfile,slipang,1,nn)
        call tstamp(time, 1)
        
!!        do irng=1,nrng
!!            temp(irng)=0
!!        enddo
!!
!!        do ielm=1,nec
!!            do iface=1,neface
!!                irng=rng(iface,ielm)
!!                if(irng.ne.0) then
!!                    temp(irng) = temp(irng)+1
!!                endif
!!            enddo
!!        enddo
!!        do irng=1,nrng
!!            write(10,*) irng,temp(irng)
!!        enddo
!!        return
!!        do ielm=1,nec
!!            write(10+my_pe,*) (ien(ie,ielm),ie=1,nen)
!!            write(20+my_pe,*) (rng(ie,ielm),ie=1,neface)
!!        enddo
!!        do inode=1,nnc
!!            write(30+my_pe,*) (xn(ii,inode),ii=1,nsd)
!!        enddo
!!        do inode=1,nnc
!!            write(40+my_pe,*) ang(inode)
!!        enddo

!!        do ielm=1,nec
!!            write(90+my_pe,*) (ien(ie,ielm),ie=1,nen)
!!            write(100+my_pe,*) (rng(ie,ielm),ie=1,neface)
!!        enddo

!!        do ie=1,nec
!!            do iface=1,neface
!!             write(110+my_pe,*) (facemap(iface,inface,elem_type),
!!     $		inface=1,nnface)
!!            enddo
!!        enddo
!!	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!	stop

!       if(my_pe.eq.4) write(6,*)'I M IN HYPO.F', xn(1,2064),xn(2,2064),xn(3,2064)
!       write(my_pe+81,122)((ien(i,j),i=1,6),j=1,nec)
! 122   format(6i8)                                     


!...CHECKPOINT ARRAYS
        dnptr = malloc(nnc*ndf*fsize)
        dn_tptr = malloc(nnc*ndf_t*fsize)

!...SETUP COMMUNICATION
        if (debug) call error("hypo: communication setup",-999,.false.)
        call treesetup(numprocs,my_pe)
!...Validating the data.
!!        do ielm=1,nec
!!            write(10+my_pe,*) (ien(ie,ielm),ie=1,nen)
!!        enddo
!!	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!	stop

        call commsetup(ien)
!...Validating the data.
!!        do ielm=1,nec
!!            write(20+my_pe,*) (ien(ie,ielm),ie=1,nen)
!!        enddo
!!	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!	stop

        cblptr = malloc(nnl*fsize)

!       write(61+my_pe,123)((xn(i,j),i=1,3),j=1,nnc)
! 123   format(3e18.9)
!       write(6,*)'WRITTEN FORT.40 '!******************************


        call pn(cbl)                      !commented for domain partioning

        call tstamp(time, 3)


!...PUT d_s and v_S out here to take care of restart problem
!     d_sptr = malloc(nsd*nn_s*fsize)
!     call fclear(d_s,nsd*nn_s)
!     v_sptr = malloc(nsd*nn_s*fsize)
!     call fclear(v_s,nsd*nn_s)
!     ds_optr = malloc(nsd*nn_s*fsize)
!     call fclear(ds_old,nsd*nn_s)

!!        if(meshmv.or.porosity) then
!!!...SETUP VOLUME SURFACE INTERFACE
!!            x_sptr = malloc(nsd*nn_s*fsize)
!!            p_sptr = malloc(nn_s*fsize)
!!            ien_sptr = malloc(ne_s*nen_s*isize)
!!            call fclear(p_s,nn_s)
!!            if(my_pe.eq.0) then
!!                call sfri(sienfile,ien_s,nen_s,ne_s)
!!                call sfrd(sxyzfile,x_s,nsd,nn_s)
!!            endif
 
!!           call f2s_setup(maxnnc,numprocs,my_pe,nn_s)

!!#if defined(UNI)
!!#else
!!            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!#endif

!!...SETUP MESH MOTION
!!            mnqm = nsd*nnc
!!            mdnptr = malloc(nsd*nnc*fsize)
!!            mvsrfptr = malloc(nnc*fsize)
!!            mdptr = malloc(nsd*nnl*fsize)
!!            mxptr = malloc(nsd*nnl*fsize)
!!            midptr= malloc(nsd*nnc*fsize)
!!            call gather(mx,xn,nsd)
!!            call fclear(mvsrf,nnc)
!!        if (debug) call error("hypo: meshmoving setup",-999,.false.)
!!        else
            midptr= malloc(1*fsize)  ! for call to genid
!!        endif

        call tstamp(time, 10)

!...SETUP BOUNDARY CONDITIONS
        if (debug) 
     $    call error("hypo: boundary condition setup",-999,.false.)
        fnptr = malloc(nnc*ndf*fsize)
        dnoptr = malloc(nnc*ndf*fsize)
        dno_globptr = malloc(nnc*ndf*fsize)
	dn1_globptr = malloc(nnc*ndf*fsize)
	dn1_transptr   = malloc(nnc*ndf*fsize)
!...SETUP INITIAL CONDITIONS
        if (debug) 
     $    call error("hypo: initial condition setup",-999,.false.)
        dn1ptr = malloc(nnc*ndf*fsize)
        d1ptr = malloc(ndf*nnl*fsize)
        doptr= malloc(ndf*nnl*fsize)
        x1ptr = malloc(nsd*nnl*fsize)
        d1_globptr = malloc(nnl*ndf*fsize)
        do_globptr = malloc(nnl*ndf*fsize)

!....INITIALIZE ARRAYS
        dn1_sumptr = malloc(nnc*ndf*fsize)                            !***************
        dn1_proptr = malloc(nnc*6*fsize)                              !***************
        d1_sumptr = malloc(nnl*ndf*fsize)                             !***************
        d1_proptr = malloc(nnl*6*fsize)                               !***************
        eddy_visptr = malloc(nnc*fsize)                               !***************
        eddy_vis1ptr=malloc(nnl*fsize)                                !***************

        call fclear(dn1_sum,ndf*nnc)                                    !***************
        call fclear(dn1_pro,6*nnc)                                      !***************
        call fclear(d1_sum,ndf*nnl)                                     !***************
        call fclear(d1_pro,6*nnl)                                       !***************
        call fclear(eddy_vis,1*nnc)                                     !***************
        call fclear(eddy_vis1,1*nnl)                                    !***************

        if (sptime) then
            dn2ptr = malloc(nnc*ndf*fsize)
            d2ptr = malloc(ndf*nnl*fsize)
            x2ptr = malloc(nsd*nnl*fsize)
        end if

        if(restart) then

            call frdf(dinfile1,dn1,ndf,nn)
!	    call equal(xn,xn1,nsd*nnc)
!!            do i=1,nnc
!!                write(30+my_pe,*)(dn1(j,i),j=1,ndf)
!!           end do
!!	   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!	   stop

            if (sptime) then
                call frdf(xinfile1,xn1,nsd,nn)
                call frdf(xinfile2,xn2,nsd,nn)
                call frdf(dinfile2,dn2,ndf,nn)
	   else
		call equal(xn,xn1,nsd*nnc)
            end if
        else
            call genic(xn,dn1)
!!	    do i=1,nnc
!!		write(30+my_pe,*)(dn1(j,i),j=1,ndf)
!!	   end do
!!	   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!	  stop
            call equal(xn,xn1,nsd*nnc)

!!            if (sptime) then
!!                call equal(xn,xn2,nsd*nnc)
!!                call equal(dn1,dn2,nnc*ndf)
!!            end if
        endif

!!!        do inode=1,nnc
!!!            write(40+my_pe,*) (dn1(ii,inode),ii=1,ndf)
!!!        enddo

 
        call tstamp(time,2)

!...LOCALIZE DATA
        if (debug) call error("hypo: localizations",-999,.false.)
        call gather (x1,xn1,nsd)
        call gather (d1,dn1,ndf)
        call gather (d1_sum,dn1_sum,ndf)
        call gather (d1_pro,dn1_pro,6)
        call gather (eddy_vis1,eddy_vis,1)

!!	do i=1,nnl
!!		write(50+my_pe,*)(d1(j,i),j=1,ndf)
!!	end do
!!	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!	stop


        if (sptime) then
            call gather (x2,xn2,nsd)	
            call gather (d2,dn2,ndf)
        end if
        call tstamp(time,4)

!...SETUP FOR TURBULENCE MODEL
        if(turb_model.eq.1.or.turb_model.eq.3) then
            hwallptr = malloc(fsize*nnc)
!	    ynodeptr = malloc(fsize*nnc)

            call fclear(hwall,nnc)
            turbptr = malloc(fsize*nnl)
            call fclear(turb,nnl)
!                call form_wall_func(ien,x1,d1,d1,rng,xn1,turb,hwall,
!     $                      viscosity,force_surface,1,nsd,ndf,nnl,nnc,
!     $                                nn,nen,nec,neface,my_pe,numprocs,
!     $                                0.5,0.5,maxnnc)
            call fclear(turb,nnl)
        endif  !setup for turbulence modelling ends

!	write(6,*)"outside form_wall_func"	
        if (debug) call error("hypo: turbulence setup",-999,.false.)

!       print *, ' im before genid'

!...SETUP idn ARRAY
        if (debug) call error("hypo: idn array setup",-999,.false.)
        idnptr = malloc(nnc*ndf*fsize)
        if (sptime) then 
            call genid(ien,idn,rng,mid)
        else
            call genid_sd(ien,idn,rng,mid)
        end if

!!	   do i=1,nnc
!!		write(60+my_pe,*) (idn(j,i),j=1,ndf)
!!	   end do	
!!	   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!	   stop


!!!        do ie=1,nec
!!!            do inen=1,nen
!!!                write(140+my_pe,*) (lm(idf,inen,ie),idf=1,ndf)
!!!            enddo
!!!        enddo
!!!        return
        counter = 0
!!       do inode=1,nnc
!!!            do ii=1,ndf
!!!                if(idn(ii,inode).ne.0) then
!!!                    write(50+my_pe,*) '1'
!!!                else
!!!                    write(50+my_pe,*) '0'
!!!                endif
!!!            enddo
!!!        enddo
!!!            write(50+my_pe,*) (idn(ii,inode),ii=1,ndf)
!!            do ii=1,ndf
!!                if(idn(ii,inode).eq.0.0) then
!!                    eqn(ii,inode) = 0
!!                else
!!                    counter = counter + 1
!!                    eqn(ii,inode) = counter
!!                endif
!!            enddo
!!           write(70+my_pe,*) (eqn(ii,inode),ii=1,ndf)
!!        enddo
!!	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!	stop

	
        write(6,*) 'mype,nnc,nec,nnl,nn',my_pe,nnc,nec,nnl,nn

        if (my_pe.eq.0) write (lookunit,2000) nq
!        if (my_pe.eq.0.and.meshmv) write (lookunit,2001) mnq

       if(my_pe.eq.0) write(6,*)'The num of equations are' ,nq


!...CONSTRUCT LOCAL SHAPE FUNCTIONS
        if (debug) call error("hypo:local shape functions",-999,.false.)
!        call mkshlocal
        call tstamp(time,1)

!...CHECK MESH DISTORTION
        heptr     = malloc(fsize*nec)
        volumeptr = malloc(fsize*nec)

!FIXME..uncomment following
!1	if(flag_he.gt.5) call ckvol(ien,x1)
!1	call ckvol(ien,x1)
!1        if (sptime) call ckvol(ien,x2)

        vloc1ptr   = malloc(ndf*nnl*fsize)
        bg1ptr     = malloc(ndf*nnc*fsize)
        p1ptr      = malloc(ndf*nnl*fsize)
        wg1ptr     = malloc(nqc*fsize)

        if (sptime) then
            bg2ptr = malloc(ndf*nnc*fsize)
            p2ptr  = malloc(ndf*nnl*fsize)
            vloc2ptr  = malloc(fsize*ndf*nnl)
            wg2ptr    = malloc(fsize*nqc)
        end if

!...INSERT BOUNDARY CONDITIONS
!!	murali
        if(.not.restart) then
            if (debug) 
     $            call error("hypo: boundary conditions",-999,.false.)
            call genf(ien,fn,rng,ang)
            call setd(dn1,fn,idn)
            if (sptime) call setd(dn2,fn,idn)
        endif
! 	murali


!!        do inode=1,nnc
!!            write(80+my_pe,*) (fn(ii,inode),ii=1,ndf)
!!        enddo
!!	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!	stop

        call gather (d1,dn1,ndf)
        call gather (d1_sum,dn1_sum,ndf)                               !***************
        call gather (d1_pro,dn1_pro,6)                                 !***************
        call gather (eddy_vis1,eddy_vis,1)                             !***************         
        if (sptime) call gather (d2,dn2,ndf)


!......SETUP FOR SPALART-ALMARAS TURBULENCE MODEL
!	write(6,*) "before turbulence model"
        if(turb_model.eq.3)then
!	write(6,*)"pointers in turbulence"
!......SETUP BOUNDARY CONDITIONS          
            fn_tptr = malloc(ndf_t*nnc*fsize)
            dno_tptr= malloc(ndf_t*nnc*fsize)
!......SETUP INITIAL CONDITIONS
            dn1_tptr= malloc(ndf_t*nnc*fsize)

!!!	can be used in block subroutine
            d1_tptr = malloc(ndf_t*nnl*fsize)
            do_tptr = malloc(ndf_t*nnl*fsize)
!......SETUP TRANSITION CONDITION
!            dist_tptr   = malloc(nnc*fsize)
!            inode_tptr  = malloc(nnc*isize)
!            ielem_tptr  = malloc(ntrans*isize)
!            vort_tptr   = malloc(nnl*fsize)
!           umag_tptr   = malloc(nnl*fsize)
!            he_tptr     = malloc(nnl*fsize)
!            nd_tptr     = malloc(ntrans*isize)
!           ielem_tptr  = malloc(ntrans*isize)
            ynoden_ptr   = malloc(fsize*nnl)
	    ynode_ptr   = malloc(fsize*nnc)

	    umag_tptr   = malloc(nnl*fsize)
	    umagn_tptr  =  malloc(nnc*fsize)

	    he_tptr     = malloc(nnl*fsize)
	    hen_tptr    = malloc(nnc*fsize)

	    dist_tptr   = malloc(nnl*fsize)
	    distn_tptr  = malloc(nnc*fsize)

	     vort_tptr   = malloc(nnl*fsize)
	     vortn_tptr  = malloc(nnc*fsize)

	    disfwall_ptr =  malloc(nnl*fsize)
	    disfwalln_ptr = malloc(nnc*fsize)

!	write(6,*)"HI"
!! fix it 
!        call fri(mwallnodefile,inode_trans,1,nn)
!	call fri(mtransnfile,nd_trans,1,ntrans)
!	call fri(mtransefile,ielem_trans,1,ntrans)
!        call fri(mrngfile,rng,neface,ne)
        call frdf(mylfile,ynode,1,nn)
	call frdf(mtransdistfile,distn_trans,1,nn)
!	call frd(mtransdistfile,dist_trans,1,nn)
!        call frd(mangfile,ang,1,nn)
!        call frd(msangfile,slipang,1,nn)
!	    mylfile,mwallnodefile
!! fix it 
!	call gather (inoden_trans,inode_trans,ndf_t)
!	murali
!	murali
!	call gather (distn_trans,dist_trans,ndf_t)
	call equal (ynoden,disfwalln,nnc)

         if(restart)then
            call frdf(dinfile_t,dn1_t,ndf_t,nn)
         else
            call genic_turb(dn1_t)
         end if

!
!          if(restart)then
!            write(6,*)"RESTART:ON,Reading RESTART file for TURB: fort.12"
!            open(unit=10,file='dataturb.in',status='old',iostat=ios)
!            if(ios.ne.0)then
!              print*,"RESTART FILE for TURBULENCE not found..."
!              print*,"ASSIGNING INITIAL CONDITIONS for TURBULENCE......."
!              print*," "
!	      close(10)
!	      call genic_turb(dn1_t)
!            else
!		call frd(dinfile_t,dn1_t,ndf_t,nn)
!              print*," "
!            endif
!          endif


!            if(restart)then
!                call frd(dinfile_t,dn1_t,ndf_t,nn)
!            else
!!		write(6,*)"before genic_turb"
!                call genic_turb(dn1_t)
!            end if

!	        write(6,*)"after genic_turb"
!......LOCALIZE DATA
            call gather (d1_t,dn1_t,ndf_t)
!......SETUP idn ARRAY
            id_tptr = malloc(ndf_t*nnc*fsize)
!	    write(6,*)"before genid_sd_turb"
            call genid_sd_turb(ien,id_t,rng)
	   if (my_pe.eq.0) write (6,*) nq_t
!	    write(6,*)"after genid_sd_turb"
            if (my_pe.eq.0) write (lookunit,1900) nq_t
!......INSERT BOUNDARY CONDITIONS
            if(.not.restart) then
!		call genf(ien,fn,rng,ang)
                call genf_turb(ien,fn_t,rng)
                call setd_turb(dn1_t,fn_t,id_t)
!		write(6,*)"assigned initial conditions"
            end if
            call gather (d1_t,dn1_t,ndf_t)
!	    write(6,*)"before trans_data"

!	call trans_data(ien,x1,d1,do,rng,
!     $         ielem_trans,umag_trans,vort_trans,
!     $         he_trans,ntrans)

!	    write(6,*)"after trans_data"
!            call gather(dist_trans,distn_trans,ndf_t)
!            call gather(ynode,ynoden,ndf_t)

        call gather (ynoden,ynode,1)
	call gather(dist_trans,distn_trans,1)
	call gather (disfwall, disfwalln,1)
        end if     !END OF LOOP FOR S-A TURBULENCE MODEL

!.......LM array isn't being generated.
!.......Just make sure why?

!.......Following must be enclosed in an if clause
        if(nsolution.le.nts) then
            call rotation(dn1,dn1_glob,slipang)
            call gather(d1_glob,dn1_glob,ndf)
            call transform(dn1_glob,dn1_trans,1)        !key = 1
        endif

!!	do i=1,nnc
!!	   write(120+my_pe,*) (dn1_trans(j,i),j=1,ndf)
!!	end do
!!	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!	stop

        its = 0

!        if (nprobe .gt. 0) then
!            do i=1,nprobe
!                flowrate(i)=0.0d0
!            enddo
!            call massflowrate(ien,d1)
!            call
!     $      MPI_REDUCE(flowrate_pe,flowrate,nprobe,MPI_DOUBLE_PRECISION,
!     $          MPI_SUM,0,MPI_COMM_WORLD,ierr)
!
!            call probe(x1,d1,flowrate)
!        endif

!!!        do inode=1,nnc
!!!            write(60+my_pe,*) (dn1(ii,inode),ii=1,ndf)
!!!            write(70+my_pe,*) (dn1_glob(ii,inode),ii=1,ndf)
!!!            write(80+my_pe,*) (dn1_trans(ii,inode),ii=1,ndf)
!!!        enddo

!.......isteady = 1 : unsteady calculation 
!.......For NEW Back Pressure, put 0.0d0 in 'time.in'
!.......For RESTART for same back pressure copy time.out to time.in
        if(restart) then
            open(unit=10,file='time.in',status='old')
            read(10,*) t_start
            close(10)
        else
            t_start = 0.0d0
        endif

!...TIMESTEP LOOP (BEGIN) ...........................................................!....NTS
                                                                        !    NTS
        if (debug) call error("hypo: timestep loop",-999,.false.)       !    NTS
                                                                        !    NTS
!!!        tt = t_start-dt

        do its = ntstart,nts                                            !    NTS
                                                                        !    NTS
!         nt = nt+dt                                                    !    NTS
!!!            tt = tt + dt                                                !    NTS
            tt = t_start + its*dt
!!!            if (mod(its,nsbref).eq.0) irefo = 1
            if (my_pe.eq.0) write(lookunit,2010) its                    !    NTS
                                                                        !    NTS
!.....................BACK PRESSURE IMPLEMENTATION......................
            if(ApplyBackPressure)then
                prs_exit_ini = BackPresInitial
                prs_exit_fin = BackPresFinal
                time_start   = 0.00
                time_exit    = 0.2
                p_down = prs_exit_ini+(prs_exit_fin-prs_exit_ini)*
     $                              (tt-time_start)/time_exit
                if(tt .gt. time_start+time_exit) p_down = prs_exit_fin

                write(6,6000)tt,p_down     ! Time , DownStream Pressure
!            else
!                write(6,6001)tt
            endif
	   if (my_pe.eq.0) write(6,6001)tt
!.......Writing time.out file
            if(iunsteady.eq.1.0)then      ! isteady = 1 for unsteady sol 
                open(unit=92,file='time.out',status='unknown')
                write(92,*)tt            ! Current time
                close(92)
            endif

!...COMPUTE ELEMENT LENGTH AND CHECK MESH DISTORTION                    !    NTS
!1            if (sptime) then
!1                if(its.eq.1.or.meshmv) call ckvol(ien,x2)               !    NTS
!1            endif                                                       !    NTS
                                   
!...UPDATE SOLUTION                                                                  !    NTS
            if (sptime) then
                call equal(d2,do,ndf*nnl)                               !    NTS
            else
                call equal(d1,do,ndf*nnl)
            end if
            call scatter(do,dno,ndf,2)

            if(turb_model.eq.3)then
                call equal(d1_t,do_t,ndf_t*nnl)
                call scatter(do_t,dno_t,ndf_t,2)
                bg_tptr = malloc(ndf_t*nnc*fsize)
                p_tptr = malloc(ndf_t*nnl*fsize)
                vloc_tptr = malloc(fsize*ndf_t*nnl)
                wg_tptr = malloc(fsize*nqc_t)
            end if
        
!...ROTATION BOUNDARY CONDITION                              !===============
!        if(irotate.eq.1)then
!               tt1=tt
!               nt1=nt*its
!           call rotatebc(ien,rng,d1,x1,nt1)
!           call scatter(d1,dn1,ndf,2)
!        end if

!...UPDATE BOUNDARY CONDITIONS FOR MESHMOTION                                        !    NTS
!......ITERATION LOOP (BEGIN) ..........................................!NIT NTS

            do iit = 1,nit                                              !NIT NTS
                if (debug) call error("hypo: nonlinear loop",-999,.false.)  !NIT NTS
                                                                        !NIT NTS

! added
!		call genf(ien,fn,rng,ang)
!		call setd(dn1,fn,idn)
! added
                call rotation(dn1,dn1_glob,slipang)
                call rotation(dno,dno_glob,slipang)
                call gather(d1_glob,dn1_glob,ndf)
                call gather(do_glob,dno_glob,ndf)
!...FORM RHS                                                                   !FLUID!NIT NTS
!                if(freezfluid) goto 1116
!       write(40+my_pe,602)((x1(i,j),i=1,3),j=1,nnl)
!       write(70+my_pe,602)((d1(i,j),i=1,2),j=1,nnl)
!       write(60+my_pe,602)((do(i,j),i=1,4),j=1,nnl)
! 602   format(2e16.9)
!1                if (sptime) then
!       write(70+my_pe,602)((x2(i,j),i=1,3),j=1,nnl)
!       write(80+my_pe,602)((d2(i,j),i=1,4),j=1,nnl)
!1                    call blkrhs(ien,x1,d1,x2,d2,do,p1,p2,vloc1,vloc2)   !FLUID!NIT NTS
!1                else
                    zlen = 0.2	!not used anywhere
                    call fclear(p1,ndf*nnl)
!                    call tempdummmy(p1)
                    call blk_iter_sd(ien,x1,d1_glob,do_glob,
     $                      p1,ang,rng,vloc1,d1_t,do_t)
!!                    call blkrhs_smag(ien,x1,d1,do,p1,vloc1,d1_t,do_t,
!!     $                  zlen,eddy_vis1)           !FLUID!NIT NTS
!1                end if


!...NEUMANN CONDITION                                                          !FLUID!NIT NTS
!             call neumann(x1,x2,d1,d2,p1,p2,ien,rng)                          !FLUID!NIT NTS
                call tstamp(time,6)                                     !FLUID!NIT NTS
!!                call mkblk(ien,block1,selm)
                call fclear(bg1,ndf*nnc)
                call fclear(wg1,ndf*nnc)

                call scatter(p1,bg1,ndf,1)                              !FLUID!NIT NTS
                call scatter(vloc1,wg1,ndf,1)

!!                do inode=1,nnc
!!                    write(85+my_pe,*) (bg1(idf,inode),idf=1,ndf)
!!                    write(80+my_pe,*)
!!     $                  ((blockg1(idf,jdf,inode),idf=1,ndf),jdf=1,ndf)
!!                enddo
!!                return 
!!                call scatter(vloc1,wg1,ndf,1)                           !FLUID!NIT NTS

                if (sptime) then
                    call scatter(p2,bg2,ndf,1)                          !FLUID!NIT NTS
                    call scatter(vloc2,wg2,ndf,1)                       !FLUID!NIT NTS
                end if

                call tstamp(time,5)                                     !FLUID!NIT NTS
                call filter(bg1,idn,ndf,nnc)                             !FLUID!NIT NTS
                if (sptime) call filter(bg2,idn,ndf,nnc)                 !FLUID!NIT NTS
                call tstamp(time,6)                                     !FLUID!NIT NTS
                call getnorm(bg1,bg1,res1,nqc)                          !FLUID!NIT NTS
                if (sptime) call getnorm(bg2,bg2,res2,nqc)              !FLUID!NIT NTS
                residue = sqrt((res1+res2)/nq)                              !FLUID!NIT NTS
!                if(my_pe.eq.0) write(6,*) 'residue',residue

                call tstamp(time, 7)                                    !FLUID!NIT NTS
!
!               if(iit.eq.1.and.its.eq.nts)then
!               call fwd("res1.out\0",bg2,ndf,nn) ! NODAL RESIDUALS
!               elseif(iit.eq.nit.and.its.eq.nts) then
!               call fwd("res2.out\0",bg2,ndf,nn) ! NODAL RESIDUALS
!               endif
!
                if(nodalres.and.mod(its,nsrest).eq.0)call resout(bg2)  !FLUID!NIT NTS

                    call fgmres_smag(bg1,x1,d1,p1,do,dn1,ien,idn,rng,
     $              wg1,d1_t, do_t,zlen,ang)     !FLUID!NIT NTS

!                write(70+my_pe,602)((bg1(i,j),i=1,2),j=1,nnl)
!       write(60+my_pe,602)((do(i,j),i=1,4),j=1,nnl)
! 602           format(2e16.9)


                call tstamp(time,8)                                     !FLUID!NIT NTS
!...UPDATE SOLUTION                                                            !FLUID!NIT NTS
                call updated(dn1,bg1,ndf,nnc,nnl)                       !FLUID!NIT NTS

!...Check for negative rho and T
                call limiter(dn1)

!               call updated_sum(dn1_sum,dn1,ndf,nnc)                          !FLUID!NIT NTS   !*************  
!               call updated_pro(dn1_pro,dn1,ndf,nnc)                          !FLUID!NIT NTS   !*************
                if (sptime) call updated(dn2,bg2,ndf,nnc)               !FLUID!NIT NTS
                call tstamp(time, 3)                                    !FLUID!NIT NTS
                call getnorm(bg1,bg1,res1,nqc)                          !FLUID!NIT NTS
                if (sptime) call getnorm(bg2,bg2,res2,nqc)              !FLUID!NIT NTS
                call tstamp(time, 7)                                    !FLUID!NIT NTS
                del = sqrt((res1+res2)/nq)                              !FLUID!NIT NTS
!		if(my_pe.eq.0) write(6,*) del
!              if(del.gt.1.0e+3)then                                      
!                       write(6,*) 'REDIDUAL IS BLOWING UP'
!                       stop
!               end if
                if(my_pe.eq.0) then                                     !FLUID!NIT NTS
                    write(lookunit,2020)its,iit,residue,del              !FLUID!NIT NTS
                    write(       6,2020)its,iit,residue,del              !FLUID!NIT NTS
                endif  
                if(del.gt.1.0e+5)then
                    write(6,*) 'RESIDUAL IS BLOWING UP'
                    stop
                end if                                                  !FLUID!NIT NTS
#if defined(UNI)
#else
                call MPI_BARRIER(MPI_COMM_WORLD,ierr)                   !FLUID!NIT NTS
#endif
                                                                        !FLUID!NIT NTS
!...GATHER UPDATED SOLUTION TO PROCESSOR LEVEL                          !FLUID!NIT NTS
                call gather (d1,dn1,ndf)                                !FLUID!NIT NTS
!             call gather (d1_sum,dn1_sum,ndf)                          !FLUID!NIT NTS    !*************
!             call gather (d1_pro,dn1_pro,6)                            !FLUID!NIT NTS    !*************
                if (sptime) call gather (d2,dn2,ndf)                    !FLUID!NIT NTS
                                                                        !FLUID!NIT NTS
!1116    if(my_pe.eq.0.and.freezfluid) print*,"Fluid is frozen!"
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF!FLUID!NIT NTS
                                                                        !NIT NTS
                if(del.lt.eps(1)) goto 29                               !NIT NTS
!250     continue
!		write(6,*)"turb_model",turb_model
                if(turb_model.eq.3)then

		if(trip) then

		call transition(x1,ien,d1,do)
#if defined(UNI)
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C...FIND TOTAL NUMBER OF EQUATIONS ON FIRST PROCESSOR AND BROADCAST TO REST
#else
	    call MPI_BCAST(vort_t,nslices,
     $		MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	    call MPI_BCAST(umag_t,nslices,
     $		MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	    call MPI_BCAST(he_t,nslices,
     $		MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif

!		write(6,*)my_pe,he_t(1),he_t(2),he_t(3),he_t(4),he_t(5),he_t(6)

 		call trans_data(ien,x1,umag_trans,vort_trans,
     $		   he_trans,ntrans)

		else 
		call fclear(umagn_trans,nnc)
		call fclear(vortn_trans,nnc)
		call fclear(hen_trans,nnc)

		call gather(umag_trans,umagn_trans,1)
		call gather(vort_trans,vortn_trans,1)
		call gather(he_trans,hen_trans,1)

		end if

!	do i=1,nnl
!		write(90+my_pe,*)disfwall(i),ynoden(i)
!	end do


                call blkrhs_turb_sd(ien,x1,d1,do,p_t,vloc_t,
     $                  d1_t,do_t,ynoden,
     $			umag_trans,he_trans,
     $                  dist_trans,vort_trans,
     $                  ntrans)
!		write(6,*)"after block"

		call fclear(bg_t,ndf_t*nnc)
                call fclear(wg_t,ndf_t*nnc)

                call scatter(p_t,bg_t,ndf_t,1)                              !FLUID!NIT NTS
                call scatter(vloc_t,wg_t,ndf_t,1)
!		call filter(bg_t,id_t,ndf_t,nnc)     
!               call scatter(p_t,bg_t,ndf_t,1)
!               call scatter(vloc_t,wg_t,ndf_t,1)
                call filter(bg_t,id_t,ndf_t,nnc)
                call getnorm (bg_t,bg_t,res1_t,nqc_t) !U HAVE TO CHECK THIS
!		write(6,*) "nq_t", nq_t
!		stop
                res1_t=dsqrt(res1_t/nq_t)
!		write(6,*)"number of turb equations",nq_t, res1_t
                call fgmres_turb_sd(bg_t,x1,d1,do,dn1,ien,id_t,rng,
     $                      d1_t,do_t,dn1_t,p_t,wg_t,ynoden,
     $                      umag_trans,he_trans,
     $                      vort_trans,dist_trans,
     $                      ntrans)

                call updated(dn1_t,bg_t,ndf_t,nnc)
		call limiter_turb(dn1_t)
                call getnorm(bg_t,bg_t,res2_t,nqc_t)!check nqc_t

                del_t=dsqrt(res2_t/nq_t)
!		write(6,*) "res2_t, nq_t del_t",res2_t, nq_t, del_t
                if(my_pe.eq.0)then
                        write(lookunit,2025)its,iit,res1_t,del_t        
                        write(       6,2025)its,iit,res1_t,del_t
                end if
#if defined(UNI)    
#else
                call MPI_BARRIER(MPI_COMM_WORLD,ierr)              
#endif

!.........GATHER UPDATED SOLUTION
                call gather(d1_t,dn1_t,ndf_t)
                end if                                                  !EDDY VISCOSITY DONE

!...ROTATION BOUNDARY CONDITION                              !===============
!        if(irotate.eq.1)then
!               nt1=nt*its
!          call rotatebc(ien,rng,d1,x1,nt1)
!           call scatter(d1,dn1,ndf,2)
!        end if

            end do
   29       continue                                                    !NIT NTS
!......ITERATION LOOP (END) .........................................................!NIT NTS

!...ROTATION BOUNDARY CONDITION                              !===============
!            if(irotate.eq.1)then
!                nt1=tt
!                call rotatebc(ien,rng,d1,x1,nt1)
!                call scatter(d1,dn1,ndf,2)
!            end if

            call updated_sum(dn1_sum,dn1,ndf,nnc)                       !************
            call updated_pro(dn1_pro,dn1,ndf,nnc)                       !************
!            call updated_vis(eddy_vis,eddy_vis1,1,nnc)                  !************
!
            call gather (d1_sum,dn1_sum,ndf)                            !************
            call gather (d1_pro,dn1_pro,6)                              !************
!            call scatter (eddy_vis1,eddy_vis,1,2)                       !************

                                                                        !    NTS
            call tstamp(time, 3)                                        !    NTS
!...COMPUTE AERO FORCES                                                 !    NTS
            if (debug) call error("hypo: drag",-999,.false.)            !    NTS
	           if (nprobe .gt. 0) call probe(xn1,dn1)
                call tstamp(time, 9)                                    !    NTS
                                                                        !    NTS
                                                                        !    NTS
            if(mod(its,nsrest).eq.0) then                              !    NTS
                if (debug) call error("OUTPUT FILES ",-999,.false.)     !    NTS
                idisk = idisk + 1                                       !    NTS
                call scatter(x1,xn1,nsd,2)                              !    NTS
                call equal(dn1,dn,ndf*nnc)                              !    NTS
!	data1.out file (restart file)
                call fwdf(doutfile1,dn,ndf,nn)                           !    NTS
!	data1.out file (restart file)
!       if (my_pe.eq.0) write(33,32)its,(dn(idf,22936),idf=1,4)
! 32            format(i3,2x,4e25.16)

!cccccc       call fwd(doutfilesum,dn1_sum,ndf,nn)                              !***************
!cccccc       call fwd(doutfilepro,dn1_pro,6,nn)                                !***************
!cccccc       call fwd(doutfilevis,eddy_vis,1,nn)                               !***************
                call equal(dn1,dn,ndf*nnc)                              !    NTS
!            call diskout_sd(dn)                                        !    NTS
                if(turb_model.eq.3)then                                 !    NTS
                    call equal(dn1_t,dn_t,ndf_t*nnc)                    !    NTS
                    call fwdf(doutfile_t,dn_t,ndf_t,nn)                  !    NTS
                    call equal(dn1_t,dn_t,ndf_t*nnc)                    !    NTS
!                 call diskout_turb_sd(dn_t)
                end if                                                  !    NTS
            endif

	    if(mod(its,nsrest).eq.0) then
                call fwdf(doutfilesum,dn1_sum,ndf,nn)
                call fwdf(doutfilepro,dn1_pro,6,nn)
            end if

!       if (my_pe.eq.4) write(33,32)its,(dn1(idf,2064),idf=1,4)
! 32            format(i3,2x,4e25.16)

            if(mod(its,nsolution).eq.0) then
                if (debug) call error("OUTPUT FILES ",-999,.false.)
                ndisk = ndisk+1                                         !    NTS 
                call ndiskout_sd(dn)
		if(turb_model.eq.3)then
		call diskout_turb_sd(dn_t)
		end if
            end if

        enddo                                                           !    NTS
!...TIMESTEP LOOP (END) .............................................................!....NTS

        call tstamp(time, 3)

!...WRITING FILES HERE TO TAKE CARE OF RESTART
!cccc      call fwd(doutfile1,dn1,ndf,nn)
!      call fwd(doutfilesum,dn1_sum,ndf,nn)                             !***************
!      call fwd(doutfilepro,dn1_pro,6,nn)                       !***************
!1        if (sptime) then
!1            call scatter(x1,xn1,nsd,2)
!1            call scatter(x2,xn2,nsd,2)
!1            call fwd(xoutfile1,xn1,nsd,nn)
!1            call fwd(xoutfile2,xn2,nsd,nn)
!1            call fwd(doutfile2,dn2,ndf,nn)
!1        end if

        if(turb_model.eq.3)then
            call fwdf(doutfile_t,dn1_t,ndf_t,nn)
        end if
        call tstamp(time, 2)

        call commdealloc

!...TIMINGS
        if (debug) call error("hypo: timings",-999,.false.)

        do i = 1,15
            time(16) = time(16)+time(i)
        enddo

        if (my_pe.eq.0) write(lookunit,2030) (time(i),i=1,16)
        call archinfo(lookunit)


        return

!       F O R M A T S
 2000 format(/"Fluid Mechanics Equations"
     &       /"-------------------------"
     &      //"Number of equations..................(nq) = ",i10)
 1900 format(/"Turbulent Viscosity Equations"
     &       /"-------------------------"
     &      //"Number of equations..................(nq_t) = ",i10)


 2001 format(/"Solid Mechanics Equations"
     &       /"-------------------------"
     &      //"Number of equations..................(nq) = ",i10)
 2010 format(/"Convergence information for time step ",i6
     &       /"-------------------------------------------"/)
 2020 format("its ",i6," iit ",i2," |residue| = ",e16.8," |del| = ",e16.8)
 2025 format("its ",i6," iit ",i2," |res_t| = ",e16.8," |del_t| = ",e16.8)
 2030 format(/"Solution time log"
     &       /"-----------------"
     &      //"I/O...................................... = ",f12.2,
     &       /"MISC..................................... = ",f12.2,
     &       /"COMMSETUP................................ = ",f12.2,
     &       /"GATHER................................... = ",f12.2,
     &       /"SCATTER.................................. = ",f12.2,
     &       /"BLOCK.................................... = ",f12.2,
     &       /"REDUCTION................................ = ",f12.2,
     &       /"SOLVER................................... = ",f12.2,
     &       /"DRAG..................................... = ",f12.2,
     &       /"CSR...................................... = ",f12.2,
     &       /"BLKMVP................................... = ",f12.2,
     &       /"CEBE(PARALLEL SOLVES).................... = ",f12.2,
     &       /"CEBE(COLORED  SOLVES).................... = ",f12.2,
     &       /"CEBE(COMM)............................... = ",f12.2,
     &       /"CEBE(FACTORIZE).......................... = ",f12.2,
     &      //"TOTAL.................................... = ",f12.2)
6000    format (/,'TIME : ',e9.4,2X,'P_down :',e11.4)
6001    format (/,1x,'TIME : ',e9.4,3x,"NO BACK_PRESSURE")

        end
