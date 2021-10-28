!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!       READ IN SIMULATION PARAMETERS (ALSO BROADCAST TO  ALL PN's)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine parseinput
        implicit none
        include "global.h"

        character*72 newname
        data newname /"fsi.in"/
        character*255 key, keyaux
        character*255 onoff
        integer i,ibc,idelta,idf,ip
        real*8 a1,a2,a3,a4,xnorm
        logical fctrl, isatty
        logical enough
        data enough /.false./


C...INPUT LOOP  

        do while (.not.enough)
C...EXIT PROGRAM  OR TERMINATE INPUT LOOP
            call getkey(key)
            if ((key.eq.'abort').or.(key.eq.'exit')) then
                !call exit(1)
                stop
            else if ((key.eq.'done').or.(key.eq.'quit')) then
                enough = .true.
C...READ  TIMESTEP AND INITIAL TIME
            else if (key.eq.'t_start') then
                call getreal(t_start)
            else if (key.eq.'dt') then 
                call getreal(dt)
C...REYNOLDS NUMBER
            else if (key.eq.'reynolds') then  
                call getreal(reynolds)
C...FLUID DENSITY
            else if (key.eq.'density') then  
                call getreal(rho)
C...READ FLAG FOR SETTING ZERO INITIAL CONDITION
            else if (key.eq.'ic_zero') then
                call getkey(onoff)
                ic_zero = (onoff.ne.'off')
C...QUADRATURE TYPE
            else if (key.eq.'iquad') then
                call getint(iquad)
C...MESH PARAMETERS
            else if (key.eq.'nn') then
                call getint(nn)
            else if (key.eq.'ne') then
                call getint(ne)
            else if (key.eq.'nrng') then
                call getint(nrng)
            else if (key.eq.'nen') then
                call getint(nen)
            else if ((key.eq.'ndf').or.(key.eq.'ndof')) then
                call getint(ndf)
            else if((key.eq.'ndf_t').or.(key.eq.'ndof_t')) then
                call getint(ndf_t)
            else if (key.eq.'nsd') then 
                call getint(nsd)
                maxnsd = nsd +1
            else if (key.eq.'iprecond') then
                call getint(iprecond)
            else if (key.eq.'isolver') then
                call getint(isolver)
            else if (key.eq.'ncolor') then
                call getint(ncolor)
            else if (key.eq.'nen_s') then
                call getint(nen_s)
            else if (key.eq.'ne_s') then
                call getint(ne_s)
            else if (key.eq.'nn_s') then
                call getint(nn_s)
            else if (key.eq.'nc') then
                call getint(nc)
            else if (key.eq.'nec') then
                call getint(nec)
            else if (key.eq.'nq') then 
                 call getint(nq)
C...LOOP SIZES
            else if (key.eq.'nit') then
                call getint(nit)
            else if (key.eq.'nts') then
                call getint(nts)
            else if (key.eq.'nslices') then
                call getint(nslices)
            else if (key.eq.'ntstart') then
                call getint(ntstart)
            else if (key.eq.'ntss') then
                call getint(ntss)
            else if (key.eq.'inner') then 
                call getint(inner)
            else if (key.eq.'outer') then 
                call getint(outer) 
            else if (key.eq.'minner') then 
                call getint(minner)
            else if (key.eq.'mouter') then 
                call getint(mouter) 
C...OUTPUT PARAMETERS
            else if (key.eq.'idisk') then
                call getint(idisk)
            else if (key.eq.'ndisk') then
                call getint(ndisk)
            else if (key.eq.'nsrest') then 
                call getint(nsrest)
            else if (key.eq.'nsolution') then
                call getint(nsolution)
            else if (key.eq.'nitbref') then 
                call getint(nitbref)
C...SURFACE FOR WHICH AERO COEFFS WILL BE COMPUTED
            else if (key.eq.'force_surface') then
                call getint(force_surface)
            else if (key.eq.'young') then
                call getreal(ym)
            else if (key.eq.'poisson') then
                call getreal(pr)
C...TIME INTEGRATION PARAMETER
            else if (key.eq.'alpha') then 
                call getreal(alpha)
C...Porosity model for membrane (v=por_coeff*dp)
            else if (key.eq.'porosity') then 
                call getreal(por_coeff); 
                porosity=.true.
C...STABILITY PARAMETERS
            else if (key.eq.'delta') then 
                call getint(idelta)
                call getreal(delta(idelta))
C...LINEAR SYSTEM SCALING OPTION
            else if (key.eq.'iscaling') then 
                call getint(iscaling)
C...LOGICAL FLAGS
            else if (key.eq.'sptime') then
                call getkey(onoff)
                sptime = (onoff.ne.'off')
            else if (key.eq.'freezstr') then
                freezstr = .true.
            else if (key.eq.'freezfluid') then
                freezfluid= .true.
            else if (key.eq.'freezmesh') then
                freezmesh = .true.
            else if (key.eq.'fsi_press_old') then
                fsi_press_current = .false.
            else if (key.eq.'fsi_press_current') then
                fsi_press_current = .true.
            else if (key.eq.'sprng') then
                call getint(sprng)
            else if (key.eq.'incdisp') then
                call getkey(onoff)
                incdisp = (onoff.ne.'off')
            else if (key.eq.'restart') then
                call getkey(onoff)
                restart = (onoff.ne.'off')
            else if (key.eq.'debug') then 
                call getkey(onoff)
                debug = (onoff.ne.'off')
            else if (key.eq.'steady') then !......steady state
                call getkey(onoff)
                steady = (onoff.ne.'off')
            else if (key.eq.'meshmv') then !.....fsi meshmotion
                call getkey(onoff)
                meshmv = (onoff.ne.'off')
            else if (key.eq.'ffmeshmv') then !.....freefall meshmotion
                call getkey(onoff)
                ffmeshmv = (onoff.ne.'off')
            else if (key.eq.'nodalres') then !.....dump nodal residuals
                call getkey(onoff)
                nodalres = (onoff.ne.'off')
            else if (key.eq.'mesh_inc') then !.....incremental meshmotion
                call getkey(onoff)
                mesh_inc = (onoff.ne.'off')
            else if (key.eq.'chi') then !.....mesh stiffening factor
                call getreal(chi)
            else if (key.eq.'stokes') then !.....stokes flow
                call getkey(onoff)
                stokes = (onoff.ne.'off')
            else if (key.eq.'part') then !......partitioned mesh
                call getkey(onoff)
                part = (onoff.ne.'off')
            else if (key.eq.'timing')then 
                call getkey(onoff)
                timing = (onoff.ne.'off')
            else if (key.eq.'title') then 
                call getkey(title)
            else if (key.eq.'eps') then
                call getint(ibc)
                call getreal(eps(ibc))
C...THERMODYNAMIC CONDITIONS
            else if (key.eq.'shock_cap') then       !...whether to include shock capturing?
                shock_cap = .true.
            else if (key.eq.'ShockCapFactor') then
                call getreal(ShockCapFactor)
            else if (key.eq.'AdiabaticWall') then
                AdiabaticWall = .true.
            else if (key.eq.'WallTemp') then
                call getreal(WallTemp)
            else if (key.eq.'ApplyBackPressure') then
                ApplyBackPressure = .true.
            else if (key.eq.'BackPresInitial') then
                call getreal(BackPresInitial)
            else if (key.eq.'BackPresFinal') then
                call getreal(BackPresFinal)
            else if (key.eq.'Sutherland_Law') then
                Sutherland_Law = .true.
            else if (key.eq.'ApplyBleed') then
                ApplyBleed = .true.
            else if (key.eq.'viscosity') then
                call getreal(viscosity)
            else if (key.eq.'mach_num') then
                call getreal(rmach)
            else if (key.eq.'flag_he') then
                call getint(flag_he)
C...MORE LOGICAL FLAGS
            else if (key.eq.'diag') then
                diag = .true.
            else if (key.eq.'block_diag') then
                block_diag = .true.
            else if (key.eq.'formatted') then
                formatted = .true.
            else if (key.eq.'explicit') then
                explicit = .true.
C...GMRES
            else if (key.eq.'gmres') then
                gmres = .true.
            else if (key.eq.'numcl') then
                call getint(numcl)
            else if (key.eq.'ngmr') then
                call getint(ngmr)
            else if (key.eq.'ikg') then
                call getint(ikg)
            else if (key.eq.'epsgm') then
                call getreal(epsgm)
            else if (key.eq.'epsrt') then
                call getreal(epsrt)
C...NEUMANN LOAD FLAG SETUP
            else if (key.eq.'neumann') then
                call getint(ibc)
                neum(ibc) = 1
C...TURBULENCE MODEL COEFF.
            else if (key.eq.'turbulent') then
                turbulent = .true.
            else if (key.eq.'turb_model') then
                call getkey(key)
                if(key.eq.'none') turb_model = 0
                if(key.eq.'prandtl') turb_model = 1
                if(key.eq.'smagor') turb_model = 5
                if(key.eq.'spalart')turb_model=3

C...BOUNDARY CONDITION TYPE SETUP(Modified)
!...1 = rho
!...2 = u
!...3 = v
!...4 = w
!...5 = e
            else if (key.eq.'pnode') then
                call getint(pnode)
            else if (key.eq.'allfixed') then
                call getint(ibc)
                bc(1,ibc) = 1
                bc(2,ibc) = 1
                bc(3,ibc) = 1
                bc(4,ibc) = 1
                bc(5,ibc) = 1
		if(turb_model.eq.3) bc_t(1,ibc)=1
            else if (key.eq.'noslip') then
                call getint(ibc)
                bc(2,ibc) = 1
                bc(3,ibc) = 1
                bc(4,ibc) = 1
		if(turb_model.eq.3) bc_t(1,ibc)=1
!		write(6,*) bc_t(1,ibc)
            else if (key.eq.'ufix') then
                call getint(ibc)
                bc(2,ibc) = 1
            else if (key.eq.'vfix') then
                call getint(ibc)
                bc(3,ibc) = 1
            else if (key.eq.'wfix') then
                call getint(ibc)
                bc(4,ibc) = 1
!!murali
!!		bc(2,ibc) = 0
!!		bc(3,ibc) = 1
	    else if (key.eq.'slip') then
		call getint(ibc)
		bc(2,ibc) = 0
                bc(3,ibc) = 1
                bc(4,ibc) = 1
!!murali
!!            else if (key.eq.'pfix') then
!!                call getint(ibc)
!!!find a proper mapping first
!!!                     bc(4,ibc) = 1
!!murali
            else if (key.eq.'ufree') then
                call getint(ibc)
                bc(3,ibc) = 1
                bc(4,ibc) = 1
            else if (key.eq.'vfree') then
                call getint(ibc)
                bc(2,ibc) = 1
                bc(4,ibc) = 1
            else if (key.eq.'wfree') then
                call getint(ibc)
                bc(2,ibc) = 1
                bc(3,ibc) = 1
!!!             else if (key.eq.'inflow') then
!!!                     call getint(ibc)
!!!                     bc(1,ibc) = 1
!!!                     bc(2,ibc) = 1
!!!                     bc(3,ibc) = 1
!!!
!!!                     if(turb_model.eq.3)then
!!!                     bc_t(1,ibc)=1
!!!                     end if
!!!             else if (key.eq.'outflow'.or.key.eq.'free') then
!!!                      call getint(ibc)
!!!                     bc(1,ibc) = 0
!!!                     bc(2,ibc) = 0
!!!                     bc(3,ibc) = 0
!!!                     bc(4,ibc) = 0
!!!
!!!             else if (key.eq.'top') then
!!!                     call getint(ibc)
!!!                     bc(1,ibc)=1
!!!                     bc(2,ibc)=1
!!!                     bc(3,ibc)=1
!!!                     bc(4,ibc)=0
!!!                     if(turb_model.eq.3)then
!!!c                        bc_t(1,ibc)=1
!!!                        end if
!!!             else if (key.eq.'mufix') then
!!!                     call getint(ibc)
!!!                     bcm(1,ibc) = 1
!!!             else if (key.eq.'mvfix') then
!!!                     call getint(ibc)
!!!                     bcm(2,ibc) = 1
!!!             else if (key.eq.'mwfix') then
!!!                     call getint(ibc)
!!!                     bcm(3,ibc) = 1
!!!             else if (key.eq.'mufree') then
!!!                     call getint(ibc)
!!!                     bcm(2,ibc) = 1
!!!                     bcm(3,ibc) = 1
!!!             else if (key.eq.'mvfree') then
!!!                     call getint(ibc)
!!!                     bcm(3,ibc) = 1
!!!                     bcm(1,ibc) = 1
!!!             else if (key.eq.'mwfree') then
!!!                     call getint(ibc)
!!!                     bcm(1,ibc) = 1
!!!                     bcm(2,ibc) = 1
!!!             else if (key.eq.'mfix') then
!!!                     call getint(ibc)
!!!                     bcm(1,ibc) = 1
!!!                     bcm(2,ibc) = 1
!!!                     bcm(3,ibc) = 1
            else if (key.eq.'freestream') then
                do idf = 1,ndf
                    call getreal(ic(idf,nrng+1))
                end do
                if(turb_model.eq.3)
     &                  call getreal(ic_t(1,nrng+1))


C...GET INFLATION PRESSURE
            else if (key.eq.'inflation') then
               call getreal(inf_pres)
C...GET DISPLACEMENTS
            else if (key.eq.'displac') then
                call getint(i)
                do idf=1,4
                   call getreal(ic(idf,i))
                end do

C...GRAVITY VECTOR      
            else if (key.eq.'gravity') then
                call getreal(gravity(1))
                call getreal(gravity(2))
                call getreal(gravity(3))
                xnorm = 
     $          sqrt(gravity(1)**2 + gravity(2)**2 + gravity(3)**2)
                if(xnorm.gt.epsilon) then
                        gravity(1) = gravity(1)/xnorm
                        gravity(2) = gravity(2)/xnorm
                        gravity(3) = gravity(3)/xnorm
                endif
C...BL THICKNESS FOR MESHMOTION
            else if (key.eq.'blyr') then
                call getreal(blyr)
!!C...PROBES SETTING FOR VELOCITY AND PRESSURE AT CERTAIN NODES
            else if (key.eq.'probe') then
                call getint(nprobe)
                do i=1,nprobe
                   call getint(node_probe(i))
                enddo
	   else if (key.eq.'trip') then
		   trip = .true.
	          call getint(ntelem)
!		 write(6,*)ntelem
		do i=1,ntelem
		     call getint(ielem_t(i))
		end do
C...ROTATION TO THE BODY
            else if (key.eq.'rotate') then
                call getint(irotate)
                call getreal(vel)
                call getreal(perc)
C...DOMAIN PARTITION
            else if (key.eq.'partition') then
                call getint(ipartition)

C...FILE NAMES
            else if (key.eq.'surf.ien') then
                call getkey(key)
                sienfile = key
            else if (key.eq.'surf.xyz') then
                call getkey(key)
                sxyzfile = key
            else if (key.eq.'PATH-mien') then
                call getkey(key)
                mienfile = key
            else if (key.eq.'PATH-mxyz') then
                call getkey(key)
                mxyzfile = key
            else if (key.eq.'PATH-mrng') then
                call getkey(key)
                mrngfile = key
            else if (key.eq.'PATH-mang') then
                call getkey(key)
                mangfile = key
            else if (key.eq.'PATH-msang') then
                call getkey(key)
                msangfile = key
!            else if (key.eq.'PATH-meshb') then
!                call getkey(key)
!                meshbfile = key
            else if (key.eq.'PATH-myl') then
                write(6,*)key
                call getkey(key)
                mylfile = key
            else if (key.eq.'PATH-mwallnode') then
                call getkey(key)
                mwallnodefile = key
!             else if (key.eq.'PATH-transnodes') then
!                call getkey(key)
!                mtransnfile = key
!            else if (key.eq.'PATH-transelem') then
!                call getkey(key) 
!                mtransefile = key
           else if (key.eq.'PATH-transdist') then
                call getkey(key)
                mtransdistfile = key 
            else if (key.eq.'color') then
                call getkey(key)
                colorfile = key
            else if (key.eq.'dinfile1') then
                call getkey(key)
                dinfile1 = key
            else if (key.eq.'dinfile2') then
                call getkey(key)
                dinfile2 = key
            else if(key.eq.'dinfile_t')then
                call getkey(key)
                dinfile_t = key
            else if (key.eq.'xinfile1') then
                call getkey(key)
                xinfile1 = key
            else if (key.eq.'xinfile2') then
                call getkey(key)
                xinfile2 = key
            else if (key.eq.'doutfile1') then
                call getkey(key)
                doutfile1 = key
            else if (key.eq.'doutfile2') then
                call getkey(key)
                doutfile2 = key
            else if(key.eq.'doutfile_t')then
                call getkey(key)
                doutfile_t= key
            else if(key.eq.'doutfilesum')then       !***********
                call getkey(key)
                doutfilesum= key
            else if(key.eq.'doutfilepro')then       !***********
                call getkey(key)
                doutfilepro= key
            else if(key.eq.'doutfilevis')then       !***********
                call getkey(key)
                doutfilevis= key
            else if (key.eq.'xoutfile1') then
                call getkey(key)
                xoutfile1 = key
            else if (key.eq.'xoutfile2') then
                call getkey(key)
                xoutfile2 = key
            else if (key.eq.'minfile') then
                call getkey(key)
                minfile = key
            else if (key.eq.'moutfile') then
                call getkey(key)
                moutfile = key
            else if ((.not.enough).and.(key.ne.' ').and.(key(1:1).ne.'#')) then
                call error("parseinput: strange keyword "//key, -999, .false.)
            end if

        end do

	write(6,*)"nen ndof",nen, ndf

C...
        reynolds=rho*reynolds

C...FURTHER DEFAULTS
        if (nsrest.eq.0) nsrest = nsrest + 1
        if (nsolution.eq.0) nsolution = nsolution + 1
        if (steady) then
           alpha = 1.0
           iunsteady=0.0
        else
           iunsteady=1.0
        endif
        if (stokes) then
           inostokes=0.0
        else
           inostokes=1.0
        endif
        if (nc.eq.0) nc = numprocs

C...ELEMENT SETUP       
        if (nen.eq.0) nen = 6

!...neface = number of elemental face
!...nnface = number of nodes per face
        if ((nsd.eq.2).and.(nen.eq.3)) then
            elem_type = elem_tria_sd
            if (sptime) elem_type = elem_tria_st
            neface = 3
            nnface = 2
        else if ((nsd.eq.2).and.(nen.eq.4)) then
            elem_type = elem_quad_sd
            if (sptime) elem_type = elem_quad_st
            neface = 4
            nnface = 2
        else if ((nsd.eq.3).and.(nen.eq.4)) then
            elem_type = elem_tetr_sd
            if (sptime) elem_type = elem_tetr_st
            neface = 4
            nnface = 3
        else if ((nsd.eq.3).and.(nen.eq.6)) then
            elem_type = elem_wedg_sd
            neface = 5
            nnface = 4
        
        else if ((nsd.eq.3).and.(nen.eq.8)) then
            elem_type = elem_hexa_sd
            if (sptime) elem_type = elem_hexa_st
            neface = 6
            nnface = 4
!	murali
	else if ((nsd.eq.3).and.(nen.eq.18)) then
            elem_type = elem_wedg_quad_sd
            neface = 5
            nnface = 9
        end if
!	murali

        if (ndf.eq.0) ndf = 5

C...NQ ESTIMATE 
        if (nq.eq.0) nq = ndf * nn


C...CLUSTER SIZES
        if (nec.eq.0) nec = (ne - 1) / nc + 1
        maxnec = nec
        if (my_pe.eq.numprocs-1) nec = ne - (numprocs - 1) * maxnec
        nnc = (nn - 1) / nc + 1
        maxnnc = nnc
        if (my_pe.eq.numprocs-1) nnc = nn - (numprocs - 1) * maxnnc

        nqc = nnc * ndf 
 
        if(turb_model.eq.3)nqc_t=nnc*ndf_t
c       write(6,*)'max turb equations ',nqc_t

C...DOMAIN PARTITION
        if (ipartition.eq.1) then
            open(unit=8,file='part_data')
            read(8,*)(ipart(ip),ip=1,numprocs)
c           read(8,*)(jpart(ip),ip=1,numprocs)
            nec = ipart(my_pe+1)
            maxnec = nec
c           nnc = jpart(my_pe+1)
c           maxnnc = nnc
c           nqc = nnc*ndf
c           if(my_pe.eq.9)write(6,*)'IN PARSEINPUT- MAXNEC',maxnec  
        end if
        close(8)
C...PROBES FOR VELOCITY AND PRESSURE OUTPUT AT SOME NODES, AND MASS
C...FLOW RATE INFORMATION

!        open(unit=10,file=meshbfile,status='unknown')
!        read(10,*) nprobe,nely,nelz
!        nmeshb = nprobe*nely*nelz
!        if(nprobe.gt.probepts) then
!            write(6,*) 'Increase the value of probepts'
!            stop
!        endif
!        if(nmeshb.gt.pmpts) then
!            write(6,*) 'Increase the value of pmpts'
!            stop
!        endif
!        do i=1,nprobe
!            read(10,*) node_probe(i)
!        enddo
!!.......cel is the containing element
!        do i=1,nmeshb
!            read(10,200) cel(i),dy(i),dz(i)
!            read(10,220) msh(1:6,i)
!        enddo
!        close(10)
        
200     format(i7,5x,2e16.9)
220     format(6e16.9)

        if (nprobe .gt. 0) then
           do i=1, nprobe
!              ! proc on which probe lies
              probe_pe(i)=(node_probe(i)-1)/maxnnc       
!              ! local node no for that probe
              node_pe(i) =node_probe(i)-probe_pe(i)*maxnnc     
             if (my_pe .eq. 0) write(6,*)i,node_probe(i),probe_pe(i),node_pe(i)
	 end do
	end if

	if (ntelem .gt. 0) then
	    do i=1, ntelem
	    ! proc on which element (transition) lies
	      itrne_pe(i)=(ielem_t(i)-1)/maxnec
!           ! local node no for that element
              ielem_pe(i) =ielem_t(i)-itrne_pe(i)*maxnec
!	     write(6,*)i,ielem_t(i),trne_pe(i),elem_pe(i)
           end do
        end if

!        if (nmeshb.gt.0) then
!            do i=1,nmeshb
!            !...proc on which containing element lies
!                cel_pe(i) = (cel(i)-1)/maxnec
!!FIXME......BASED ON THE ASSUMPTION THAT NODES ARE NUMBERED
!!CONSEQUTIVELY
!            !...local element number
!                el_pe(i)  = cel(i) - cel_pe(i)*maxnec
!            enddo
!        endif
                


C...INITIAL CONDITIONS SETUP

        dtinv = 1.0/dt
        if(steady) dtinv = 0.0
        oma  = 1.0 - alpha
        nee  = ndf*nen

!	write(6,*)"nee",nee

        if(debug.and.my_pe.eq.0) then
            write(6,*)'INTEGERS',nn,ne,nen,nsd,ndf
            write(6,*)'REALS',t_start,dt
            write(6,*)'LOGICALS',stokes,steady,debug
            write(6,*) (ic(idf,nrng+1),idf=1,ndf)
        endif

        return
        end
