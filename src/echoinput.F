c       cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine echoinput

        implicit none
        include "global.h"

        character*3 yon
        integer io, idelta, irng, isd, idf, i

        if (my_pe.ne.0) return

        io = lookunit

        write(io,'("include all the new parameters in here")')
        write(io,'("3D compressible semi-discrete code")')
        write(io,*) title
!	write(io,'("3D incompressible space-time code")')
        write(io,'(/"Control information")')
        write(io,'( "-------------------"/)')
        write(io,'(" Nodes......................(nn) = ",i9)') nn
        write(io,'(" Elements...................(ne) = ",i9)') ne
        write(io,'(" Time steps....................(nts) = ",i5)') nts
        write(io,'(" Starting time step........(ntstart) = ",i5)') ntstart
        write(io,'(" Solution saved after steps (nsolution) = ",i5)') nsolution
        write(io,'(" Iterations....................(nit) = ",i5)') nit
        write(io,'(" Time steps b/output.......(nsrest) = ",i5)') nsrest
        write(io,'(" Space dimensions..............(nsd) = ",i5)') nsd
        write(io,'(" Degrees of freedom............(ndf) = ",i5)') ndf
        write(io,'(" Number of element nodes.......(nen) = ",i5)') nen
        write(io,'(" Number of element face.....(neface) = ",i5)') neface
        write(io,'(" Integration points..........(iquad) = ",i5)') iquad
        write(io,'(" Scaling..................(iscaling) = ",i5)') iscaling
        write(io,'(" Inner GMRES iterations......(inner) = ",i5)') inner
        write(io,'(" Outer GMRES iterations......(outer) = ",i5)') outer
        write(io,'(" Debug mode..................(debug) = ",a5)') yon(debug)
        write(io,'(" Performance ouput..........(timing) = ",a5)') yon(timing)
        write(io,'(" Restart...................(restart) = ",a5)') yon(restart)
        write(io,'(" FSI Meshmotion.............(meshmv) = ",a5)') yon(meshmv)
        write(io,'(" Freefall Meshmotion......(ffmeshmv) = ",a5)') 
     .          yon(ffmeshmv)
        write(io,'(" Nodal residuals..........(nodalres) = ",a5)') yon(nodalres)
        write(io,'(" Sptime.....................(sptime) = ",a5)') yon(sptime)
        write(io,'(" Steady state...............(steady) = ",a5)') yon(steady)
        write(io,'(" Stokes problem.............(stokes) = ",a5)') yon(stokes)
        write(io,'(" Initial condition.........(ic_zero) = ",a5)') yon(ic_zero)
        write(io,'(" Block diagonal.........(block_diag) = ",a5)') yon(block_diag)
        write(io,'(" Diagonal.....................(diag) = ",a5)') yon(diag)
        write(io,'(" GMRES.......................(gmres) = ",a5)') yon(gmres)
        write(io,'(" Time step................ = ",e15.8)') dt
        write(io,'(" Initial time............. = ",e15.8)') t_start
        write(io,'(" Reynolds number.......... = ",e15.8)') reynolds
        write(io,'(" Wall temprature.......... = ",e15.8)') WallTemp
        write(io,'(" Shock capturing.........(shock_cap) = ",a5)') yon(shock_cap)
        if(shock_cap) then
            write(io,'(" Shock Capturig factor.... = ",e15.8)') ShockCapFactor
        endif
        write(io,'(" ApplyBackPressure.................. = ",a5)') yon(ApplyBackPressure)
        if(ApplyBackPressure) then
            write(io,'(" Back pressure initial.... = ",e15.8)') BackPresInitial
            write(io,'(" Back pressure final...... = ",e15.8)') BackPresFinal
        endif
        write(io,'(" AdiabaticWall...................... = ",a5)') yon(AdiabaticWall)
        write(io,'(" ApplyBleed......................... = ",a5)') yon(ApplyBleed)
        write(io,'(" Sutherland_Law..................... = ",a5)') yon(Sutherland_Law)
        write(io,'(" viscosity................ = ",e15.8)') viscosity
        write(io,'(" rmach.................... = ",e15.8)') rmach
        write(io,'(" epsgm.................... = ",e15.8)') epsgm
        write(io,'(" epsrt.................... = ",e15.8)') epsrt
        write(io,'(" formatted.......................... = ",a5)') yon(formatted)
        write(io,'(" explicit........................... = ",a5)') yon(explicit)
        write(io,'(" turbulent.......................... = ",a5)') yon(turbulent)
        write(io,'(" flag_he............................ = ",i5)') flag_he
        write(io,'(" numcl.............................. = ",i5)') numcl
        write(io,'(" ngmr............................... = ",i5)') ngmr
        write(io,'(" ikg................................ = ",i5)') ikg
        if (.not. sptime) 
     .  write(io,'(" Crank-nick param(alpha) . = ",e15.8)') alpha
        do idelta=0,ndelpad
            if (delta(idelta).ne.0.0) then
                write(io,'(" Delta_",i1,".................. = ",e15.8)')
     &                  idelta, delta(idelta)
            end if
        end do

        if (turb_model.eq.1)  then
            write(io,'(/"Baldwin Lomax Turbulence Model")')
            write(io,'( "------------------------------------"/)')
            write(io,'( " Kappa........(turb_kappa) = ",e15.8)') 0.41
        else if (turb_model.eq.5)  then
            write(io,'(/"Hitachi Smagorinsky Turbulence Model")')
            write(io,'( "------------------------------------"/)')
            write(io,'( " Kappa........(turb_kappa) = ",e15.8)') 0.15
        else if (turb_model.eq.3)then
           write(io,'(/"Spalart-Allamras  Turbulence Model")')
           write(io,'( "------------------------------------"/)')
        else    
            write(io,'(/"No turbulence modeling")')
        end if

        do isd=1,nsd
            if (gravity(isd).ne.0.0) then
                write(io,'(" Gravity_",i1,".......(gravity) = ",e15.8)') 
     .          isd, gravity(isd)
            end if
        end do

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        write(io,'(/"Boundary condition information" )')
        write(io,'( "------------------------------"/)')
        if(nsd.eq.3) then
            write(io,12)
 12   format('                     rho     u      v      w      e')
            write(io,13)
 13   format('                     ---    ---    ---    ---    ---')
            do i=1,nrng
                write(io,14) i,bc(1,i),bc(2,i),bc(3,i),bc(4,i),bc(5,i)
 14   format (' boundary no. ',i2,3x,5(3x,i1,3x))
            enddo
            do i=1,nrng
                write(io,15) i,ic(1,i),ic(2,i),ic(3,i),ic(4,i),ic(5,i)
 15   format (' boundary no. ',i2,3x,5(f6.3,1x))
            enddo
        else
            write(io,16)
 16   format('                      u      v      p ')
            write(io,17)
 17   format('                     ---    ---    ---')
            do i=1,nrng
                write(io,18) i,bc(1,i),bc(2,i),bc(3,i)
 18   format (' boundary no. ',i2,3x,3(3x,i1,3x))
            enddo
            do i=1,nrng
                write(io,19) i,ic(1,i),ic(2,i),ic(3,i)
 19   format (' boundary no. ',i2,3x,3(f6.3,1x))
            enddo
        endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        return
        end

