CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   
C       VAYU: INCOMPRESSIBLE NAVIER STOKES SOLVER       
C       Message Passing Version using MPI
C       Author: Vinay Kalro
C
C       HOOAH: PARACHUTE FLUID-STRUCTURE INTERACTION SOLVER
C                 -SPTIM3D: Main driver for FSI code with stabilized space-time
C                           formulation for the incompressible fluid (f90)
C                 -PARMEM:  Main subroutines for the cable-membrane parachute
C                           structures code
C                 -FEMLIB:  Library of f90 subroutines
C                 -PARLIB:  Library of c subroutines
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   
        program main
        implicit none
        include "global.h"
        real*8 rc
        integer irc
        
!...SETUP MPI MACHINE
        numprocs = 1
        my_pe= 0
!        call system ("rm -r PROBE ")
!        call system ("mkdir PROBE ")
#if defined(UNI)
#else
        call MPI_INIT(ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,my_pe,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
#endif

!       WRITE(6,*)'im in processor num  ', my_pe ,'  out of  ', numprocs
        if(my_pe.eq.0) then
            write(6,*)'NUMBER OF PNs:',numprocs
        endif


#if defined(UNI)
#else
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

!...INITIALIZE ALL PARAMETERS
        call initialize

!...READ INPUT
        call parseinput


C...WRITE DATA TO LOG
        call echoinput

C...RUN DRIVER
        call hypo

C...SHUTDOWN MPI MACHINE
#if defined(UNI)
#else
        call MPI_FINALIZE(ierr)
#endif

        stop
        end
