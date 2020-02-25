!...IDENTIFY NODES ON CLUSTER BOUNDARIES 
        subroutine pn(cbl) 
        implicit none
        include "global.h"
        include "fmalloc.h"

        real*8  cbl(nnl)
        real*8  cbg(nnc)
        pointer (cbgptr,cbg)

        integer in,one

        one = 1
!                       write(41+my_pe,*)nnc,nnl 

        cbgptr  = malloc(nnc*fsize)

        do in = 1,nnl
                cbl(in) = 1.0   
        enddo

        call scatter(cbl,cbg,one,1)

        call tstamp(time, 5)
!                if(my_pe.eq.0)write(6,*)'IIIIIIIIIIIIIIIIII am at 1 '!*****************************
        call gather(cbl,cbg,one)
!                if(my_pe.eq.0)write(6,*)'IIIIIIIIIIIIIIIIII am at 2 '!**************************

        call tstamp(time, 4)

!       if(my_pe.eq.0) then
!       write(6,*) 'PNTEST'
!       write(6,100) (cbg(in),in=1,nnc)
! 100   format(8f8.4)
!       endif



        call free(cbgptr)


        return
        end


