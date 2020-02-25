        subroutine filter(bg,id,ndf,nn)
        implicit none
        integer inl,idf,ndf,nn
        real*8 bg(ndf,nn),id(ndf,nn)
        
        do inl = 1,nn
                do idf = 1,ndf
                         bg(idf,inl) = id(idf,inl)*bg(idf,inl)
                enddo
        enddo

        endsubroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subroutine limiter(dn)
        implicit none
        include "global.h"
        real*8 dn(ndf,nnc)
        real*8 ro,uvel,vvel,wvel,e,temp,fac
        integer i

        fac = 1.0d0/(gamma*(gamma-1.0d0))

        do i=1,nnc
            ro    = dn(1,i)
            uvel  = dn(2,i)/dn(1,i)
            vvel  = dn(3,i)/dn(1,i)
            wvel  = dn(4,i)/dn(1,i)
            e     = dn(5,i)/dn(1,i)
            temp  = (e - 0.5d0*(uvel*uvel+vvel*vvel+wvel*wvel))
     $            / fac

            if(ro.lt.0.1d0) then
                write(6,*) 'ro at',i,'is',dn(1,i)
                dn(1,i) = 0.1d0
                dn(2,i) = dn(1,i) * uvel
                dn(3,i) = dn(1,i) * vvel
                dn(4,i) = dn(1,i) * wvel
                dn(5,i) = dn(1,i) * e
            endif

            if(temp.lt.0.1d0) then
                write(6,*) 'temp at',i,'is',temp
                temp = 0.1d0
                e = fac * temp + 0.5 *
     $              (dn(2,i)*dn(2,i)+dn(3,i)*dn(3,i)+dn(4,i)*dn(4,i))
     $            / (dn(1,i)*dn(1,i))
                dn(5,i) = dn(1,i) * e
            endif
        enddo

        endsubroutine
!******************************************************************************
 	subroutine limiter_turb(dn1)
        implicit none
        include "global.h"
        real*8 dn1(ndf_t,nnc)
        integer j

        do  j=1,nnc
        if(dn1(ndf_t,j).lt.0.0) then
            write(*,*) 'turb viscosity at',j,'is',dn1(ndf_t,j)
             dn1(ndf_t,j) = -dn1(ndf_t,j)/100.0d0
!            d(1,j) = 1.0e-10  !3.0e-07 !1.0e-08
         end if
        end do

        endsubroutine

