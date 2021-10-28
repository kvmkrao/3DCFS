        subroutine probe(xn,dn)

        implicit none
        include "global.h"

        integer i,j,idf,isd
        real*8  xn(nsd,nnc)
        real*8  dn(ndf,nnc),d2(7,nprobe)
        integer ibase,i1,i2
	real*8 cv,R,gma,velc
        character*15 FileName,FN(nprobe)

        FileName='PROBE/probe.00'
        ibase=ichar('0')
        cv = 1.0d0/((1.4d0)*(0.4d0))
	gma= 1.40d0
        R    = 1.0d0/gma

!	do i=1,nnc
!                  d2(1,i) = dn(2,i)/dn(1,i) ! u
!                  d2(2,i) = dn(3,i)/dn(1,i) ! v
!                  d2(3,i) = dn(4,i)/dn(1,i) ! w
!                  d2(4,i) = dn(1,i)         !rho
!                  d2(5,i) =(dn(5,i)/dn(1,i) !temp
!     $          - (d2(1,i)**2+d2(2,i)**2+d2(3,i)**2) *
!     $             0.5e0)/cv
!		  velc = dsqrt(d2(1,i)**2+d2(2,i)**2+d2(3,i)**2)
!		  d2(6,i) = velc/sqrt(d2(5,i))
!		  d2(7,i) = dn(1,i)*R* d2(5,i)
!	end do


        if (its .eq. 0) then
        do i=1, nprobe
           if (my_pe .eq. probe_pe(i))
     .    write(50,102)my_pe,i,node_pe(i),(xn(isd,node_pe(i)),isd=1,nsd)
        
        end do

        else

        do i=1, nprobe
           if (my_pe .eq. probe_pe(i)) then 
          
          i2=i/10
          i1=(i-(i2)*10)+ibase
          i2=i2+ibase
          FileName(13:13)=char(i2)
          FileName(14:14)=char(i1)
          FN(i)=FileName
          open(unit=61,file=FN(i),status='unknown',access='append')

		   j = node_pe(i)
                  d2(1,i) = dn(2,j)/dn(1,j) ! u
                  d2(2,i) = dn(3,j)/dn(1,j) ! v
                  d2(3,i) = dn(4,j)/dn(1,j) ! w
                  d2(4,i) = dn(1,j)         !rho
                  d2(5,i) =(dn(5,j)/dn(1,j) !temp
     $          - (d2(1,i)**2+d2(2,i)**2+d2(3,i)**2) *
     $             0.5e0)/cv
                  velc = dsqrt(d2(1,i)**2+d2(2,i)**2+d2(3,i)**2)
                  d2(6,i) = velc/sqrt(d2(5,i))
                  d2(7,i) = dn(1,j)*R* d2(5,i)

         write(61,101)its,tt,(d2(idf,i),idf=1,7) !,flowrate(i)
        end if
        end do
        end if
 101    format(i7,2x,8e18.10)
 102    format(i3,2x,i7,5x,4e19.10)

        return
        end

!***********************************************************************
        subroutine massflowrate(ien,d1)

        implicit none
        include "global.h"

        integer ien(nen,nec)
        integer ipn,iely,ielz,i,iel,n1,n2,n3,n4,n5,n6
        real*8  d1(ndf,nnl)
        real*8 r,cv,dm

        r = 1.0d0/gamma
        cv = r/(gamma-1.0d0)

        do ipn=1,nprobe
            flowrate_pe(ipn) = 0.0d0
!!            node=node_pe(ipn)         ! probe node
!!            rho = d1(1,node)
!!            u   = d1(2,node)
!!            v   = d1(3,node)
!!            w   = d1(4,node)
!!            e   = d1(5,node)
!!            temp= (e - 0.5*(u**2+v**2+w**2))/cv
!!            mach= dsqrt((u**2+v**2+w**2)/temp)

!
            do iely=1,nely
                do ielz=1,nelz
                    i = (ipn-1)*nely*nelz+(iely-1)*nelz+ielz
                    if(my_pe.eq.cel_pe(i)) then
                        iel = el_pe(i)
                        n1 = ien(1,iel)
                        n2 = ien(2,iel)
                        n3 = ien(3,iel)
                        n4 = ien(4,iel)
                        n5 = ien(5,iel)
                        n6 = ien(6,iel)

                        dm = msh(1,i)*d1(2,n1)
     $                     + msh(2,i)*d1(2,n2)
     $                     + msh(3,i)*d1(2,n3)
     $                     + msh(4,i)*d1(2,n4)
     $                     + msh(5,i)*d1(2,n5)
     $                     + msh(6,i)*d1(2,n6)

                        flowrate_pe(ipn) = flowrate_pe(ipn) 
     $                                   + dm*dy(i)*dz(i)
                    endif
                enddo
            enddo
        enddo

        endsubroutine
