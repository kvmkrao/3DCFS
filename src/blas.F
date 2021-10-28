        subroutine saxpy(n,a,x,incx,y,incy)
        implicit none
        integer n,incx,incy
        real*8 a
        real*8 x(1),y(1)
        integer i

        do i = 1,n
                y(i) = y(i) + a*x(i)
        enddo

        return

        end

        function sdot(n,x,incx,y,incy)
        integer n,incx,incy
        real*8 sdot
        real*8 x(1),y(1)

        sdot = 0.0

        do i = 1,n
                sdot = sdot + x(i)*y(i)
        enddo

        return
        end



        subroutine equal(f,d,m)
        implicit none
        integer i,m
        real* 8 f(m),d(m)

        do i=1,m
                d(i) = f(i)
        enddo

        return
        end


         subroutine equal1(f,d,m)
        implicit none
        integer i,m
        integer  f(m),d(m)
                                      
                                                                        
        do i=1,m
                d(i) = f(i)
                write(44,*)d(i),f(i),i
        enddo
                                                                                                                   
        return
        end

        subroutine fclear(f,m)
        implicit none
        integer i,m
        real* 8 f(m)

        do i=1,m
                f(i) = 0.0
        enddo

        return
        end

!**********************************************************************
        subroutine mkblk(ien,block,selm)

        implicit none
        include "global.h"

        integer ien(nen,nec)
        integer iel,i,node,j,k,lt1,lt2
        real*8  selm(nee,nee,nec), block(ndf,ndf,nnl)

        do iel = 1,nec
            do i = 1,nen
                node = ien(i,iel)
                do j = 1,ndf
                    do k = 1,ndf
                        lt1 = (i-1)*ndf + k
                        lt2 = (i-1)*ndf + j
                        block(k,j,node)=block(k,j,node)+selm(lt1,lt2,iel)
                    end do
                end do
            end do
        end do
!!        do iel=1,nec
!!            do i=1,nen
!!                node=ien(i,iel)
!!                write(80+my_pe,*) ((block(k,j,node),k=1,ndf),j=1,ndf)
!!            enddo
!!        enddo
        return
        end

!**********************************************************************
!!!        subroutine facblk(id,block)
!!!
!!!        IMPLICIT REAL*8 (A-H,O-Z)
!!!        common /cinfo/ ndof, nsd, numnp, numel
!!!
!!!        integer id(ndof,numnp)
!!!        real*8 block(ndof,ndof,numnp)
!!!
!!!        do i = 1,numnp
!!!            do j = 1,ndof
!!!
!!!                if(id(j,i).eq.0) then
!!!                    do k = 1,ndof
!!!                        block(j,k,i) = 0.0
!!!                        block(k,j,i) = 0.0
!!!                    end do
!!!                    block(j,j,i) = 1.0
!!!                end if
!!!            end do
!!!
!!!            do j = 1,ndof-1
!!!                do k = j+1,ndof
!!!                    block(k,j,i) = block(k,j,i)/block(j,j,i)
!!!                        do l = j+1,ndof
!!!                            block(k,l,i) = block(k,l,i) - 
!!!     $                          block(k,j,i)*block(j,l,i)
!!!                        end do
!!!                end do
!!!            end do
!!!
!!!            do j = ndof-1,1,-1
!!!                do k = j+1,ndof
!!!                    block(j,k,i) = block(j,k,i)/block(j,j,i)
!!!                end do
!!!            end do
!!!        end do
!!!
!!!        return
!!!        end
