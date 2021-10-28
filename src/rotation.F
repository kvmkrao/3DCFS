        subroutine rotation(d,d_glob,ang)
        implicit none
        include 'global.h'

        real*8  d(ndf,nnc),d_glob(ndf,nnc),slipang(nnc)
        real*8  ang(nnc)
!        real*8  sin_a,cos_a
        integer i
!.......d() = local coordinate (rotated one)
!.......d_glob() = Global X-Y coordinate 
!.......angle(i) = slope in radians at i th node.
!.......if node is NOT on NOSLIP Solid Bnd then slope = 0.0

        do i=1,nnc
!            cos_a       = cos(ang(i))
!            sin_a       = sin(ang(i))
            d_glob(1,i) = d(1,i)
            d_glob(4,i) = d(4,i)
            d_glob(5,i) = d(5,i)
            d_glob(2,i) = d(2,i) !cos_a*d(2,i)-sin_a*d(3,i)
            d_glob(3,i) = d(3,i) !sin_a*d(2,i)+cos_a*d(3,i)
        end do
        return
        endsubroutine
