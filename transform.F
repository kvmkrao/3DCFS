        subroutine transform(d,d_trans,key)

!.......key=1===> go from d to d_trans
!.......key=2===> go from d_trans to d
!.......key=3===> rhoo, Entropy, Pressure, Local Mach no.

        implicit none
        include 'global.h'

        real*8  d(ndf,nnc),d_trans(ndf,nnc)
        real*8  rhoo,u,v,w,temp,pres,c_sound,rmach_loc,roi,entropy
        real*8  cv,gas_c,ub

        integer key,i

        gas_c    = 1.0d0/gamma
        cv       = gas_c/(gamma-1.0d0)

        if(key.eq.1)then
            do i=1,nnc
                d_trans(1,i) = d(2,i)/d(1,i)
                d_trans(2,i) = d(3,i)/d(1,i)
                d_trans(3,i) = d(4,i)/d(1,i)
                d_trans(4,i) = d(1,i)
                d_trans(5,i) =(d(5,i)/d(1,i) 
     $          -     ((d(2,i)/d(1,i))**2+(d(3,i)/d(1,i))**2
     $                +(d(4,i)/d(1,i))**2) * 0.5e0)/cv
            end do
        endif

        if(key.eq.2)then
            do i=1,nnc
                d(1,i) = d_trans(4,i)
                d(2,i) = d_trans(1,i)*d_trans(4,i)
                d(3,i) = d_trans(2,i)*d_trans(4,i)
                d(4,i) = d_trans(3,i)*d_trans(4,i)
                d(5,i) = d_trans(4,i)*(d_trans(5,i)*cv +
     $          0.5e0*(d_trans(1,i)**2+d_trans(2,i)**2+d_trans(3,i)**2))
            end do
        endif
      
        if (key.eq.3)then
            do i=1,nnc
                rhoo = d(1,i)
                u   = d(2,i)/rhoo
                v   = d(3,i)/rhoo
                w   = d(4,i)/rhoo
                temp= (d(5,i)/rhoo-(u*u+v*v+w*w)/2.0e0)/cv
                pres= rhoo*gas_c*temp
                c_sound = sqrt(gamma*gas_c*temp)
                rmach_loc = sqrt(u*u+v*v+w*w)/c_sound
                roi = d(5,i) - 0.5e0*(d(2,i)*d(2,i)+d(3,i)*d(3,i)
     $                             +d(4,i)*d(4,i))/d(1,i)
                entropy = log((gamma-1)*roi/d(1,i)**gamma)

                d_trans(1,i) = rhoo
                d_trans(2,i) = temp
                d_trans(3,i) = pres
                d_trans(4,i) = rmach_loc
                d_trans(5,i) = entropy
            end do
        endif
        return
        endsubroutine
