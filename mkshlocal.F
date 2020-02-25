        subroutine mkshlocal
        implicit none
        include "global.h"
        integer iq,i,j,ii
        real*8 sm1,sp1,tm1,tp1,nm1,np1,rm1,rp1
	real*8 t1,t2,t3,t4,t5,t6,t7,t8,t9

C...SET UP FOR Semi-Discrete and Spacetime-Jump terms
C...TETS
        if (elem_type.eq.elem_tetr_sd .or. elem_type.eq.elem_tetr_st) then
                call quad3d4n(iquad, nquad3, xq3, wq3, nsdpad, nquadpad)
                do iq = 1,nquad3
                sq3(0,1,iq) = xq3(1,iq)
                sq3(0,2,iq) = xq3(2,iq)
                sq3(0,3,iq) = xq3(3,iq)
                sq3(0,4,iq) = 1 - xq3(1,iq) - xq3(2,iq) - xq3(3,iq)



c               local first derivatives
                sq3(1,1,iq) =  1.0
                sq3(1,2,iq) =  0.0
                sq3(1,3,iq) =  0.0
                sq3(1,4,iq) = -1.0

                sq3(2,1,iq) =  0.0
                sq3(2,2,iq) =  1.0
                sq3(2,3,iq) =  0.0
                sq3(2,4,iq) = -1.0

                sq3(3,1,iq) =  0.0
                sq3(3,2,iq) =  0.0
                sq3(3,3,iq) =  1.0
                sq3(3,4,iq) = -1.0

                enddo

C Prism
        
        elseif(elem_type.eq.elem_wedg_sd)then
!	call quad3d15n(iquad,nquad3,xq3,wq3,nsdpad,nquadpad)
        call quad3d6n(iquad,nquad3,xq3,wq3,nsdpad,nquadpad)
                do iq=1,nquad3
                sq3(0,1,iq) = (1.0 - xq3(1,iq) - xq3(2,iq))
     $                                  * (1.0 - xq3(3,iq)) / 2.0
                sq3(0,2,iq) = xq3(1,iq) * (1.0 - xq3(3,iq)) / 2.0
                sq3(0,3,iq) = xq3(2,iq) * (1.0 - xq3(3,iq)) / 2.0
                sq3(0,4,iq) = (1.0 - xq3(1,iq) - xq3(2,iq))
     $                                  * (1.0 + xq3(3,iq)) / 2.0
                sq3(0,5,iq) = xq3(1,iq) * (1.0 + xq3(3,iq)) / 2.0
                sq3(0,6,iq) = xq3(2,iq) * (1.0 + xq3(3,iq)) / 2.0
                                                                                                     
c       Local Derivatives of shape functions
                                                                                                                             
                sq3(1,1,iq) = -1.0 * (1.0 - xq3(3,iq)) / 2.0
                sq3(1,2,iq) = +1.0 * (1.0 - xq3(3,iq)) / 2.0
                sq3(1,3,iq) = +0.0
                sq3(1,4,iq) = -1.0 * (1.0 + xq3(3,iq)) / 2.0
                sq3(1,5,iq) = +1.0 * (1.0 + xq3(3,iq)) / 2.0
                sq3(1,6,iq) = +0.0
                                                                                                                             
                sq3(2,1,iq) = -1.0 * (1.0 - xq3(3,iq)) / 2.0
                sq3(2,2,iq) = +0.0
                sq3(2,3,iq) = +1.0 * (1.0 - xq3(3,iq)) / 2.0
                sq3(2,4,iq) = -1.0 * (1.0 + xq3(3,iq)) / 2.0
                sq3(2,5,iq) = +0.0
                sq3(2,6,iq) = +1.0 * (1.0 + xq3(3,iq)) / 2.0
                                                                                                                             
                sq3(3,1,iq) = -1.0 * (1.0 - xq3(1,iq)-xq3(2,iq))/2.0
                sq3(3,2,iq) = -1.0 * xq3(1,iq) / 2.0
                sq3(3,3,iq) = -1.0 * xq3(2,iq) / 2.0
                sq3(3,4,iq) = +1.0 * (1.0 - xq3(1,iq)-xq3(2,iq))/2.0
                sq3(3,5,iq) = +1.0 * xq3(1,iq) / 2.0
                sq3(3,6,iq) = +1.0 * xq3(2,iq) / 2.0

                end do


        elseif(elem_type.eq.elem_wedg_quad_sd)then
	call quad3d15n(iquad,nquad3,xq3,wq3,nsdpad,nquadpad)
!        call quad3d6n(iquad,nquad3,xq3,wq3,nsdpad,nquadpad)

!	write(6,*) "Hi"
        do iq=1,nquad3

         t1 = (2.0*xq3(1,iq)-1.0)/2.0d0
         t2 = (2.0*xq3(2,iq)-1.0)/2.0d0
         t3 = (1.0-xq3(1,iq)-xq3(2,iq))/2.0d0
         t4 = (1.0-xq3(3,iq))
         t5 = (1.0-xq3(3,iq)*xq3(3,iq))
         t6 = (1.0+xq3(3,iq))
	 t7 = (4.0*xq3(1,iq)-1.0)/2.0d0
	 t8 = (-3.0+4.0*xq3(1,iq)+4.0*xq3(2,iq))/2.0d0
	 t9 = (4.0*xq3(2,iq)-1.0)/2.0


	 sq3(0,1,iq) = xq3(1,iq)*t1*t4 - 0.5*xq3(1,iq)*t5
	 sq3(0,2,iq) = xq3(2,iq)*t2*t4 - 0.5*xq3(2,iq)*t5
	 sq3(0,3,iq) = (1.0-2.0*xq3(1,iq)-2.0*xq3(2,iq))*t3*t4 - t3*t5
         sq3(0,4,iq) = 2.0*xq3(1,iq)*xq3(2,iq)*t4
	 sq3(0,5,iq) = 4.0*xq3(1,iq)*t3*t4
	 sq3(0,6,iq) = 4.0*xq3(2,iq)*t3*t4
         sq3(0,7,iq) = xq3(1,iq)*t1*t6 - 0.5*xq3(1,iq)*t5
	 sq3(0,8,iq) = xq3(2,iq)*t2*t6 - 0.5*xq3(2,iq)*t5
	 sq3(0,9,iq) = (1.0-2.0*xq3(1,iq)-2.0*xq3(2,iq))*t3*t6 - t3*t5
         sq3(0,10,iq) = 2.0*xq3(1,iq)*xq3(2,iq)*t6
	 sq3(0,11,iq) = 4.0*xq3(1,iq)*t3*t6
	 sq3(0,12,iq) = 4.0*xq3(2,iq)*t3*t6
         sq3(0,13,iq) = xq3(1,iq)*t5
	 sq3(0,14,iq) = xq3(2,iq)*t5
	 sq3(0,15,iq) = 2.0*t3*t5

!	,xi 
         sq3(1,1,iq) = t4*t7 - 0.5*t5
	 sq3(1,2,iq) = 0.0d0
	 sq3(1,3,iq) = t4*t8 + 0.5*t5
         sq3(1,4,iq) =  2.0*xq3(2,iq)*t4
	 sq3(1,5,iq) = (2.0-4.0*xq3(1,iq)-2.0*xq3(2,iq))*t4
	 sq3(1,6,iq) = -2.0*xq3(2,iq)*t4
         sq3(1,7,iq) = t6*t7 - 0.5*t5
	 sq3(1,8,iq) = 0.0d0
	 sq3(1,9,iq) = t6*t8 + 0.5*t5
         sq3(1,10,iq) =  2.0*xq3(2,iq)*t6
	 sq3(1,11,iq) = (2.0-4.0*xq3(1,iq)-2.0*xq3(2,iq))*t6
	 sq3(1,12,iq) = -2.0*xq3(2,iq)*t6
         sq3(1,13,iq) =  t5
	 sq3(1,14,iq) = 0.0d0
	 sq3(1,15,iq) = -t5

!	,eta	
         sq3(2,1,iq) = 0.0d0
	 sq3(2,2,iq) = t4*t9 - 0.5*t5 
	 sq3(2,3,iq) = t4*t8 + 0.5*t5
         sq3(2,4,iq) =  2.0*xq3(1,iq)*t4
	 sq3(2,5,iq) = -2.0*xq3(1,iq)*t4
	 sq3(2,6,iq) = (2.0-4.0*xq3(2,iq)-2.0*xq3(1,iq))*t4
         sq3(2,7,iq) = 0.0d0
	 sq3(2,8,iq) = t6*t9 - 0.5*t5
	 sq3(2,9,iq) = t6*t8 + 0.5*t5
         sq3(2,10,iq) =  2.0*xq3(1,iq)*t6
	 sq3(2,11,iq) = -2.0*xq3(1,iq)*t6
	 sq3(2,12,iq) = (2.0-4.0*xq3(2,iq)-2.0*xq3(1,iq))*t6 
         sq3(2,13,iq) = 0.0d0
	 sq3(2,14,iq) =  t5
	 sq3(2,15,iq) = -t5

! 	,zeta 
         sq3(3,1,iq) = -xq3(1,iq)*t1 + xq3(1,iq)*xq3(3,iq)
	 sq3(3,2,iq) = -xq3(2,iq)*t2 + xq3(2,iq)*xq3(3,iq)
	 sq3(3,3,iq) = -(1.0-2.0*xq3(1,iq)-2.0*xq3(2,iq))*t3
     $				      +2.0*xq3(3,iq)*t3
         sq3(3,4,iq) = -2.0*xq3(1,iq)*xq3(2,iq)
	 sq3(3,5,iq) = -4.0*xq3(1,iq)*t3
	 sq3(3,6,iq) = -4.0*xq3(2,iq)*t3
         sq3(3,7,iq) = xq3(1,iq)*t1 + xq3(1,iq)*xq3(3,iq)
	 sq3(3,8,iq) = xq3(2,iq)*t2 + xq3(2,iq)*xq3(3,iq)
	 sq3(3,9,iq) = (1.0-2.0*xq3(1,iq)-2.0*xq3(2,iq))*t3 
     $				     +2.0*xq3(3,iq)*t3
         sq3(3,10,iq) = 2.0*xq3(1,iq)*xq3(2,iq)
	 sq3(3,11,iq) = 4.0*xq3(1,iq)*t3
	 sq3(3,12,iq) = 4.0*xq3(2,iq)*t3
         sq3(3,13,iq) = -2.0*xq3(1,iq)*xq3(3,iq)
	 sq3(3,14,iq) = -2.0*xq3(2,iq)*xq3(3,iq)
	 sq3(3,15,iq) = -4.0*xq3(3,iq)*t3
	end do

C...BRICKS
        else if (elem_type.eq.elem_hexa_sd .or. elem_type.eq.elem_hexa_st) then
        call quad3d8n(2, nquad3, xq3, wq3, nsdpad, nquadpad)
                do iq = 1,nquad3
                sq3(0,1,iq) = (1 - xq3(1,iq))
     &                     * (1 - xq3(2,iq)) * (1 - xq3(3,iq)) / 8
                sq3(0,2,iq) = (1 + xq3(1,iq))
     &                     * (1 - xq3(2,iq)) * (1 - xq3(3,iq)) / 8
                sq3(0,3,iq) = (1 + xq3(1,iq))
     &                     * (1 + xq3(2,iq)) * (1 - xq3(3,iq)) / 8
                sq3(0,4,iq) = (1 - xq3(1,iq))
     &                     * (1 + xq3(2,iq)) * (1 - xq3(3,iq)) / 8
                sq3(0,5,iq) = (1 - xq3(1,iq))
     &                     * (1 - xq3(2,iq)) * (1 + xq3(3,iq)) / 8
                sq3(0,6,iq) = (1 + xq3(1,iq))
     &                     * (1 - xq3(2,iq)) * (1 + xq3(3,iq)) / 8
                sq3(0,7,iq) = (1 + xq3(1,iq))
     &                     * (1 + xq3(2,iq)) * (1 + xq3(3,iq)) / 8
                sq3(0,8,iq) = (1 - xq3(1,iq))
     &                     * (1 + xq3(2,iq)) * (1 + xq3(3,iq)) / 8

c               local first derivatives
                sq3(1,1,iq) = - (1 - xq3(2,iq)) * (1 - xq3(3,iq)) / 8
                sq3(1,2,iq) = + (1 - xq3(2,iq)) * (1 - xq3(3,iq)) / 8
                sq3(1,3,iq) = + (1 + xq3(2,iq)) * (1 - xq3(3,iq)) / 8
                sq3(1,4,iq) = - (1 + xq3(2,iq)) * (1 - xq3(3,iq)) / 8
                sq3(1,5,iq) = - (1 - xq3(2,iq)) * (1 + xq3(3,iq)) / 8
                sq3(1,6,iq) = + (1 - xq3(2,iq)) * (1 + xq3(3,iq)) / 8
                sq3(1,7,iq) = + (1 + xq3(2,iq)) * (1 + xq3(3,iq)) / 8
                sq3(1,8,iq) = - (1 + xq3(2,iq)) * (1 + xq3(3,iq)) / 8
        
                sq3(2,1,iq) = - (1 - xq3(1,iq)) * (1 - xq3(3,iq)) / 8
                sq3(2,2,iq) = - (1 + xq3(1,iq)) * (1 - xq3(3,iq)) / 8
                sq3(2,3,iq) = + (1 + xq3(1,iq)) * (1 - xq3(3,iq)) / 8
                sq3(2,4,iq) = + (1 - xq3(1,iq)) * (1 - xq3(3,iq)) / 8
                sq3(2,5,iq) = - (1 - xq3(1,iq)) * (1 + xq3(3,iq)) / 8
                sq3(2,6,iq) = - (1 + xq3(1,iq)) * (1 + xq3(3,iq)) / 8
                sq3(2,7,iq) = + (1 + xq3(1,iq)) * (1 + xq3(3,iq)) / 8
                sq3(2,8,iq) = + (1 - xq3(1,iq)) * (1 + xq3(3,iq)) / 8
        
                sq3(3,1,iq) = - (1 - xq3(1,iq)) * (1 - xq3(2,iq)) / 8
                sq3(3,2,iq) = - (1 + xq3(1,iq)) * (1 - xq3(2,iq)) / 8
                sq3(3,3,iq) = - (1 + xq3(1,iq)) * (1 + xq3(2,iq)) / 8
                sq3(3,4,iq) = - (1 - xq3(1,iq)) * (1 + xq3(2,iq)) / 8
                sq3(3,5,iq) = + (1 - xq3(1,iq)) * (1 - xq3(2,iq)) / 8
                sq3(3,6,iq) = + (1 + xq3(1,iq)) * (1 - xq3(2,iq)) / 8
                sq3(3,7,iq) = + (1 + xq3(1,iq)) * (1 + xq3(2,iq)) / 8
                sq3(3,8,iq) = + (1 - xq3(1,iq)) * (1 + xq3(2,iq)) / 8
        
                enddo
C...TRIANGLES
        else if (elem_type.eq.elem_tria_sd) then
        call quad2d3n(iquad, nquad3, xq3, wq3, nsdpad, nquadpad)
                do iq = 1,nquad3
                        sq3(0,1,iq) = xq3(1,iq)
                        sq3(0,2,iq) = xq3(2,iq)
                        sq3(0,3,iq) = 1.0 -xq3(1,iq) -xq3(2,iq)
                enddo
C...QUADS
        else if (elem_type.eq.elem_quad_sd) then
        call quad2d4n(iquad, nquad3, xq3, wq3, nsdpad, nquadpad)
                  do iq = 1,nquad3
                         sm1     = (1-xq3(1,iq))
                         sp1     = (1+xq3(1,iq))
                         tm1     = (1-xq3(2,iq))
                         tp1     = (1+xq3(2,iq))
C...SHAPE FUNCTIONS VALUES
                        sq3(0,1,iq) = 0.25*sm1*tm1
                        sq3(0,2,iq) = 0.25*sp1*tm1
                        sq3(0,3,iq) = 0.25*sp1*tp1
                        sq3(0,4,iq) = 0.25*sm1*tp1

C...LOCAL DERIVATIVES
                        sq3(1,1,iq) = - 0.25*tm1
                        sq3(1,2,iq) = + 0.25*tm1
                        sq3(1,3,iq) = + 0.25*tp1
                        sq3(1,4,iq) = - 0.25*tp1

                        sq3(2,1,iq) = - 0.25*sm1
                        sq3(2,2,iq) = - 0.25*sp1
                        sq3(2,3,iq) = + 0.25*sp1
                        sq3(2,4,iq) = + 0.25*sm1
                enddo

        end if


C...SETUP SPACETIME SHAPEFUNCTIONS
C...TETS
        if (elem_type.eq.elem_tetr_st) then
                call quad4d4n(iquad,nquad4,xq4,wq4,4,16)

                do iq = 1,nquad4
                        rm1     = (1-xq4(4,iq))
                        rp1     = (1+xq4(4,iq))
C...SHAPE FUNCTIONS VALUES
                        sq4(0,1,iq) = 0.5*xq4(1,iq)*rm1
                        sq4(0,2,iq) = 0.5*xq4(2,iq)*rm1
                        sq4(0,3,iq) = 0.5*xq4(3,iq)*rm1
                        sq4(0,4,iq) =
     .                     0.5*(1.0-xq4(1,iq)-xq4(2,iq)-xq4(3,iq))*rm1
                        sq4(0,5,iq) = 0.5*xq4(1,iq)*rp1
                        sq4(0,6,iq) = 0.5*xq4(2,iq)*rp1
                        sq4(0,7,iq) = 0.5*xq4(3,iq)*rp1
                        sq4(0,8,iq) =
     .                     0.5*(1.0-xq4(1,iq)-xq4(2,iq)-xq4(3,iq))*rp1

C...LOCAL DERIVATIVES
                        sq4(1,1,iq) =  0.5*rm1
                        sq4(1,2,iq) =  0.0
                        sq4(1,3,iq) =  0.0
                        sq4(1,4,iq) = -0.5*rm1
                        sq4(1,5,iq) =  0.5*rp1
                        sq4(1,6,iq) =  0.0
                        sq4(1,7,iq) =  0.0
                        sq4(1,8,iq) = -0.5*rp1

                        sq4(2,1,iq) =  0.0
                        sq4(2,2,iq) =  0.5*rm1
                        sq4(2,3,iq) =  0.0
                        sq4(2,4,iq) = -0.5*rm1
                        sq4(2,5,iq) =  0.0
                        sq4(2,6,iq) =  0.5*rp1
                        sq4(2,7,iq) =  0.0
                        sq4(2,8,iq) = -0.5*rp1

                        sq4(3,1,iq) =  0.0
                        sq4(3,2,iq) =  0.0
                        sq4(3,3,iq) =  0.5*rm1
                        sq4(3,4,iq) = -0.5*rm1
                        sq4(3,5,iq) =  0.0
                        sq4(3,6,iq) =  0.0
                        sq4(3,7,iq) =  0.5*rp1
                        sq4(3,8,iq) = -0.5*rp1

                        sq4(4,1,iq) = -0.5*xq4(1,iq)
                        sq4(4,2,iq) = -0.5*xq4(2,iq)
                        sq4(4,3,iq) = -0.5*xq4(3,iq)
                        sq4(4,4,iq) = -0.5*(1.0-xq4(1,iq)-xq4(2,iq)-xq4(3,iq))
                        sq4(4,5,iq) =  0.5*xq4(1,iq)
                        sq4(4,6,iq) =  0.5*xq4(2,iq)
                        sq4(4,7,iq) =  0.5*xq4(3,iq)
                        sq4(4,8,iq) =  0.5*(1.0-xq4(1,iq)-xq4(2,iq)-xq4(3,iq))

                enddo

C...BRICKS
        else if (elem_type.eq.elem_hexa_st) then
                call quad4d8n(iquad,nquad4,xq4,wq4,4,16)

                do iq = 1,nquad4
                         sm1     = (1-xq4(1,iq))
                         sp1     = (1+xq4(1,iq))
                         tm1     = (1-xq4(2,iq))
                         tp1     = (1+xq4(2,iq))
                         nm1     = (1-xq4(3,iq))
                         np1     = (1+xq4(3,iq))
                         rm1     = (1-xq4(4,iq))
                         rp1     = (1+xq4(4,iq))
C...SHAPE FUNCTIONS VALUES
                        sq4(0, 1,iq) = 0.0625*sm1*tm1*nm1*rm1
                        sq4(0, 2,iq) = 0.0625*sp1*tm1*nm1*rm1
                        sq4(0, 3,iq) = 0.0625*sp1*tp1*nm1*rm1
                        sq4(0, 4,iq) = 0.0625*sm1*tp1*nm1*rm1
                        sq4(0, 5,iq) = 0.0625*sm1*tm1*np1*rm1
                        sq4(0, 6,iq) = 0.0625*sp1*tm1*np1*rm1
                        sq4(0, 7,iq) = 0.0625*sp1*tp1*np1*rm1
                        sq4(0, 8,iq) = 0.0625*sm1*tp1*np1*rm1
                        sq4(0, 9,iq) = 0.0625*sm1*tm1*nm1*rp1
                        sq4(0,10,iq) = 0.0625*sp1*tm1*nm1*rp1
                        sq4(0,11,iq) = 0.0625*sp1*tp1*nm1*rp1
                        sq4(0,12,iq) = 0.0625*sm1*tp1*nm1*rp1
                        sq4(0,13,iq) = 0.0625*sm1*tm1*np1*rp1
                        sq4(0,14,iq) = 0.0625*sp1*tm1*np1*rp1
                        sq4(0,15,iq) = 0.0625*sp1*tp1*np1*rp1
                        sq4(0,16,iq) = 0.0625*sm1*tp1*np1*rp1

C...LOCAL DERIVATIVES
                        sq4(1, 1,iq) = - 0.0625*tm1*nm1*rm1
                        sq4(1, 2,iq) = + 0.0625*tm1*nm1*rm1
                        sq4(1, 3,iq) = + 0.0625*tp1*nm1*rm1
                        sq4(1, 4,iq) = - 0.0625*tp1*nm1*rm1
                        sq4(1, 5,iq) = - 0.0625*tm1*np1*rm1
                        sq4(1, 6,iq) = + 0.0625*tm1*np1*rm1
                        sq4(1, 7,iq) = + 0.0625*tp1*np1*rm1
                        sq4(1, 8,iq) = - 0.0625*tp1*np1*rm1
                        sq4(1, 9,iq) = - 0.0625*tm1*nm1*rp1
                        sq4(1,10,iq) = + 0.0625*tm1*nm1*rp1
                        sq4(1,11,iq) = + 0.0625*tp1*nm1*rp1
                        sq4(1,12,iq) = - 0.0625*tp1*nm1*rp1
                        sq4(1,13,iq) = - 0.0625*tm1*np1*rp1
                        sq4(1,14,iq) = + 0.0625*tm1*np1*rp1
                        sq4(1,15,iq) = + 0.0625*tp1*np1*rp1
                        sq4(1,16,iq) = - 0.0625*tp1*np1*rp1

                        sq4(2, 1,iq) = - 0.0625*sm1*nm1*rm1
                        sq4(2, 2,iq) = - 0.0625*sp1*nm1*rm1
                        sq4(2, 3,iq) = + 0.0625*sp1*nm1*rm1
                        sq4(2, 4,iq) = + 0.0625*sm1*nm1*rm1
                        sq4(2, 5,iq) = - 0.0625*sm1*np1*rm1
                        sq4(2, 6,iq) = - 0.0625*sp1*np1*rm1
                        sq4(2, 7,iq) = + 0.0625*sp1*np1*rm1
                        sq4(2, 8,iq) = + 0.0625*sm1*np1*rm1
                        sq4(2, 9,iq) = - 0.0625*sm1*nm1*rp1
                        sq4(2,10,iq) = - 0.0625*sp1*nm1*rp1
                        sq4(2,11,iq) = + 0.0625*sp1*nm1*rp1
                        sq4(2,12,iq) = + 0.0625*sm1*nm1*rp1
                        sq4(2,13,iq) = - 0.0625*sm1*np1*rp1
                        sq4(2,14,iq) = - 0.0625*sp1*np1*rp1
                        sq4(2,15,iq) = + 0.0625*sp1*np1*rp1
                        sq4(2,16,iq) = + 0.0625*sm1*np1*rp1

                        sq4(3, 1,iq) = - 0.0625*sm1*tm1*rm1
                        sq4(3, 2,iq) = - 0.0625*sp1*tm1*rm1
                        sq4(3, 3,iq) = - 0.0625*sp1*tp1*rm1
                        sq4(3, 4,iq) = - 0.0625*sm1*tp1*rm1
                        sq4(3, 5,iq) = + 0.0625*sm1*tm1*rm1
                        sq4(3, 6,iq) = + 0.0625*sp1*tm1*rm1
                        sq4(3, 7,iq) = + 0.0625*sp1*tp1*rm1
                        sq4(3, 8,iq) = + 0.0625*sm1*tp1*rm1
                        sq4(3, 9,iq) = - 0.0625*sm1*tm1*rp1
                        sq4(3,10,iq) = - 0.0625*sp1*tm1*rp1
                        sq4(3,11,iq) = - 0.0625*sp1*tp1*rp1
                        sq4(3,12,iq) = - 0.0625*sm1*tp1*rp1
                        sq4(3,13,iq) = + 0.0625*sm1*tm1*rp1
                        sq4(3,14,iq) = + 0.0625*sp1*tm1*rp1
                        sq4(3,15,iq) = + 0.0625*sp1*tp1*rp1
                        sq4(3,16,iq) = + 0.0625*sm1*tp1*rp1

                        sq4(4, 1,iq) = - 0.0625*sm1*tm1*nm1
                        sq4(4, 2,iq) = - 0.0625*sp1*tm1*nm1
                        sq4(4, 3,iq) = - 0.0625*sp1*tp1*nm1
                        sq4(4, 4,iq) = - 0.0625*sm1*tp1*nm1
                        sq4(4, 5,iq) = - 0.0625*sm1*tm1*np1
                        sq4(4, 6,iq) = - 0.0625*sp1*tm1*np1
                        sq4(4, 7,iq) = - 0.0625*sp1*tp1*np1
                        sq4(4, 8,iq) = - 0.0625*sm1*tp1*np1
                        sq4(4, 9,iq) = + 0.0625*sm1*tm1*nm1
                        sq4(4,10,iq) = + 0.0625*sp1*tm1*nm1
                        sq4(4,11,iq) = + 0.0625*sp1*tp1*nm1
                        sq4(4,12,iq) = + 0.0625*sm1*tp1*nm1
                        sq4(4,13,iq) = + 0.0625*sm1*tm1*np1
                        sq4(4,14,iq) = + 0.0625*sp1*tm1*np1
                        sq4(4,15,iq) = + 0.0625*sp1*tp1*np1
                        sq4(4,16,iq) = + 0.0625*sm1*tp1*np1

                enddo

        endif   

c........SHAPE FUNCTION FOR WEDGE ELEMENT
        
c       if(elem_type.eq.elem_wedg_sd)then
c               call quad3d6n(iquad,nquad3,xq3,wq3,nsdpad,nquadpad)
c               do iq=1,nquad3
c               sq3(0,1,iq) = (1.0 - xq3(1,iq) - xq3(2,iq))
c     $                                  * (1.0 - xq3(3,iq)) / 2.0
c               sq3(0,2,iq) = xq3(1,iq) * (1.0 - xq3(3,iq)) / 2.0
c               sq3(0,3,iq) = xq3(2,iq) * (1.0 - xq3(3,iq)) / 2.0
c               sq3(0,4,iq) = (1.0 - xq3(1,iq) - xq3(2,iq))
c     $                                  * (1.0 + xq3(3,iq)) / 2.0
c                       sq3(0,5,iq) = xq3(1,iq) * (1.0 + xq3(3,iq)) / 2.0
c               sq3(0,6,iq) = xq3(2,iq) * (1.0 + xq3(3,iq)) / 2.0
c
cc      Local Derivatives of shape functions
c
c               sq3(1,1,iq) = -1.0 * (1.0 - xq(3,iq)) / 2.0
c               sq3(1,2,iq) = +1.0 * (1.0 - xq(3,iq)) / 2.0
c                sq3(1,3,iq) = +0.0
c                sq3(1,4,iq) = -1.0 * (1.0 + xq(3,iq)) / 2.0
c               sq3(1,5,iq) = +1.0 * (1.0 + xq(3,iq)) / 2.0
c               sq3(1,6,iq) = +0.0
c
c               sq3(2,1,iq) = -1.0 * (1.0 - xq(3,iq)) / 2.0
c               sq3(2,2,iq) = +0.0
c               sq3(2,3,iq) = +1.0 * (1.0 - xq(3,iq)) / 2.0
c               sq3(2,4,iq) = -1.0 * (1.0 + xq(3,iq)) / 2.0
c               sq3(2,5,iq) = +0.0
c cc                    sq3(2,6,iq) = +1.0 * (1.0 + xq(3,iq)) / 2.0
c
c               sq3(3,1,iq) = -1.0 * (1.0 - xq(1,iq)-xq(2,iq))/2.0
c               sq3(3,2,iq) = -1.0 * xq(1,iq) / 2.0
c               sq3(3,3,iq) = -1.0 * xq(2,iq) / 2.0
c               sq3(3,4,iq) = +1.0 * (1.0 - xq(1,iq)-xq(2,iq))/2.0
c               sq3(3,5,iq) = +1.0 * xq(1,iq) / 2.0
c               sq3(3,6,iq) = +1.0 * xq(2,iq) / 2.0
c               
c               end do          
c

        return
        end
