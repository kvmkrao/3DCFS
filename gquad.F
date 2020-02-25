c       cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       quad3d4n.fcm                                                         c
c       ---------------------------------------------------------------------c
c     quadrature definitions in three dimensional tetrahedral ref domain   c
c       ---------------------------------------------------------------------c
c       930329 - converted from quad2d3n.fcm                                 c
c                rules taken from Hughes Table 3.I.2                         c
c       cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine quad3d4n(iquad, nquad, xquad, wquad, maxnsd, maxnquad)



        real* 8 xquad(maxnsd,maxnquad)
        real* 8 wquad(maxnquad)

        if ((iquad.eq.1).or.(iquad.eq.313131)) then

c               1 point quadrature - code 3.1,3.1,3.1 - old code 1 - precision 2
                nquad = 1

                xquad(1,1) = + 0.250000000000000
                xquad(2,1) = + 0.250000000000000
                xquad(3,1) = + 0.250000000000000

                wquad(1) = 0.166666666666667
                return

        else if ((iquad.eq.2).or.(iquad.eq.323232)) then
c               4 point quadrature - code 3.2,3.2,3.2 - old code 2 - precision 3
                nquad = 4

                xquad(1,1) = + 0.585410200000000
                xquad(2,1) = + 0.138196600000000
                xquad(3,1) = + 0.138196600000000
                xquad(1,2) = + 0.138196600000000
                xquad(2,2) = + 0.585410200000000
                xquad(3,2) = + 0.138196600000000
                xquad(1,3) = + 0.138196600000000
                xquad(2,3) = + 0.138196600000000
                xquad(3,3) = + 0.585410200000000
                xquad(1,4) = + 0.138196600000000
                xquad(2,4) = + 0.138196600000000
                xquad(3,4) = + 0.138196600000000

                wquad(1) = 0.041666666666667
                wquad(2) = 0.041666666666667
                wquad(3) = 0.041666666666667
                wquad(4) = 0.041666666666667
                return

        else if ((iquad.eq.3).or.(iquad.eq.333333)) then
c               5 point quadrature - code 3.3,3.3,3.3 - old code 3 - precision 4
                nquad = 5

                xquad(1,1) = + 0.250000000000000
                xquad(2,1) = + 0.250000000000000
                xquad(3,1) = + 0.250000000000000
                xquad(1,2) = + 0.333333333333333
                xquad(2,2) = + 0.166666666666667
                xquad(3,2) = + 0.166666666666667
                xquad(1,3) = + 0.166666666666667
                xquad(2,3) = + 0.333333333333333
                xquad(3,3) = + 0.166666666666667
                xquad(1,4) = + 0.166666666666667
                xquad(2,4) = + 0.166666666666667
                xquad(3,4) = + 0.333333333333333
                xquad(1,5) = + 0.166666666666667
                xquad(2,5) = + 0.166666666666667
                xquad(3,5) = + 0.166666666666667

                wquad(1) = - 0.133333333333333
                wquad(2) = + 0.075000000000000
                wquad(3) = + 0.075000000000000
                wquad(4) = + 0.075000000000000
                wquad(5) = + 0.075000000000000
                return

        else if ((iquad.eq.4).or.(iquad.eq.222222)) then
c       4 point Lobatto quadrature - element corners - code 2.2,2.2,2.2 - old 4
                nquad = 4

                xquad(1,1) = + 0.0
                xquad(2,1) = + 0.0
                xquad(3,1) = + 0.0
                xquad(1,2) = + 1.0
                xquad(2,2) = + 0.0
                xquad(3,2) = + 0.0
                xquad(1,3) = + 0.0
                xquad(2,3) = + 1.0
                xquad(3,3) = + 0.0
                xquad(1,4) = + 0.0
                xquad(2,4) = + 0.0
                xquad(3,4) = + 1.0
        
                wquad(1) = 0.041666666666667
                wquad(2) = 0.041666666666667
                wquad(3) = 0.041666666666667
                wquad(4) = 0.041666666666667
                return

        else
                call error("quad3d4n: unknown quadrature code ", iquad, .true.)
        end if

        end


c       cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       quad3d15n                                                            c
c       ---------------------------------------------------------------------c
c       quadrature definitions in three dimensional tetrahedral ref domain   c
c       ---------------------------------------------------------------------c
c       930329 - converted from quad3d6n                                     c
c                rules taken from Zienkiewicz Table 9.2                      c
c       cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine quad3d15n(igaus,ngaus,gaus,weis,maxnsd, maxnquad)

        real* 8 gaus(maxnsd,maxnquad)
        real* 8 weis(maxnquad)

        integer igaus,ngaus
!
        goto (1,2,3,4),igaus
!
 1      continue
        ngaus = 1
!
        gaus(1,1) = + 1.0d0/3.0d0
        gaus(2,1) = + 1.0d0/3.0d0
        gaus(3,1) = + 1.0d0/3.0d0
!
        weis(1) = 1.0d0
        return
!
 2      continue
	ngaus = 3

        gaus(1,1) = + 0.5d0
        gaus(2,1) = + 0.5d0
        gaus(3,1) = + 0.0d0

        gaus(1,2) = + 0.0d0
        gaus(2,2) = + 0.5d0
        gaus(3,2) = + 0.5d0

        gaus(1,3) = + 0.5d0
        gaus(2,3) = + 0.0d0
        gaus(3,3) = + 0.5d0

	weis(1)   = + 1.0d0/3.0d0
	weis(2)   = + 1.0d0/3.0d0
	weis(3)   = + 1.0d0/3.0d0

	return

3      continue
        ngaus = 4
!
        gaus(1,1) = + 1.0d0/3.0d0
        gaus(2,1) = + 1.0d0/3.0d0 
        gaus(3,1) = + 1.0d0/3.0d0

        gaus(1,2) = + 0.60d0
        gaus(2,2) = + 0.20d0
        gaus(3,2) = + 0.20d0

        gaus(1,3) = + 0.20d0
        gaus(2,3) = + 0.60d0
        gaus(3,3) = + 0.20d0

        gaus(1,4) = + 0.20d0
        gaus(2,4) = + 0.20d0
        gaus(3,4) = + 0.60d0

	weis(1)   = -27.0d0/48.0d0
	weis(2)   =  25.0d0/48.0d0
	weis(3)   =  25.0d0/48.0d0
	weis(4)   =  25.0d0/48.0d0
        return
!
!.....  lumped at the nodes
!
 4      continue
        ngaus = 7

        gaus(1,1) = 1.0d0/3.0d0
        gaus(2,1) = 1.0d0/3.0d0 
        gaus(3,1) = 1.0d0/3.0d0

        gaus(1,2) = 0.0597158717d0
        gaus(2,2) = 0.4701420641d0
        gaus(3,2) = 0.4701420641d0

        gaus(1,3) = 0.4701420641d0
        gaus(2,3) = 0.0597158717d0
        gaus(3,3) = 0.4701420641d0

        gaus(1,4) = 0.4701420641d0
        gaus(2,4) = 0.4701420641d0
        gaus(3,4) = 0.0597158717d0

        gaus(1,5) = 0.7974269853d0
        gaus(2,5) = 0.1012865073d0
        gaus(3,5) = 0.1012865073d0

        gaus(1,6) = 0.1012865073d0
        gaus(2,6) = 0.7974269853d0
        gaus(3,6) = 0.1012865073d0

        gaus(1,7) = 0.1012865073d0
        gaus(2,7) = 0.1012865073d0
        gaus(3,7) = 0.7974269853d0

        weis(1) = 0.2250000000d0
        weis(2) = 0.1323941527d0
        weis(3) = 0.1323941527d0
        weis(4) = 0.1323941527d0
        weis(5) = 0.1259391805d0
        weis(6) = 0.1259391805d0
	weis(7) = 0.1259391805d0
!
        return
        end

***************************************************************************
c       cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       quad3d8n.fcm                                                         c
c       ---------------------------------------------------------------------c
c     quadrature definitions in three dimensional quad reference domain    c
c       ---------------------------------------------------------------------c
c       920422 - written                                                     c
c       930311 - removed initialization of xquad                             c
c       cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine quad3d8n(iquad, nquad, xquad, wquad, maxnsd, maxnquad)

        real* 8 xquad(maxnsd,maxnquad)
        real* 8 wquad(maxnquad)

        if ((iquad.eq.1).or.(iquad.eq.131313)) then
c               1 x 1 x 1 Gaussian quadrature - code 1.1,1.1,1.1 - old code 1
                nquad = 1

                xquad(1,1) = + 0.000000000000000
                xquad(2,1) = + 0.000000000000000
                xquad(3,1) = + 0.000000000000000

                wquad(1) = 8.0
                return

        else if ((iquad.eq.2).or.(iquad.eq.121212)) then
c               2 x 2 x 2 Gaussian quadrature - code 1.2,1.2,1.2 - old code 2
                nquad = 8

                xquad(1,1) = - 0.577350269189626
                xquad(2,1) = - 0.577350269189626
                xquad(3,1) = - 0.577350269189626
                xquad(1,2) = + 0.577350269189626
                xquad(2,2) = - 0.577350269189626
                xquad(3,2) = - 0.577350269189626
                xquad(1,3) = - 0.577350269189626
                xquad(2,3) = + 0.577350269189626
                xquad(3,3) = - 0.577350269189626
                xquad(1,4) = + 0.577350269189626
                xquad(2,4) = + 0.577350269189626
                xquad(3,4) = - 0.577350269189626
                xquad(1,5) = - 0.577350269189626
                xquad(2,5) = - 0.577350269189626
                xquad(3,5) = + 0.577350269189626
                xquad(1,6) = + 0.577350269189626
                xquad(2,6) = - 0.577350269189626
                xquad(3,6) = + 0.577350269189626
                xquad(1,7) = - 0.577350269189626
                xquad(2,7) = + 0.577350269189626
                xquad(3,7) = + 0.577350269189626
                xquad(1,8) = + 0.577350269189626
                xquad(2,8) = + 0.577350269189626
                xquad(3,8) = + 0.577350269189626

                wquad(1) = 1.0
                wquad(2) = 1.0
                wquad(3) = 1.0
                wquad(4) = 1.0
                wquad(5) = 1.0
                wquad(6) = 1.0
                wquad(7) = 1.0
                wquad(8) = 1.0
                return

        else if (iquad.eq.3) then
c               3 x 3 x 3 Gaussian quadrature - code 1.3,1.3,1.3 - old code 3
                nquad = 6
                xquad(1,1) = 1.0
                xquad(2,1) = 0.0
                xquad(3,1) = 0.0
                xquad(1,2) = 0.0
                xquad(2,2) = 1.0
                xquad(3,2) = 0.0
                xquad(1,3) = 0.0
                xquad(2,3) = 0.0
                xquad(3,3) = 1.0
                xquad(1,4) =-1.0
                xquad(2,4) = 0.0
                xquad(3,4) = 0.0
                xquad(1,5) = 0.0
                xquad(2,5) =-1.0
                xquad(3,5) = 0.0
                xquad(1,6) = 0.0
                xquad(2,6) = 0.0
                xquad(3,6) =-1.0
                wquad(1) = 1.6666667
                wquad(2) = 1.6666667
                wquad(3) = 1.6666667
                wquad(4) = 1.6666667
                wquad(5) = 1.6666667
                wquad(6) = 1.6666667

                return

        else if ((iquad.eq.4).or.(iquad.eq.222222)) then
c       2 x 2 x 2 Lobatto quadrature - element corners - code 2.2,2.2,2.2 - old 4
                nquad = 8

                xquad(1,1) = - 1.0
                xquad(2,1) = - 1.0
                xquad(3,1) = - 1.0
                xquad(1,2) = + 1.0
                xquad(2,2) = - 1.0
                xquad(3,2) = - 1.0
                xquad(1,3) = - 1.0
                xquad(2,3) = + 1.0
                xquad(3,3) = - 1.0
                xquad(1,4) = + 1.0
                xquad(2,4) = + 1.0
                xquad(3,4) = - 1.0
                xquad(1,5) = - 1.0
                xquad(2,5) = - 1.0
                xquad(3,5) = + 1.0
                xquad(1,6) = + 1.0
                xquad(2,6) = - 1.0
                xquad(3,6) = + 1.0
                xquad(1,7) = - 1.0
                xquad(2,7) = + 1.0
                xquad(3,7) = + 1.0
                xquad(1,8) = + 1.0
                xquad(2,8) = + 1.0
                xquad(3,8) = + 1.0
        
                wquad(1) = 1.0
                wquad(2) = 1.0
                wquad(3) = 1.0
                wquad(4) = 1.0
                wquad(5) = 1.0
                wquad(6) = 1.0
                wquad(7) = 1.0
                wquad(8) = 1.0
        
                return

        else if (iquad.eq.221212) then
c               2 x 2 x 2 mixed quadrature - code 2.2,1.2,1.2
                nquad = 8

                xquad(1,1) = - 0.577350269189626
                xquad(2,1) = - 0.577350269189626
                xquad(3,1) = - 1.000000000000000
                xquad(1,2) = + 0.577350269189626
                xquad(2,2) = - 0.577350269189626
                xquad(3,2) = - 1.000000000000000
                xquad(1,3) = - 0.577350269189626
                xquad(2,3) = + 0.577350269189626
                xquad(3,3) = - 1.000000000000000
                xquad(1,4) = + 0.577350269189626
                xquad(2,4) = + 0.577350269189626
                xquad(3,4) = - 1.000000000000000
                xquad(1,5) = - 0.577350269189626
                xquad(2,5) = - 0.577350269189626
                xquad(3,5) = + 1.000000000000000
                xquad(1,6) = + 0.577350269189626
                xquad(2,6) = - 0.577350269189626
                xquad(3,6) = + 1.000000000000000
                xquad(1,7) = - 0.577350269189626
                xquad(2,7) = + 0.577350269189626
                xquad(3,7) = + 1.000000000000000
                xquad(1,8) = + 0.577350269189626
                xquad(2,8) = + 0.577350269189626
                xquad(3,8) = + 1.000000000000000

                wquad(1) = 1.0
                wquad(2) = 1.0
                wquad(3) = 1.0
                wquad(4) = 1.0
                wquad(5) = 1.0
                wquad(6) = 1.0
                wquad(7) = 1.0
                wquad(8) = 1.0
                return

        else
                call error("quad3d8n: unknown quadrature code ", iquad, .true.)
        end if

        end

        subroutine quad2d4n(igaus, ngaus, gaus, weis,maxnsd,nen)
        implicit none

c*****************************************************
c       INPUT:  igaus = 2
c       ngaus = 2x2 = 4
c       OUTPUT:
c       gaus(maxnsd,ngauss) = gauss pt. "natural" coors
c       weis(ngauss)     = gauss pt. weight 
c*****************************************************
        integer  nsd
        parameter( nsd = 2)
c       real*8 gaus(maxnsd,1), weis(1)
        integer  ngaus,igaus,maxnsd,nen
        real*8 gaus(maxnsd,9), weis(9)

        goto (1,2,3,4), igaus

 1    continue
        ngaus = 1

        gaus(1,1) = + 0.000000000000000
        gaus(2,1) = + 0.000000000000000

        weis(1) = 4.0
        return

 2    continue
        ngaus = 4

        gaus(1,1) = - 0.577350269189626
        gaus(2,1) = - 0.577350269189626
        gaus(1,2) = + 0.577350269189626
        gaus(2,2) = - 0.577350269189626
        gaus(1,3) = - 0.577350269189626
        gaus(2,3) = + 0.577350269189626
        gaus(1,4) = + 0.577350269189626
        gaus(2,4) = + 0.577350269189626

        weis(1) = 1.0
        weis(2) = 1.0
        weis(3) = 1.0
        weis(4) = 1.0
        return

 3    continue
        ngaus = 9

        gaus(1,1) = - 0.774596669241483
        gaus(2,1) = - 0.774596669241483
        gaus(1,2) = + 0.000000000000000
        gaus(2,2) = - 0.774596669241483
        gaus(1,3) = + 0.774596669241483
        gaus(2,3) = - 0.774596669241483
        gaus(1,4) = - 0.774596669241483
        gaus(2,4) = + 0.000000000000000
        gaus(1,5) = + 0.000000000000000
        gaus(2,5) = + 0.000000000000000
        gaus(1,6) = + 0.774596669241483
        gaus(2,6) = + 0.000000000000000
        gaus(1,7) = - 0.774596669241483
        gaus(2,7) = + 0.774596669241483
        gaus(1,8) = + 0.000000000000000
        gaus(2,8) = + 0.774596669241483
        gaus(1,9) = + 0.774596669241483
        gaus(2,9) = + 0.774596669241483

        weis(1) = 0.3086419753
        weis(2) = 0.4938271605
        weis(3) = 0.3086419753
        weis(4) = 0.4938271605
        weis(5) = 0.7901234568
        weis(6) = 0.4938271605
        weis(7) = 0.3086419753
        weis(8) = 0.4938271605
        weis(9) = 0.3086419753
        return

 4    continue
        ngaus = 4

        gaus(1,1) = - 1.0
        gaus(2,1) = - 1.0
        gaus(1,2) = + 1.0
        gaus(2,2) = - 1.0
        gaus(1,3) = - 1.0
        gaus(2,3) = + 1.0
        gaus(1,4) = + 1.0
        gaus(2,4) = + 1.0

        weis(1) = 1.0
        weis(2) = 1.0
        weis(3) = 1.0
        weis(4) = 1.0

        return
        end

c---------------------------------------------------------------------c
c quadrature definitions in two dimensional traingle reference domain c
c---------------------------------------------------------------------c

        subroutine quad2d3n(iquad, nquad, xquad, wquad, maxnsd, maxnquad)

        real*8  xquad(maxnsd,maxnquad) 
        real*8  wquad(maxnquad) 

        goto(10,20,30,40),iquad

c       1 point quadrature - code 3.1 - old code 1
 10             nquad = 1

                xquad(1,1) = + 0.333333333333333
                xquad(2,1) = + 0.333333333333333

                wquad(1) = 0.5
                return

c       3 point quadrature - code 3.2 - old code 2 - precision 2
 20             nquad = 3

                xquad(1,1) = + 0.166666666666667
                xquad(2,1) = + 0.166666666666667
                xquad(1,2) = + 0.166666666666667
                xquad(2,2) = + 0.666666666666667
                xquad(1,3) = + 0.666666666666667
                xquad(2,3) = + 0.166666666666667

                wquad(1) = 0.166666666666667
                wquad(2) = 0.166666666666667
                wquad(3) = 0.166666666666667
                return

c       4 point quadrature - code 3.4 - old code 3 - precision 3
 30             nquad = 4

                xquad(1,1) = + 0.333333333333333
                xquad(2,1) = + 0.333333333333333
                xquad(1,2) = + 0.600000000000000
                xquad(2,2) = + 0.200000000000000
                xquad(1,3) = + 0.200000000000000
                xquad(2,3) = + 0.600000000000000
                xquad(1,4) = + 0.200000000000000
                xquad(2,4) = + 0.200000000000000

                wquad(1) = - 0.281250000000000
                wquad(2) = + 0.260416666666667
                wquad(3) = + 0.260416666666667
                wquad(4) = + 0.260416666666667
                return

c       3 point Lobatto quadrature - element corners - code 2.2,2.2 - old code 4
 40             nquad = 3

                xquad(1,1) = + 0.0
                xquad(2,1) = + 0.0
                xquad(1,2) = + 1.0
                xquad(2,2) = + 0.0
                xquad(1,3) = + 0.0
                xquad(2,3) = + 1.0
        
                wquad(1) = 0.166666666666667
                wquad(2) = 0.166666666666667
                wquad(3) = 0.166666666666667
                return


        end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine quad4d8n(iquad, nquad, xquad, wquad, nsd, nquadpad)

        real* 8 xquad(nsd,nquadpad)
        real* 8 wquad(nquadpad)

        do ig=1,nquadpad
        do isd=1,nsd
        xquad(isd,ig) = 0.0
        end do
        end do
cccccccccccccccccccccccccccccccccccccccccccccc
            if (iquad.eq.1) then
                nquad = 1
                wquad(1) = 16.0
                return
cccccccccccccccccccccccccccccccccccccccccccccc
       else if (iquad.eq.3.or.iquad.eq.131313) then
                nquad = 8

                xquad(1,1) = - 0.577350269189626
                xquad(2,1) = - 0.577350269189626
                xquad(3,1) = - 0.577350269189626
                xquad(4,1) =   0.000000000000000

                xquad(1,2) = + 0.577350269189626
                xquad(2,2) = - 0.577350269189626
                xquad(3,2) = - 0.577350269189626
                xquad(4,2) =   0.000000000000000

                xquad(1,3) = - 0.577350269189626
                xquad(2,3) = + 0.577350269189626
                xquad(3,3) = - 0.577350269189626
                xquad(4,3) =   0.000000000000000

                xquad(1,4) = + 0.577350269189626
                xquad(2,4) = + 0.577350269189626
                xquad(3,4) = - 0.577350269189626
                xquad(4,4) =   0.000000000000000

                xquad(1,5) = - 0.577350269189626
                xquad(2,5) = - 0.577350269189626
                xquad(3,5) = + 0.577350269189626
                xquad(4,5) =   0.000000000000000

                xquad(1,6) = + 0.577350269189626
                xquad(2,6) = - 0.577350269189626
                xquad(3,6) = + 0.577350269189626
                xquad(4,6) =   0.000000000000000

                xquad(1,7) = - 0.577350269189626
                xquad(2,7) = + 0.577350269189626
                xquad(3,7) = + 0.577350269189626
                xquad(4,7) =   0.000000000000000

                xquad(1,8) = + 0.577350269189626
                xquad(2,8) = + 0.577350269189626
                xquad(3,8) = + 0.577350269189626
                xquad(4,8) =   0.000000000000000

                wquad(1)   = 2.0                  
                wquad(2)   = 2.0                 
                wquad(3)   = 2.0                 
                wquad(4)   = 2.0                 
                wquad(5)   = 2.0                 
                wquad(6)   = 2.0                 
                wquad(7)   = 2.0                 
                wquad(8)   = 2.0                 

ccccccccccccccccccccccccccccccccccccccccccccccc
        else if (iquad.eq.2) then
                nquad = 16
                xquad(1,1) = - 0.577350269189626
                xquad(2,1) = - 0.577350269189626
                xquad(3,1) = - 0.577350269189626
                xquad(4,1) = - 0.577350269189626

                xquad(1,2) = + 0.577350269189626
                xquad(2,2) = - 0.577350269189626
                xquad(3,2) = - 0.577350269189626
                xquad(4,2) = - 0.577350269189626

                xquad(1,3) = - 0.577350269189626
                xquad(2,3) = + 0.577350269189626
                xquad(3,3) = - 0.577350269189626
                xquad(4,3) = - 0.577350269189626

                xquad(1,4) = + 0.577350269189626
                xquad(2,4) = + 0.577350269189626
                xquad(3,4) = - 0.577350269189626
                xquad(4,4) = - 0.577350269189626

                xquad(1,5) = - 0.577350269189626
                xquad(2,5) = - 0.577350269189626
                xquad(3,5) = + 0.577350269189626
                xquad(4,5) = - 0.577350269189626

                xquad(1,6) = + 0.577350269189626
                xquad(2,6) = - 0.577350269189626
                xquad(3,6) = + 0.577350269189626
                xquad(4,6) = - 0.577350269189626

                xquad(1,7) = - 0.577350269189626
                xquad(2,7) = + 0.577350269189626
                xquad(3,7) = + 0.577350269189626
                xquad(4,7) = - 0.577350269189626

                xquad(1,8) = + 0.577350269189626
                xquad(2,8) = + 0.577350269189626
                xquad(3,8) = + 0.577350269189626
                xquad(4,8) = - 0.577350269189626

                xquad(1,9) = - 0.577350269189626
                xquad(2,9) = - 0.577350269189626
                xquad(3,9) = - 0.577350269189626
                xquad(4,9) = + 0.577350269189626

                xquad(1,10) = + 0.577350269189626
                xquad(2,10) = - 0.577350269189626
                xquad(3,10) = - 0.577350269189626
                xquad(4,10) = + 0.577350269189626

                xquad(1,11) = - 0.577350269189626
                xquad(2,11) = + 0.577350269189626
                xquad(3,11) = - 0.577350269189626
                xquad(4,11) = + 0.577350269189626

                xquad(1,12) = + 0.577350269189626
                xquad(2,12) = + 0.577350269189626
                xquad(3,12) = - 0.577350269189626
                xquad(4,12) = + 0.577350269189626

                xquad(1,13) = - 0.577350269189626
                xquad(2,13) = - 0.577350269189626
                xquad(3,13) = + 0.577350269189626
                xquad(4,13) = + 0.577350269189626

                xquad(1,14) = + 0.577350269189626
                xquad(2,14) = - 0.577350269189626
                xquad(3,14) = + 0.577350269189626
                xquad(4,14) = + 0.577350269189626

                xquad(1,15) = - 0.577350269189626
                xquad(2,15) = + 0.577350269189626
                xquad(3,15) = + 0.577350269189626
                xquad(4,15) = + 0.577350269189626

                xquad(1,16) = + 0.577350269189626
                xquad(2,16) = + 0.577350269189626
                xquad(3,16) = + 0.577350269189626
                xquad(4,16) = + 0.577350269189626

                wquad(1) = 1.0
                wquad(2) = 1.0
                wquad(3) = 1.0
                wquad(4) = 1.0
                wquad(5) = 1.0
                wquad(6) = 1.0
                wquad(7) = 1.0
                wquad(8) = 1.0
                wquad(9) = 1.0
                wquad(10) = 1.0
                wquad(11) = 1.0
                wquad(12) = 1.0
                wquad(13) = 1.0
                wquad(14) = 1.0
                wquad(15) = 1.0
                wquad(16) = 1.0
 
                return
                endif

                return
                end

        subroutine quad3d6n(iquad, nquad, xquad, wquad, nsd, nquadpad)
        real* 8 xquad(nsd,nquadpad)
        real* 8 wquad(nquadpad)

        if (iquad.eq.1) then

                nquad = 2

                xquad(1,1) = + 1.0/3.0          
                xquad(2,1) = + 1.0/3.0            
                xquad(3,1) = - 0.577350269189626

                xquad(1,1) = + 1.0/3.0          
                xquad(2,1) = + 1.0/3.0            
                xquad(3,1) = + 0.577350269189626

                wquad(1) = 0.5 
                wquad(2) = 0.5 

        return

        else if ((iquad.eq.2)) then
                nquad = 6

                xquad(1,1) = + 0.166666666666667
                xquad(2,1) = + 0.166666666666667
                xquad(3,1) = - 0.577350269189626
                xquad(1,2) = + 0.166666666666667
                xquad(2,2) = + 0.666666666666667
                xquad(3,2) = - 0.577350269189626
                xquad(1,3) = + 0.666666666666667
                xquad(2,3) = + 0.166666666666667
                xquad(3,3) = - 0.577350269189626

                xquad(1,4) = + 0.166666666666667
                xquad(2,4) = + 0.166666666666667
                xquad(3,4) = + 0.577350269189626
                xquad(1,5) = + 0.166666666666667
                xquad(2,5) = + 0.666666666666667
                xquad(3,5) = + 0.577350269189626
                xquad(1,6) = + 0.666666666666667
                xquad(2,6) = + 0.166666666666667
                xquad(3,6) = + 0.577350269189626

                wquad(1) = 0.166666666666667
                wquad(2) = 0.166666666666667
                wquad(3) = 0.166666666666667
                wquad(4) = 0.166666666666667
                wquad(5) = 0.166666666666667
                wquad(6) = 0.166666666666667

        return

        else
                call error("quad3d6n: unknown quadrature code ", iquad, .true.)
        end if

        end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine quad4d4n(iquad, nquad, xquad, wquad, nsd, nquadpad)

        real* 8 xquad(nsd,nquadpad)
        real* 8 wquad(nquadpad)

        do ig=1,nquadpad
        do isd=1,nsd
        xquad(isd,ig) = 0.0
        end do
        end do
cccccccccccccccccccccccccccccccccccccccccccccc
            if (iquad.eq.1) then
                nquad = 1
                wquad(1) = 1. / 3.
                return
cccccccccccccccccccccccccccccccccccccccccccccc
       else if (iquad.eq.2.or.iquad.eq.323232) then
                nquad = 4

                xquad(1,1) = 0.585410200000000
                xquad(2,1) = 0.138196600000000
                xquad(3,1) = 0.138196600000000
                xquad(4,1) = 0.000000000000000

                xquad(1,2) = 0.138196600000000
                xquad(2,2) = 0.138196600000000
                xquad(3,2) = 0.138196600000000
                xquad(4,2) = 0.000000000000000

                xquad(1,3) = 0.138196600000000
                xquad(2,3) = 0.138196600000000
                xquad(3,3) = 0.585410200000000
                xquad(4,3) = 0.000000000000000

                xquad(1,4) = 0.138196600000000
                xquad(2,4) = 0.585410200000000
                xquad(3,4) = 0.138196600000000
                xquad(4,4) = 0.000000000000000

                wquad(1)   = 1. / 12.
                wquad(2)   = 1. / 12.
                wquad(3)   = 1. / 12.
                wquad(4)   = 1. / 12.
ccccccccccccccccccccccccccccccccccccccccccccccc
                endif

                return
                end
