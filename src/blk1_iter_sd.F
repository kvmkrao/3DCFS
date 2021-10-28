        subroutine blk_iter_sd(ien,x1,d1,do,bp1,
     $          angle,rng,wg1,d1_t,do_t)

        implicit none
        include 'global.h'
	include 'turbpara.h'

!.......External variables.
        integer ien(nen,nec)
        integer rng(neface,nec)
        real*8  selm(nee,nee)
	integer mgau2 
	parameter(mgau2=21)
!        real*8  b(ndf,nnl)
        real*8  bp1(ndf,nnl)
        real*8  wg1(ndf,nnl)
        real*8  x1(nsd,nnl)
        real*8  d1(ndf,nnl)
        real*8  do(ndf,nnl)
        real*8  angle(nnc)
        real*8  un(ndf), unp1(ndf)
        real*8  uc(ndf), ux(ndf), uy(ndf), uz(ndf)
        real*8  uxnp1(ndf), uynp1(ndf), uznp1(ndf)
        real*8  uxn(ndf), uyn(ndf), uzn(ndf)
        real*8  ut(ndf), up(ndf),uch(ndf)


!.......Internal variables.
        real*8  ident(ndf,ndf)
        real*8  xl1(nsd,nen)
        real*8  dlo(ndf,nen),dl1(ndf,nen)
        real*8  a1(ndf,ndf), a2(ndf,ndf), a3(ndf,ndf)
        real*8  a1a1(ndf,ndf), a1a2(ndf,ndf), a1a3(ndf,ndf)
        real*8  a2a1(ndf,ndf), a2a2(ndf,ndf), a2a3(ndf,ndf)
        real*8  a3a1(ndf,ndf), a3a2(ndf,ndf), a3a3(ndf,ndf)
        real*8  a1m(ndf,ndf), a2m(ndf,ndf), a3m(ndf,ndf)
        real*8  a0i(ndf,ndf), res(ndf), dv(nsd,ndf)
        real*8  k11(ndf,ndf), k12(ndf,ndf), k13(ndf,ndf)
        real*8  k21(ndf,ndf), k22(ndf,ndf), k23(ndf,ndf)
        real*8  k31(ndf,ndf), k32(ndf,ndf), k33(ndf,ndf)
        real*8  sh(0:3,nen),xs(nsd,nsd),det
	real*8  gaus(nsd,mgau2)
	real*8  weis(mgau2)
        real*8  s0,s1,s2,s3,b0,b1,b2,b3
        real*8  s(nee,nee),p(nee)
        real*8  rotate(nee,nee),p_rot(nee)

!.......Counter variables.
        integer idof,jdof,iee,jee
        integer i,j,na,nb,nai,nbj,irng
        integer ie,inen,inode,jnode,isd,ig,iface
        integer ia,ia1,ia2,ia3,ia4,ia5
        integer na1,na2,na3,na4,na5
        integer ib
        integer ip,iq

!.......Temporary variables
        real*8  tmp,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8
!	murali
	real*8  temp1,temp2,temp3,temp4,temp5,temp6,dis1,dis2,dis3
	real*8 he_max,he_min,he_ugn,he_rgn
	real*8 umag0,umagx,umagy,umagz,grad_umag
	real*8 rx,ry,rz,sx,sy,sz
!	murali
        real*8  cf11,cf21,cf31,cf41,cf12,cf22,cf32,cf42
        real*8  cf13,cf23,cf33,cf43,cf14,cf24,cf34,cf44
        real*8  xr11,xr21,xr31,xr41,xr12,xr22,xr32,xr42
        real*8  xr13,xr23,xr33,xr43,xr14,xr24,xr34,xr44
        real*8  g,gm1,gm3,gm1b2,tmgb2
        real*8  gm1u,gm1v,gm1w,gm1uv,gm1uw,gm1vw
        real*8  uu,vv,ww,uv,uw,vw
        real*8  u,v,w,e,ge
        real*8  roi,roii,entropy,temp
        real*8  temp_down,prs_down,temp_surf
        real*8  v1,v2,v3,v4,v5
        real*8  xk1,xd1,xd2,xd3,e1,e2,e3
        real*8  rnum,rden,torr
        real*8  effs,eff0,effa,effsa,effdt
        real*8  dif,dif_mol
        real*8  gm,ro,et,xi
        real*8  gas_const,cv,t,ts,xmu
        real*8  xla,r,xkcv,prandtl,prandtl_turb
        real*8  xvv,umag,pres,c_sound,umagc
        real*8  h1,h2,hel
        real*8  tau(ndf),tau1,tau2,tau_delta,delt
        real*8  cos_a, sin_a, angle_a

!.......Matrices and Variables Related to Turbulence
	real*8 d1_t(ndf_t,nnl),do_t(ndf_t,nnl)
	real*8 dl1_t(ndf_t,nen),dlo_t(ndf_t,nen)
	real*8 nut_n,nut_np1,nut
	real*8 ki,fv1,kin_eddy,dyn_eddy
!.......Turbulence Variables END

        integer countf(22)

        logical flag
        integer icode

        parameter(prandtl_turb=0.90d0)
        parameter(torr=0.0000000001d0)

!.......Thermodynamic variables
        prandtl = 0.72d0

        dif = viscosity       ! free stream dynamic viscosity (mu)
!        dif_mol = dif
        gas_const=1.0/gamma
        cv = gas_const/(gamma-1.0)

!.......Define an identity matrix
        do idof =1,ndf
            do jdof =1,ndf
                ident(idof,jdof) = 0.0d0
            enddo
            ident(idof,idof) = 1.0d0
        enddo

        call fclear(bp1,ndf*nnl)
        call fclear(wg1,ndf*nnl)

        call lgau3t(iquad,nquad3,gaus,weis)   !Quadrature Points

!.......ELEMENT LOOP BEGINS.............................!
        do ie=1,nec

!.......LOCALIZE DATA...................................!
            do inen=1,nen
                inode=ien(inen,ie)
!.......COORDINATES
                do isd=1,nsd
                    xl1(isd,inen)=x1(isd,inode)
                enddo
!.......DOFs
                do idof=1,ndf
                    dl1(idof,inen)=d1(idof,inode)
                    dlo(idof,inen)=do(idof,inode)
                enddo
	       if(turbulent)then
	      do idof =1,ndf_t
	        dlo_t(idof,inen)=do_t(idof,inode)
	        dl1_t(idof,inen)=d1_t(idof,inode)
	      end do  
	    endif

            enddo


!.......Initialize to zero for each element
            do iee=1,nee
                p(iee) = 0.0d0
                p_rot(iee) = 0.0d0
                do jee=1,nee
                    s(iee,jee) = 0.0d0
                    selm(iee,jee) = 0.0d0
                    rotate(iee,jee) = 0.0d0
                end do
            end do

!...........QUADRATURE LOOP BEGINS.........................!
            do ig=1,nquad3
                call sf36o(gaus(1,ig),xl1,det,sh,1,xs)
                det=dabs(det)

!...........Conservation variables at integration point
                do idof=1,ndf
                    unp1(idof) = 
     &                  sh(0,1)*dl1(idof,1) + sh(0,2)*dl1(idof,2)
     &                + sh(0,3)*dl1(idof,3) + sh(0,4)*dl1(idof,4)
     &                + sh(0,5)*dl1(idof,5) + sh(0,6)*dl1(idof,6)

                    un(idof)   = 
     &                  sh(0,1)*dlo(idof,1) + sh(0,2)*dlo(idof,2)
     &                + sh(0,3)*dlo(idof,3) + sh(0,4)*dlo(idof,4)
     &                + sh(0,5)*dlo(idof,5) + sh(0,6)*dlo(idof,6)

                    uc(idof) = alpha * unp1(idof) + oma* un(idof)

                    uxnp1(idof) = 
     &                  sh(1,1)*dl1(idof,1) + sh(1,2)*dl1(idof,2)
     &                + sh(1,3)*dl1(idof,3) + sh(1,4)*dl1(idof,4)
     &                + sh(1,5)*dl1(idof,5) + sh(1,6)*dl1(idof,6)

                    uynp1(idof) = 
     &                  sh(2,1)*dl1(idof,1) + sh(2,2)*dl1(idof,2)
     &                + sh(2,3)*dl1(idof,3) + sh(2,4)*dl1(idof,4)
     &                + sh(2,5)*dl1(idof,5) + sh(2,6)*dl1(idof,6)

                    uznp1(idof) = 
     &                  sh(3,1)*dl1(idof,1) + sh(3,2)*dl1(idof,2)
     &                + sh(3,3)*dl1(idof,3) + sh(3,4)*dl1(idof,4)
     &                + sh(3,5)*dl1(idof,5) + sh(3,6)*dl1(idof,6)

                    uxn(idof)   = 
     &                  sh(1,1)*dlo(idof,1) + sh(1,2)*dlo(idof,2)
     &                + sh(1,3)*dlo(idof,3) + sh(1,4)*dlo(idof,4)
     &                + sh(1,5)*dlo(idof,5) + sh(1,6)*dlo(idof,6)

                    uyn(idof)   = 
     &                  sh(2,1)*dlo(idof,1) + sh(2,2)*dlo(idof,2)
     &                + sh(2,3)*dlo(idof,3) + sh(2,4)*dlo(idof,4)
     &                + sh(2,5)*dlo(idof,5) + sh(2,6)*dlo(idof,6)

                    uzn(idof)   = 
     &                  sh(3,1)*dlo(idof,1) + sh(3,2)*dlo(idof,2)
     &                + sh(3,3)*dlo(idof,3) + sh(3,4)*dlo(idof,4)
     &                + sh(3,5)*dlo(idof,5) + sh(3,6)*dlo(idof,6)

                    ux(idof) = alpha*uxnp1(idof) + oma*uxn(idof)
                    uy(idof) = alpha*uynp1(idof) + oma*uyn(idof)
                    uz(idof) = alpha*uznp1(idof) + oma*uzn(idof)
                    if (steady) then
                        ut(idof)= 0.0d0
                    else
                        ut(idof)= (unp1(idof)-un(idof))/dt
                    endif
                end do
	  if(turbulent)then
	    nut_np1= sh(0,1)* dl1_t(1,1) + sh(0,2)* dl1_t(1,2)
     &             + sh(0,3)* dl1_t(1,3) + sh(0,4)* dl1_t(1,4)
     &             + sh(0,5)* dl1_t(1,5) + sh(0,6)* dl1_t(1,6)

	    nut_n  = sh(0,1)* dlo_t(1,1) + sh(0,2)* dlo_t(1,2)
     &             + sh(0,3)* dlo_t(1,3) + sh(0,4)* dlo_t(1,4)
     &             + sh(0,5)* dlo_t(1,5) + sh(0,6)* dlo_t(1,6)

	    nut    = alpha * nut_np1 + oma * nut_n 
	  endif

!.......Primitive variables at integration point
                up(1) = uc(1)              !density
                up(2) = uc(2)/uc(1)        !x-velocity
                up(3) = uc(3)/uc(1)        !y-velocity
                up(4) = uc(4)/uc(1)        !w-velocity
                up(5) = uc(5)/uc(1)        !total energy

                g    =  gamma
                gm1  =  g-1.0d0
                gm3  =  g-3.0d0
                gm1b2=  gm1/2.0d0
                tmgb2= -gm3/2.0d0
                u    = up(2)
                v    = up(3)
                w    = up(4)
                e    = up(5)
                gm1u = gm1*u
                gm1v = gm1*v
                gm1w = gm1*w
                gm1uv= gm1u*v
                gm1uw= gm1u*w
                gm1vw= gm1v*w
                ge   = g*e
                uu   = u*u
                vv   = v*v
                ww   = w*w
                uv   = u*v
                uw   = u*w
                vw   = v*w

                a1(1,1) =  0.e0
                a1(1,2) =  1.e0
                a1(1,3) =  0.e0
                a1(1,4) =  0.e0
                a1(1,5) =  0.e0

                a1(2,1) = -tmgb2*uu + gm1b2*(vv+ww)
                a1(2,2) = -gm3*u
                a1(2,3) = -gm1v
                a1(2,4) = -gm1w
                a1(2,5) =  gm1

                a1(3,1) = -uv
                a1(3,2) =  v
                a1(3,3) =  u
                a1(3,4) =  0.e0
                a1(3,5) =  0.e0

                a1(4,1) = -uw
                a1(4,2) =  w
                a1(4,3) =  0.e0
                a1(4,4) =  u
                a1(4,5) =  0.e0

                a1(5,1) = -ge*u + gm1u*(uu+vv+ww)
                a1(5,2) =  ge - gm1b2*(uu+uu+uu+vv+ww)
                a1(5,3) = -gm1uv
                a1(5,4) = -gm1uw
                a1(5,5) =  g*u

                a2(1,1) =  0.e0
                a2(1,2) =  0.e0
                a2(1,3) =  1.e0
                a2(1,4) =  0.e0
                a2(1,5) =  0.e0

                a2(2,1) = -uv
                a2(2,2) =  v
                a2(2,3) =  u
                a2(2,4) =  0.e0
                a2(2,5) =  0.e0

                a2(3,1) = -tmgb2*vv + gm1b2*(uu+ww)
                a2(3,2) = -gm1u
                a2(3,3) = -gm3*v
                a2(3,4) = -gm1w
                a2(3,5) =  gm1

                a2(4,1) = -vw
                a2(4,2) =  0.e0
                a2(4,3) =  w
                a2(4,4) =  v
                a2(4,5) =  0.e0

                a2(5,1) = -ge*v + gm1v*(uu+vv+ww)
                a2(5,2) = -gm1uv
                a2(5,3) =  ge - gm1b2*(vv+vv+vv+uu+ww)
                a2(5,4) = -gm1vw
                a2(5,5) =  g*v

                a3(1,1) =  0.e0
                a3(1,2) =  0.e0
                a3(1,3) =  0.e0
                a3(1,4) =  1.e0
                a3(1,5) =  0.e0

                a3(2,1) = -uw
                a3(2,2) =  w
                a3(2,3) =  0.e0
                a3(2,4) =  u
                a3(2,5) =  0.e0

                a3(3,1) = -vw
                a3(3,2) =  0.e0
                a3(3,3) =  w
                a3(3,4) =  v
                a3(3,5) =  0.e0

                a3(4,1) = -tmgb2*ww + gm1b2*(uu+vv)
                a3(4,2) = -gm1u
                a3(4,3) = -gm1v
                a3(4,4) = -gm3*w
                a3(4,5) =  gm1

                a3(5,1) = -ge*w + gm1w*(uu+vv+ww)
                a3(5,2) = -gm1uw
                a3(5,3) = -gm1vw
                a3(5,4) =  ge - gm1b2*(ww+ww+ww+uu+vv)
                a3(5,5) =  g*w


!.......Shock capturing coefficient
                if (shock_cap) then
                    roi = uc(5)-
     &                0.5e0*(uc(2)*uc(2)+uc(3)*uc(3)+uc(4)*uc(4))/uc(1)
                    if(roi.le.0.0)then
                        print *,"-ve roi Element No:",ie,"  igaus:",ig
                    endif  
                    entropy = log(gm1*roi/uc(1)**gamma)
                    roii = 1.0d0/roi

                    v1 =  (-uc(5) + roi*(gamma+1.0e0-entropy)) * roii
                    v2 =    uc(2)*roii
                    v3 =    uc(3)*roii
                    v4 =    uc(4)*roii
                    v5 =   -uc(1)*roii

                    xk1=(v2*v2+v3*v3+v4*v4)/2.0d0/v5

                    xd1=-v2*v3
                    xd2=-v2*v4
                    xd3=-v3*v4

                    e1 = v2*v5
                    e2 = v3*v5
                    e3 = v4*v5

                    temp = -1.0d0/roi/v5

!.......Form A0 inverse
                    a0i(1,1)=temp*(xk1*xk1+g)
                    a0i(1,2)=temp*xk1*v2
                    a0i(1,3)=temp*xk1*v3
                    a0i(1,4)=temp*xk1*v4
                    a0i(1,5)=temp*(xk1+1.0e0)*v5

                    a0i(2,1)=     a0i(1,2)
                    a0i(2,2)=temp*(v2*v2-v5)
                    a0i(2,3)=temp*(-xd1)
                    a0i(2,4)=temp*(-xd2)
                    a0i(2,5)=temp*e1

                    a0i(3,1)=     a0i(1,3)
                    a0i(3,2)=     a0i(2,3) 
                    a0i(3,3)=temp*(v3*v3-v5)
                    a0i(3,4)=temp*(-xd3)
                    a0i(3,5)=temp*e2

                    a0i(4,1)=     a0i(1,4)
                    a0i(4,2)=     a0i(2,4) 
                    a0i(4,3)=     a0i(3,4)
                    a0i(4,4)=temp*(v4*v4-v5)
                    a0i(4,5)=temp*e3

                    a0i(5,1)=     a0i(1,5)
                    a0i(5,2)=     a0i(2,5)
                    a0i(5,3)=     a0i(3,5)
                    a0i(5,4)=     a0i(4,5)
                    a0i(5,5)=temp*v5*v5


!.......Residual computation

                    res(1) =
     +                  ut(1) +
     +                  a1(1,1)*ux(1) + a1(1,2)*ux(2) +
     +                  a1(1,3)*ux(3) + a1(1,4)*ux(4) + a1(1,5)*ux(5) +
     +                  a2(1,1)*uy(1) + a2(1,2)*uy(2) +
     +                  a2(1,3)*uy(3) + a2(1,4)*uy(4) + a2(1,5)*uy(5) +
     +                  a3(1,1)*uz(1) + a3(1,2)*uz(2) +
     +                  a3(1,3)*uz(3) + a3(1,4)*uz(4) + a3(1,5)*uz(5)

                    res(2) = 
     +                  ut(2) +
     +                  a1(2,1)*ux(1) + a1(2,2)*ux(2) +
     +                  a1(2,3)*ux(3) + a1(2,4)*ux(4) + a1(2,5)*ux(5) +
     +                  a2(2,1)*uy(1) + a2(2,2)*uy(2) +
     +                  a2(2,3)*uy(3) + a2(2,4)*uy(4) + a2(2,5)*uy(5) +
     +                  a3(2,1)*uz(1) + a3(2,2)*uz(2) +
     +                  a3(2,3)*uz(3) + a3(2,4)*uz(4) + a3(2,5)*uz(5) 

                    res(3) = 
     +                  ut(3) +
     +                  a1(3,1)*ux(1) + a1(3,2)*ux(2) +
     +                  a1(3,3)*ux(3) + a1(3,4)*ux(4) + a1(3,5)*ux(5) +
     +                  a2(3,1)*uy(1) + a2(3,2)*uy(2) +
     +                  a2(3,3)*uy(3) + a2(3,4)*uy(4) + a2(3,5)*uy(5) +
     +                  a3(3,1)*uz(1) + a3(3,2)*uz(2) +
     +                  a3(3,3)*uz(3) + a3(3,4)*uz(4) + a3(3,5)*uz(5) 

                    res(4) = 
     +                  ut(4) +
     +                  a1(4,1)*ux(1) + a1(4,2)*ux(2) +
     +                  a1(4,3)*ux(3) + a1(4,4)*ux(4) + a1(4,5)*ux(5) +
     +                  a2(4,1)*uy(1) + a2(4,2)*uy(2) +
     +                  a2(4,3)*uy(3) + a2(4,4)*uy(4) + a2(4,5)*uy(5) +
     +                  a3(4,1)*uz(1) + a3(4,2)*uz(2) +
     +                  a3(4,3)*uz(3) + a3(4,4)*uz(4) + a3(4,5)*uz(5) 

                    res(5) = 
     +                  ut(5) +
     +                  a1(5,1)*ux(1) + a1(5,2)*ux(2) +
     +                  a1(5,3)*ux(3) + a1(5,4)*ux(4) + a1(5,5)*ux(5) +
     +                  a2(5,1)*uy(1) + a2(5,2)*uy(2) +
     +                  a2(5,3)*uy(3) + a2(5,4)*uy(4) + a2(5,5)*uy(5) +
     +                  a3(5,1)*uz(1) + a3(5,2)*uz(2) +
     +                  a3(5,3)*uz(3) + a3(5,4)*uz(4) + a3(5,5)*uz(5) 

                    rnum = 
     +                  res(1) * ( a0i(1,1)*res(1) + a0i(1,2)*res(2) +
     +                    a0i(1,3)*res(3) + a0i(1,4)*res(4) +
     +                    a0i(1,5)*res(5)                    )
     +                + res(2) * ( a0i(2,1)*res(1) + a0i(2,2)*res(2) +
     +                    a0i(2,3)*res(3) + a0i(2,4)*res(4) +
     +                    a0i(2,5)*res(5)                    )
     +                + res(3) * ( a0i(3,1)*res(1) + a0i(3,2)*res(2) +
     +                    a0i(3,3)*res(3) + a0i(3,4)*res(4) +
     +                    a0i(3,5)*res(5)                    )
     +                + res(4) * ( a0i(4,1)*res(1) + a0i(4,2)*res(2) +
     +                    a0i(4,3)*res(3) + a0i(4,4)*res(4) +
     +                    a0i(4,5)*res(5)                    )
     +                + res(5) * ( a0i(5,1)*res(1) + a0i(5,2)*res(2) +
     +                    a0i(5,3)*res(3) + a0i(5,4)*res(4) +
     +                    a0i(5,5)*res(5)                    )

                    dv(1,1) = 
     &                  xs(1,1)*ux(1) + xs(1,2)*uy(1) + xs(1,3)*uz(1)
                    dv(1,2) = 
     &                  xs(1,1)*ux(2) + xs(1,2)*uy(2) + xs(1,3)*uz(2)
                    dv(1,3) = 
     &                  xs(1,1)*ux(3) + xs(1,2)*uy(3) + xs(1,3)*uz(3)
                    dv(1,4) = 
     &                  xs(1,1)*ux(4) + xs(1,2)*uy(4) + xs(1,3)*uz(4)
                    dv(1,5) = 
     &                  xs(1,1)*ux(5) + xs(1,2)*uy(5) + xs(1,3)*uz(5)

                    dv(2,1) = 
     &                  xs(2,1)*ux(1) + xs(2,2)*uy(1) + xs(2,3)*uz(1)
                    dv(2,2) = 
     &                  xs(2,1)*ux(2) + xs(2,2)*uy(2) + xs(2,3)*uz(2)
                    dv(2,3) = 
     &                  xs(2,1)*ux(3) + xs(2,2)*uy(3) + xs(2,3)*uz(3)
                    dv(2,4) = 
     &                  xs(2,1)*ux(4) + xs(2,2)*uy(4) + xs(2,3)*uz(4)
                    dv(2,5) = 
     &                  xs(2,1)*ux(5) + xs(2,2)*uy(5) + xs(2,3)*uz(5)

                    dv(3,1) = 
     &                  xs(3,1)*ux(1) + xs(3,2)*uy(1) + xs(3,3)*uz(1)
                    dv(3,2) = 
     &                  xs(3,1)*ux(2) + xs(3,2)*uy(2) + xs(3,3)*uz(2)
                    dv(3,3) = 
     &                  xs(3,1)*ux(3) + xs(3,2)*uy(3) + xs(3,3)*uz(3)
                    dv(3,4) = 
     &                  xs(3,1)*ux(4) + xs(3,2)*uy(4) + xs(3,3)*uz(4)
                    dv(3,5) = 
     &                  xs(3,1)*ux(5) + xs(3,2)*uy(5) + xs(3,3)*uz(5)


                    rden = torr + 
     +                  dv(1,1) * (a0i(1,1)*dv(1,1) + a0i(1,2)*dv(1,2)
     +                          + a0i(1,3)*dv(1,3) + a0i(1,4)*dv(1,4) 
     +                          + a0i(1,5)*dv(1,5)                    )
     +                + dv(1,2) * (a0i(2,1)*dv(1,1) + a0i(2,2)*dv(1,2)
     +                          + a0i(2,3)*dv(1,3) + a0i(2,4)*dv(1,4) 
     +                          + a0i(2,5)*dv(1,5)                    )
     +                + dv(1,3) * (a0i(3,1)*dv(1,1) + a0i(3,2)*dv(1,2)
     +                          + a0i(3,3)*dv(1,3) + a0i(3,4)*dv(1,4) 
     +                          + a0i(3,5)*dv(1,5)                    )
     +                + dv(1,4) * (a0i(4,1)*dv(1,1) + a0i(4,2)*dv(1,2)
     +                          + a0i(4,3)*dv(1,3) + a0i(4,4)*dv(1,4) 
     +                          + a0i(4,5)*dv(1,5)                    )
     +                + dv(1,5) * (a0i(5,1)*dv(1,1) + a0i(5,2)*dv(1,2)
     +                          + a0i(5,3)*dv(1,3) + a0i(5,4)*dv(1,4) 
     +                          + a0i(5,5)*dv(1,5)                    )

     +                + dv(2,1) * (a0i(1,1)*dv(2,1) + a0i(1,2)*dv(2,2)
     +                          + a0i(1,3)*dv(2,3) + a0i(1,4)*dv(2,4) 
     +                          + a0i(1,5)*dv(2,5)                    )
     +                + dv(2,2) * (a0i(2,1)*dv(2,1) + a0i(2,2)*dv(2,2)
     +                          + a0i(2,3)*dv(2,3) + a0i(2,4)*dv(2,4) 
     +                          + a0i(2,5)*dv(2,5)                    )
     +                + dv(2,3) * (a0i(3,1)*dv(2,1) + a0i(3,2)*dv(2,2)
     +                          + a0i(3,3)*dv(2,3) + a0i(3,4)*dv(2,4) 
     +                          + a0i(3,5)*dv(2,5)                    )
     +                + dv(2,4) * (a0i(4,1)*dv(2,1) + a0i(4,2)*dv(2,2)
     +                          + a0i(4,3)*dv(2,3) + a0i(4,4)*dv(2,4) 
     +                          + a0i(4,5)*dv(2,5)                    )
     +                + dv(2,5) * (a0i(5,1)*dv(2,1) + a0i(5,2)*dv(2,2)
     +                          + a0i(5,3)*dv(2,3) + a0i(5,4)*dv(2,4) 
     +                          + a0i(5,5)*dv(2,5)                    )

     +                + dv(3,1) * (a0i(1,1)*dv(3,1) + a0i(1,2)*dv(3,2)
     +                          + a0i(1,3)*dv(3,3) + a0i(1,4)*dv(3,4) 
     +                          + a0i(1,5)*dv(3,5)                    )
     +                + dv(3,2) * (a0i(2,1)*dv(3,1) + a0i(2,2)*dv(3,2)
     +                          + a0i(2,3)*dv(3,3) + a0i(2,4)*dv(3,4) 
     +                          + a0i(2,5)*dv(3,5)                    )
     +                + dv(3,3) * (a0i(3,1)*dv(3,1) + a0i(3,2)*dv(3,2)
     +                          + a0i(3,3)*dv(3,3) + a0i(3,4)*dv(3,4) 
     +                          + a0i(3,5)*dv(3,5)                    )
     +                + dv(3,4) * (a0i(4,1)*dv(3,1) + a0i(4,2)*dv(3,2)
     +                          + a0i(4,3)*dv(3,3) + a0i(4,4)*dv(3,4) 
     +                          + a0i(4,5)*dv(3,5)                    )
     +                + dv(3,5) * (a0i(5,1)*dv(3,1) + a0i(5,2)*dv(3,2)
     +                          + a0i(5,3)*dv(3,3) + a0i(5,4)*dv(3,4) 
     +                          + a0i(5,5)*dv(3,5)                    )


                    delt = ShockCapFactor*sqrt(rnum/rden)  !ShockCapFactor=0.5 
                else
                    delt = 0.0d0
                end if

                eff0 = weis(ig) * det
                effs = eff0 * delt


                gm = gamma
                ro = uc(1)
                u  = uc(2)/uc(1)
                v  = uc(3)/uc(1)
                w  = uc(4)/uc(1)
                et = uc(5)/uc(1)
                vv = u*u + v*v + w*w
                xi = et - 0.5d0*vv     
                t  = xi / cv

                ts = 110.0d0   !CONSTANT in SUTHERLAND's LAW of VISCOSITY

!.......viscosity constant (independent of temp)
!                xmu= dif
!.......SutherLands Law.
!.......Ref Temp: 288.16 K (equal to freestream temp)
!.......Formulae: mu/mu_ref=((T/ref_temp)^3/2)*((ref_temp+110)/(T+110))
!.......Afer Non-Dimensionaling above, we get:
!.......dividing the num and denom. by T_ref
!                if(Sutherland_Law)then
                    ts=ts/288.16d0
                    xmu= dif * ((t)**1.5d0)*((1.0d0+ts)/(t+ts))
!                endif
c.......Exponential Law
!         xmu= dif*(t**0.768)

                dyn_eddy=0.0d0

	  if(turbulent)then
	    ki = nut/(xmu/ro)
	    fv1= (ki**3)/(ki**3+cv1**3)
	    kin_eddy = nut*fv1
	    dyn_eddy = kin_eddy*ro
	  endif

                xla =-2.0d0/3.0d0*(xmu + dyn_eddy)  
                r   = 1.0d0/ro
                xkcv= gm * (xmu/prandtl + dyn_eddy/prandtl_turb)
                xmu = xmu + dyn_eddy

                xvv = (xkcv*0.5d0 -xmu)*vv
                
                k11(1,1) = 0.0e0
                k11(1,2) = 0.0e0
                k11(1,3) = 0.0e0
                k11(1,4) = 0.0e0
                k11(1,5) = 0.0e0
                k11(2,1) =-u*4.0d0/3.0d0*xmu*r
                k11(2,2) =   4.0d0/3.0d0*xmu*r
                k11(2,3) = 0.0e0
                k11(2,4) = 0.0e0
                k11(2,5) = 0.0e0
                k11(3,1) =-v*xmu            *r
                k11(3,2) = 0.0e0
                k11(3,3) =   xmu            *r
                k11(3,4) = 0.0e0
                k11(3,5) = 0.0e0
                k11(4,1) =-w*xmu            *r
                k11(4,2) = 0.0e0
                k11(4,3) = 0.0e0
                k11(4,4) =   xmu            *r
                k11(4,5) = 0.0e0
                k11(5,1) =(-xkcv*xi + xvv - xmu/3.0d0*u*u) * r
                k11(5,2) =(4.0d0/3.0d0*xmu-xkcv)*u*r
                k11(5,3) =(xmu - xkcv)*v    *r
                k11(5,4) =(xmu - xkcv)*w    *r
                k11(5,5) = xkcv             *r
        
                k12(1,1) = 0.0e0
                k12(1,2) = 0.0e0
                k12(1,3) = 0.0e0
                k12(1,4) = 0.0e0
                k12(1,5) = 0.0e0
                k12(2,1) = v*2.0d0/3.0d0*xmu *r
                k12(2,2) = 0.0e0
                k12(2,3) =  -2.0d0/3.0d0*xmu *r
                k12(2,4) = 0.0e0
                k12(2,5) = 0.0e0
                k12(3,1) = -u*xmu            *r
                k12(3,2) =  xmu              *r
                k12(3,3) = 0.0e0
                k12(3,4) = 0.0e0
                k12(3,5) = 0.0e0
                k12(4,1) = 0.0e0
                k12(4,2) = 0.0e0
                k12(4,3) = 0.0e0
                k12(4,4) = 0.0e0
                k12(4,5) = 0.0e0
                k12(5,1) = -xmu/3.0d0 *u*v   *r
                k12(5,2) =  xmu * v          *r
                k12(5,3) = -xmu*2.0d0/3.0d0*u*r
                k12(5,4) = 0.0e0
                k12(5,5) = 0.0e0
        
                k13(1,1) = 0.0e0
                k13(1,2) = 0.0e0
                k13(1,3) = 0.0e0
                k13(1,4) = 0.0e0
                k13(1,5) = 0.0e0
                k13(2,1) = w*2.0d0/3.0d0*xmu *r
                k13(2,2) = 0.0e0
                k13(2,3) = 0.0e0
                k13(2,4) =  -2.0d0/3.0d0*xmu *r
                k13(2,5) = 0.0e0
                k13(3,1) = 0.0e0
                k13(3,2) = 0.0e0
                k13(3,3) = 0.0e0
                k13(3,4) = 0.0e0
                k13(3,5) = 0.0e0
                k13(4,1) = -u*xmu            *r
                k13(4,2) =  xmu              *r
                k13(4,3) = 0.0e0
                k13(4,4) = 0.0e0
                k13(4,5) = 0.0e0
                k13(5,1) = -xmu/3.0d0 *u*w   *r
                k13(5,2) =  xmu * w          *r
                k13(5,3) = 0.0e0
                k13(5,4) = -xmu*2.0d0/3.0d0*u*r
                k13(5,5) = 0.0e0
        
                k21(1,1) = 0.0e0
                k21(1,2) = 0.0e0
                k21(1,3) = 0.0e0
                k21(1,4) = 0.0e0
                k21(1,5) = 0.0e0
                k21(2,1) = -v*xmu            *r
                k21(2,2) = 0.0e0
                k21(2,3) =  xmu              *r
                k21(2,4) = 0.0e0
                k21(2,5) = 0.0e0
                k21(3,1) = u*2.0d0/3.0d0*xmu *r
                k21(3,2) =  -2.0d0/3.0d0*xmu *r
                k21(3,3) = 0.0e0
                k21(3,4) = 0.0e0
                k21(3,5) = 0.0e0
                k21(4,1) = 0.0e0
                k21(4,2) = 0.0e0
                k21(4,3) = 0.0e0
                k21(4,4) = 0.0e0
                k21(4,5) = 0.0e0
                k21(5,1) = -xmu/3.0d0 *u*v   *r
                k21(5,2) = -xmu*2.0d0/3.0d0*v*r
                k21(5,3) =  xmu * u          *r
                k21(5,4) = 0.0e0
                k21(5,5) = 0.0e0
        
                k22(1,1) = 0.0e0
                k22(1,2) = 0.0e0
                k22(1,3) = 0.0e0
                k22(1,4) = 0.0e0
                k22(1,5) = 0.0e0
                k22(2,1) =-u*xmu            *r
                k22(2,2) =   xmu            *r
                k22(2,3) = 0.0e0
                k22(2,4) = 0.0e0
                k22(2,5) = 0.0e0
                k22(3,1) =-v*4.0d0/3.0d0*xmu*r
                k22(3,2) = 0.0e0
                k22(3,3) =   4.0d0/3.0d0*xmu*r
                k22(3,4) = 0.0e0
                k22(3,5) = 0.0e0
                k22(4,1) =-w*xmu            *r
                k22(4,2) = 0.0e0
                k22(4,3) = 0.0e0
                k22(4,4) =   xmu            *r
                k22(4,5) = 0.0e0
                k22(5,1) =(-xkcv*xi + xvv - xmu/3.0d0*v*v) * r
                k22(5,2) =(xmu - xkcv)*u    *r
                k22(5,3) =(4.0d0/3.0d0*xmu-xkcv)*v*r
                k22(5,4) =(xmu - xkcv)*w    *r
                k22(5,5) = xkcv             *r
        
                k23(1,1) = 0.0e0
                k23(1,2) = 0.0e0
                k23(1,3) = 0.0e0
                k23(1,4) = 0.0e0
                k23(1,5) = 0.0e0
                k23(2,1) = 0.0e0
                k23(2,2) = 0.0e0
                k23(2,3) = 0.0e0
                k23(2,4) = 0.0e0
                k23(2,5) = 0.0e0
                k23(3,1) = w*2.0d0/3.0d0*xmu *r
                k23(3,2) = 0.0e0
                k23(3,3) = 0.0e0
                k23(3,4) =  -2.0d0/3.0d0*xmu *r
                k23(3,5) = 0.0e0
                k23(4,1) = -v*xmu            *r
                k23(4,2) = 0.0e0
                k23(4,3) =  xmu              *r
                k23(4,4) = 0.0e0
                k23(4,5) = 0.0e0
                k23(5,1) = -xmu/3.0d0 *v*w   *r
                k23(5,2) = 0.0e0
                k23(5,3) =  xmu * w          *r
                k23(5,4) = -xmu*2.0d0/3.0d0*v*r
                k23(5,5) = 0.0e0
        
                k31(1,1) = 0.0e0
                k31(1,2) = 0.0e0
                k31(1,3) = 0.0e0
                k31(1,4) = 0.0e0
                k31(1,5) = 0.0e0
                k31(2,1) = -w*xmu            *r
                k31(2,2) = 0.0e0
                k31(2,3) = 0.0e0
                k31(2,4) =  xmu              *r
                k31(2,5) = 0.0e0
                k31(3,1) = 0.0e0
                k31(3,2) = 0.0e0
                k31(3,3) = 0.0e0
                k31(3,4) = 0.0e0
                k31(3,5) = 0.0e0
                k31(4,1) = u*2.0d0/3.0d0*xmu *r
                k31(4,2) =  -2.0d0/3.0d0*xmu *r
                k31(4,3) = 0.0e0
                k31(4,4) = 0.0e0
                k31(4,5) = 0.0e0
                k31(5,1) = -xmu/3.0d0 *u*w   *r
                k31(5,2) = -xmu*2.0d0/3.0d0*w*r
                k31(5,3) = 0.0e0
                k31(5,4) =  xmu * u          *r
                k31(5,5) = 0.0e0
        
                k32(1,1) = 0.0e0
                k32(1,2) = 0.0e0
                k32(1,3) = 0.0e0
                k32(1,4) = 0.0e0
                k32(1,5) = 0.0e0
                k32(2,1) = 0.0e0
                k32(2,2) = 0.0e0
                k32(2,3) = 0.0e0
                k32(2,4) = 0.0e0
                k32(2,5) = 0.0e0
                k32(3,1) = -w*xmu            *r
                k32(3,2) = 0.0e0
                k32(3,3) = 0.0e0
                k32(3,4) =  xmu              *r
                k32(3,5) = 0.0e0
                k32(4,1) = v*2.0d0/3.0d0*xmu *r
                k32(4,2) = 0.0e0
                k32(4,3) =  -2.0d0/3.0d0*xmu *r
                k32(4,4) = 0.0e0
                k32(4,5) = 0.0e0
                k32(5,1) = -xmu/3.0d0 *v*w   *r
                k32(5,2) = 0.0e0
                k32(5,3) = -xmu*2.0d0/3.0d0*w*r
                k32(5,4) =  xmu * v          *r
                k32(5,5) = 0.0e0
        
                k33(1,1) = 0.0e0
                k33(1,2) = 0.0e0
                k33(1,3) = 0.0e0
                k33(1,4) = 0.0e0
                k33(1,5) = 0.0e0
                k33(2,1) =-u*xmu            *r
                k33(2,2) =   xmu            *r
                k33(2,3) = 0.0e0
                k33(2,4) = 0.0e0
                k33(2,5) = 0.0e0
                k33(3,1) =-v*xmu            *r
                k33(3,2) = 0.0e0
                k33(3,3) =   xmu            *r
                k33(3,4) = 0.0e0
                k33(3,5) = 0.0e0
                k33(4,1) =-w*4.0d0/3.0d0*xmu*r
                k33(4,2) = 0.0e0
                k33(4,3) = 0.0e0
                k33(4,4) =   4.0d0/3.0d0*xmu*r
                k33(4,5) = 0.0e0
                k33(5,1) =(-xkcv*xi + xvv - xmu/3.0d0*w*w) * r
                k33(5,2) =(xmu - xkcv)*u    *r
                k33(5,3) =(xmu - xkcv)*v    *r
                k33(5,4) =(4.0d0/3.0d0*xmu-xkcv)*w*r
                k33(5,5) = xkcv             *r

                umag = sqrt(up(2)*up(2) + up(3)*up(3) + up(4)*up(4))
                temp = (up(5)-0.5e0*umag*umag)/cv
                pres = up(1)*gas_const*temp
                c_sound = sqrt(gamma*gas_const*temp)
                umagc = umag + c_sound

!.......hel = Element Length (min side len / max side len / in the dir of vel)
!.......USE ONE OF THE OPTIONS
!	murali
		temp1=sqrt(dl1(1,1)*dl1(1,1)+dl1(2,1)*dl1(2,1)
     $			  +dl1(3,1)*dl1(3,1))
	        temp2=sqrt(dl1(1,2)*dl1(1,2)+d1(2,2)*dl1(2,2)
     $			  +dl1(3,2)*dl1(3,2))
	        temp3=sqrt(dl1(1,3)*dl1(1,3)+dl1(2,3)*dl1(2,3)
     $			  +dl1(3,3)*dl1(3,3))
         	temp4=sqrt(dl1(1,4)*dl1(1,4)+dl1(2,4)*dl1(2,4)
     $			  +dl1(3,4)*dl1(3,4))
         	temp5=sqrt(dl1(1,5)*dl1(1,5)+dl1(2,5)*dl1(2,5)
     $			  +dl1(3,5)*dl1(3,5))
         	temp6=sqrt(dl1(1,6)*dl1(1,6)+dl1(2,6)*dl1(2,6)
     $			  +dl1(3,6)*dl1(3,6))

        	umag0=sh(0,1)*temp1+sh(0,2)*temp2+sh(0,3)*temp3
     $		      	   +sh(0,4)*temp4
     $                     +sh(0,5)*temp5+sh(0,6)*temp6
        	umagx=sh(1,1)*temp1+sh(1,2)*temp2+sh(1,3)*temp3
     $			   +sh(1,4)*temp4
     $                     +sh(1,5)*temp5+sh(1,6)*temp6
        	umagy=sh(2,1)*temp1+sh(2,2)*temp2+sh(2,3)*temp3
     $			   +sh(2,4)*temp4
     $                     +sh(2,5)*temp5+sh(2,6)*temp6
        	umagz=sh(3,1)*temp1+sh(3,2)*temp2+sh(3,3)*temp3
     $			   +sh(3,4)*temp4
     $                     +sh(3,5)*temp5+sh(3,6)*temp6

		grad_umag=dsqrt(umagx*umagx+umagy*umagy+umagz*umagz)

	        temp1=xl1(1,1)-xl1(1,2)
	        temp2=xl1(2,1)-xl1(2,2)
	        temp3=xl1(3,1)-xl1(3,2)
	        dis1 = sqrt(temp1*temp1+temp2*temp2+temp3*temp3)
	        temp1=xl1(1,1)-xl1(1,3)
	        temp2=xl1(2,1)-xl1(2,3)
	        temp3=xl1(3,1)-xl1(3,3)
	        dis2 = sqrt(temp1*temp1+temp2*temp2+temp3*temp3)
	        temp1=xl1(1,2)-xl1(1,3)
	        temp2=xl1(2,2)-xl1(2,3)
	        temp3=xl1(3,2)-xl1(3,3)
	        dis3 = sqrt(temp1*temp1+temp2*temp2+temp3*temp3)

        	he_max=0.2*max(dis1,dis2,dis3)
		he_min=min(dis1,dis2,dis3)

         if (grad_umag .gt. 1.0e-9) then
             rx = umagx/grad_umag
             ry = umagy/grad_umag
             rz = umagz/grad_umag

             he_rgn = dabs(rx*sh(1,1)+ry*sh(2,1)+rz*sh(3,1))
     &              + dabs(rx*sh(1,2)+ry*sh(2,2)+rz*sh(3,2))
     &              + dabs(rx*sh(1,3)+ry*sh(2,3)+rz*sh(3,3))
     &              + dabs(rx*sh(1,4)+ry*sh(2,4)+rz*sh(3,4))
     &              + dabs(rx*sh(1,5)+ry*sh(2,5)+rz*sh(3,5))
     &              + dabs(rx*sh(1,6)+ry*sh(2,6)+rz*sh(3,6))

             he_rgn = 2.0/he_rgn
            else
             he_rgn = he_min
        end if

          if (umag .gt. 1.0d-6) then
            sx=un(2)/umag
            sy=un(3)/umag
            sz=un(4)/umag

            he_ugn = dabs(sx*sh(1,1)+sy*sh(2,1)+sz*sh(3,1))
     &             + dabs(sx*sh(1,2)+sy*sh(2,2)+sz*sh(3,2))
     &             + dabs(sx*sh(1,3)+sy*sh(2,3)+sz*sh(3,3))
     &             + dabs(sx*sh(1,4)+sy*sh(2,4)+sz*sh(3,4))
     &             + dabs(sx*sh(1,5)+sy*sh(2,5)+sz*sh(3,5))
     &             + dabs(sx*sh(1,6)+sy*sh(2,6)+sz*sh(3,6))

            he_ugn = 2.0/he_ugn
            else
            he_ugn = he_min
            end if

           he_ugn = he_min
           he_rgn = he_min


!1                if(flag_he.eq.1)then
!.......Element Length(hel) based on Direction of vel.
!1                    if(umag.gt.1e-8)then
!1                        h1 = up(2)/umag  !cos(theta), theta=angle of vel vector with hor.
!1                        h2 = up(3)/umag  !sin(theta)
!1
!1                        hel = abs(h1*sh(1,1)+h2*sh(2,1))
!1     $                     + abs(h1*sh(1,2)+h2*sh(2,2))
!1     $                     + abs(h1*sh(1,3)+h2*sh(2,3))
!1                        hel = 0.5d0/hel
!1                    else
!1                        hel=min(
!1     &              sqrt((xl1(1,2)-xl1(1,1))**2+(xl1(2,2)-xl1(2,1))**2),
!1     $              sqrt((xl1(1,3)-xl1(1,1))**2+(xl1(2,3)-xl1(2,1))**2),
!1     $              sqrt((xl1(1,3)-xl1(1,2))**2+(xl1(2,3)-xl1(2,2))**2))
!1                    endif
!1                else if(flag_he.eq.2) then
!.......MY Suggestion: Length based on AVERAGE of lengths of sides
!1                    hel=(
!!     &              sqrt((xl1(1,2)-xl1(1,1))**2+(xl1(2,2)-xl1(2,1))**2)+
!1     $              sqrt((xl1(1,3)-xl1(1,1))**2+(xl1(2,3)-xl1(2,1))**2)+
!1     $        sqrt((xl1(1,3)-xl1(1,2))**2+(xl1(2,3)-xl1(2,2))**2))/3.0
!1                else if(flag_he.eq.3) then
!.......Element Length(hel) based on Length of Edge of Triangle
!.......Based on MINIMUM length
!1                    hel=dmin1(
!1     &              sqrt((xl1(1,2)-xl1(1,1))**2+(xl1(2,2)-xl1(2,1))**2),
!1     $              sqrt((xl1(1,3)-xl1(1,1))**2+(xl1(2,3)-xl1(2,1))**2),
!1     $              sqrt((xl1(1,3)-xl1(1,2))**2+(xl1(2,3)-xl1(2,2))**2))
!1                else
!.......Based on MAXIMUM length
!.......NOTE: For Fine MESH max length may be used
!1                    hel=dmax1(
!1     &              sqrt((xl1(1,2)-xl1(1,1))**2+(xl1(2,2)-xl1(2,1))**2),
!1     $              sqrt((xl1(1,3)-xl1(1,1))**2+(xl1(2,3)-xl1(2,1))**2),
!1     $              sqrt((xl1(1,3)-xl1(1,2))**2+(xl1(2,3)-xl1(2,2))**2))
!1                endif
!1
!	murali
                if (dtintau) then
                   tau1=1.0d0/dsqrt((2.0d0/dt)**2+(2.0d0*umagc/he_ugn)**2 
     $           + (4.0d0*(dif/ro)/he_ugn/he_ugn)**2)
                else
!	murali
!1           	   tau1 = 1.0d0/dsqrt((2.0d0*umagc/hel)**2      
!1     $           + (12.0d0*(dif/ro)/hel/hel)**2)
		   tau1 = 1.0d0/dsqrt((2.0d0*umagc/he_ugn)**2
     $            + (12.0d0*dif/he_rgn/he_rgn)**2)	  
!	murali
!                    tau1 = hel/(2.0d0*umagc)        !2nd term of tau1 (above) dropped
                endif

                tau_delta = delt/umagc/umagc
!       	tau_delta = 0.0d0                !Change Suggested for Finer MESH
                tau2 = max(0.0e0, tau1-tau_delta)

                tau(1) = tau2
                tau(2) = tau2
                tau(3) = tau2
                tau(4) = tau2
                tau(5) = tau2

                do i=1,ndf 
                    do j=1,ndf 

!.......Nevermind the indentation
        a1a1(i,j) = a1(i,1)*tau(1)*a1(1,j) + a1(i,2)*tau(2)*a1(2,j)
     +            + a1(i,3)*tau(3)*a1(3,j) + a1(i,4)*tau(4)*a1(4,j)
     +            + a1(i,5)*tau(5)*a1(5,j)
        a1a2(i,j) = a1(i,1)*tau(1)*a2(1,j) + a1(i,2)*tau(2)*a2(2,j)
     +            + a1(i,3)*tau(3)*a2(3,j) + a1(i,4)*tau(4)*a2(4,j)
     +            + a1(i,5)*tau(5)*a2(5,j)
        a1a3(i,j) = a1(i,1)*tau(1)*a3(1,j) + a1(i,2)*tau(2)*a3(2,j)
     +            + a1(i,3)*tau(3)*a3(3,j) + a1(i,4)*tau(4)*a3(4,j)
     +            + a1(i,5)*tau(5)*a3(5,j)
        a2a1(i,j) = a2(i,1)*tau(1)*a1(1,j) + a2(i,2)*tau(2)*a1(2,j)
     +            + a2(i,3)*tau(3)*a1(3,j) + a2(i,4)*tau(4)*a1(4,j)
     +            + a2(i,5)*tau(5)*a1(5,j)
        a2a2(i,j) = a2(i,1)*tau(1)*a2(1,j) + a2(i,2)*tau(2)*a2(2,j)
     +            + a2(i,3)*tau(3)*a2(3,j) + a2(i,4)*tau(4)*a2(4,j)
     +            + a2(i,5)*tau(5)*a2(5,j)
        a2a3(i,j) = a2(i,1)*tau(1)*a3(1,j) + a2(i,2)*tau(2)*a3(2,j)
     +            + a2(i,3)*tau(3)*a3(3,j) + a2(i,4)*tau(4)*a3(4,j)
     +            + a2(i,5)*tau(5)*a3(5,j)
        a3a1(i,j) = a3(i,1)*tau(1)*a1(1,j) + a3(i,2)*tau(2)*a1(2,j)
     +            + a3(i,3)*tau(3)*a1(3,j) + a3(i,4)*tau(4)*a1(4,j)
     +            + a3(i,5)*tau(5)*a1(5,j)
        a3a2(i,j) = a3(i,1)*tau(1)*a2(1,j) + a3(i,2)*tau(2)*a2(2,j)
     +            + a3(i,3)*tau(3)*a2(3,j) + a3(i,4)*tau(4)*a2(4,j)
     +            + a3(i,5)*tau(5)*a2(5,j)
        a3a3(i,j) = a3(i,1)*tau(1)*a3(1,j) + a3(i,2)*tau(2)*a3(2,j)
     +            + a3(i,3)*tau(3)*a3(3,j) + a3(i,4)*tau(4)*a3(4,j)
     +            + a3(i,5)*tau(5)*a3(5,j)
        a1m (i,j) = a1(i,1)*tau(1)*ident(1,j) 
     +            + a1(i,2)*tau(2)*ident(2,j)
     +            + a1(i,3)*tau(3)*ident(3,j) 
     +            + a1(i,4)*tau(4)*ident(4,j)
     +            + a1(i,5)*tau(5)*ident(5,j)
        a2m (i,j) = a2(i,1)*tau(1)*ident(1,j) 
     +            + a2(i,2)*tau(2)*ident(2,j)
     +            + a2(i,3)*tau(3)*ident(3,j) 
     +            + a2(i,4)*tau(4)*ident(4,j)
     +            + a2(i,5)*tau(5)*ident(5,j)
        a3m (i,j) = a3(i,1)*tau(1)*ident(1,j) 
     +            + a3(i,2)*tau(2)*ident(2,j)
     +            + a3(i,3)*tau(3)*ident(3,j) 
     +            + a3(i,4)*tau(4)*ident(4,j)
     +            + a3(i,5)*tau(5)*ident(5,j)

                    end do
                end do

                effa  = eff0 * alpha
                effsa = effs * alpha
                effdt = eff0 / dt

                do 640 na=1,nen
                    s0 = sh(0,na)
                    s1 = sh(1,na)
                    s2 = sh(2,na)
                    s3 = sh(3,na)
                    do 640 i=1, ndf
                        nai = (na-1)*ndf + i
                        do 610 nb=1,nen
                            b0 = sh(0,nb)
                            b1 = sh(1,nb)
                            b2 = sh(2,nb)
                            b3 = sh(3,nb)
                            do 610 j=1, ndf
                                nbj = (nb-1)*ndf + j
                  
!...........Nevermind the indentation
                                s(nai,nbj) = s(nai,nbj)
     &   + (s0*ident(i,j)+s1*a1m(i,j)+s2*a2m(i,j)+s3*a3m(i,j))*b0*effdt	  !s3*a3m(i,j)
     &   + (s0*(a1(i,j)*b1+a2(i,j)*b2+a3(i,j)*b3)                         !a3(i,j)*b3
     &   +  s1*((k11(i,j)+a1a1(i,j))*b1
     &         +(k12(i,j)+a1a2(i,j))*b2
     &         +(k13(i,j)+a1a3(i,j))*b3)				  !(k13(i,j)+a1a3(i,j))*b3)
     &   +  s2*((k21(i,j)+a2a1(i,j))*b1
     &         +(k22(i,j)+a2a2(i,j))*b2
     &         +(k23(i,j)+a2a3(i,j))*b3)				  !(k23(i,j)+a2a3(i,j))*b3)
     &   +  s3*((k31(i,j)+a3a1(i,j))*b1					  !s3*((k31(i,j)+a3a1(i,j))*b1
     &         +(k32(i,j)+a3a2(i,j))*b2					  !(k32(i,j)+a3a2(i,j))*b2
     &         +(k33(i,j)+a3a3(i,j))*b3)) * effa			  !(k33(i,j)+a3a3(i,j))*b3)) * effa
     &   +  ident(i,j)*(s1*b1+s2*b2+s3*b3)*effsa			  !(s3*b3)*effsa

610                     continue

                        do 640 j=1,ndf
                            p(nai) = p(nai)
                ! first term
     &          -(s0*ident(i,j))*ut(j)*eff0
                ! second term
     &          -(s1*a1m(i,j) + s2*a2m(i,j) + s3*a3m(i,j))*ut(j)*eff0 	   !-s3*a3m(i,j))*ut(j)*eff0
                ! third term
     &          -(s0*(a1(i,j)*ux(j)+a2(i,j)*uy(j)+a3(i,j)*uz(j)))*eff0     !a3(i,j)*uz(j)*eff0
                ! fourth term
     &          -( s1*((k11(i,j)+a1a1(i,j))*ux(j)
     &                +(k12(i,j)+a1a2(i,j))*uy(j)
     &                +(k13(i,j)+a1a3(i,j))*uz(j))			   !-k13(i,j)+a1a3(i,j))*uz(j)*s1*eff0
     &            +s2*((k21(i,j)+a2a1(i,j))*ux(j)
     &                +(k22(i,j)+a2a2(i,j))*uy(j)
     &                +(k23(i,j)+a2a3(i,j))*uz(j))			   !-k23(i,j)+a2a3(i,j))*uz(j)*s2*eff0
     &            +s3*((k31(i,j)+a3a1(i,j))*ux(j)                          !-k31(i,j)+a3a1(i,j))*ux(j)*s3*eff0
     &                +(k32(i,j)+a3a2(i,j))*uy(j)			   !-k32(i,j)+a3a2(i,j))*uy(j)*s3*eff0
     &                +(k33(i,j)+a3a3(i,j))*uz(j)))     *eff0		   !-k33(i,j)+a3a3(i,j))*uz(j)*s3*eff0
                ! where is the fifth term
                ! last term
     &          - ident(i,j)*(s1*ux(j)+s2*uy(j)+s3*uz(j))*effs             !-ident(i,j)*s3*uz(j)*effs

640             continue

200         continue
            enddo       !quadrature loop ends


!.......Implementing NonLinear BND condition...........................
!.......Non-Linear BND Condition Implementation END....................


!.......Implementing NonLinear BND condition for Wall Temperature
!.......IF AdiabaticWall=.FALSE. then Wall Temp is specified.
            if(.not.AdiabaticWall)then
                do iface=1,neface
                    irng=rng(iface,ie)
                    write(19,*) 'ie,iface,irng',ie,iface,irng
!.......check if it is Wall (Bnd 6-18)
!.......Condition on Wall Temperature Implemented
!.......WallTemp = Read from input file
                    if(irng.eq.1) then
                        do jnode=1,nnface
                            inode=facemap(iface,jnode,elem_type)
                            if(inode.ne.0) then
                                idof=ndf
                                nai=(inen-1)*ndf+idof
!.......Wall Temp is Specified .....use it directly for replacing equations
                                p(nai)=-d1(idof,inode)+
     $                              cv*WallTemp*d1(1,inode)
     $                        +0.5d0/d1(1,inode)*(d1(2,inode)*d1(2,inode)
     $                        +d1(3,inode)*d1(3,inode)
     $                        +d1(4,inode)*d1(4,inode))
                                do nb=1,nen
                                    do j =1,ndf
                                        nbj = (nb-1)*ndf + j
                                        s(nai,nbj)=0.0d0
                                    end do
                                end do
                                na1=(inen-1)*ndf+1
                                na2=(inen-1)*ndf+2
                                na3=(inen-1)*ndf+3
                                na4=(inen-1)*ndf+4
                                na5=(inen-1)*ndf+5
                               s(nai,na1)=-cv*WallTemp +
     $                              0.5d0/d1(1,inode)/d1(1,inode)*
     $                              (d1(2,inode)*d1(2,inode)
     $                              +d1(3,inode)*d1(3,inode)
     $                              +d1(4,inode)*d1(4,inode))
                                s(nai,na2)=-d1(2,inode)/d1(1,inode)
                                s(nai,na3)=-d1(3,inode)/d1(1,inode)
                                s(nai,na4)=-d1(4,inode)/d1(1,inode)
                                s(nai,na5)=1.0d0
                            end if 
                        end do
                    endif
                enddo
            endif
!.......Non-Linear BND Condition for Wall Temperature Implementation END

!.......Rotating the flow LOCALLY to implement SLIP BND CONDITION.
!.......REFER: Thesis by Aliabadi or LE BEAU
!.......Rotate the entire Element Level Matrix.
!.......Implementing SLIP BND CONDITION START..........................
!.......FIXME this is a little tricky
            flag=.FALSE.
            do iface=1,neface
                irng=rng(iface,ie)
                if(irng.eq.9) then
                    flag=.true.
                    exit
                endif
            enddo

!!            flag=.FALSE.

            if(flag)then
                do ia=1,nen
                    angle_a = 0.0 !angle(ien(ia,ie))
!!                    write(125,*) ien(ia,ie),angle_a
                    cos_a = cos(angle_a)
                    sin_a = sin(angle_a)
                    ia1 = (ia-1)*ndf+1
                    ia2 = ia1+1
                    ia3 = ia2+1
                    ia4 = ia3+1
                    ia5 = ia4+1
                    rotate(ia1,ia1) = 1.0d0
                    rotate(ia2,ia2) = cos_a
                    rotate(ia2,ia3) =-sin_a
                    rotate(ia3,ia2) = sin_a
                    rotate(ia3,ia3) = cos_a
                    rotate(ia4,ia4) = 1.0d0
                    rotate(ia5,ia5) = 1.0d0
                end do
                do ia=1,nee
                    do ip=1,nee
                        do ib=1,nee
                            do iq=1,nee
                                selm(ia,ib) = selm(ia,ib)
     $                          + rotate(ip,ia)*s(ip,iq)*rotate(iq,ib)
                            enddo
                        enddo
                        p_rot(ia) = p_rot(ia) + rotate(ip,ia)*p(ip)
                    enddo
                enddo
            else
                do ia=1,nee
                    do ib=1,nee
                        selm(ia,ib) = s(ia,ib)
                    enddo
                    p_rot(ia) = p(ia)
                enddo
            end if

!...........assemble the arrays.

!...........block diagonal assembly
!!            do inen=1,nen
!!                inode=ien(inen,ie)
!!                do idof=1,ndf
!!                    ia=(inen-1)*ndf+idof
!!                    bp1(idof,inode) = bp1(idof,inode) + p_rot(ia)
!!                    do jdof=1,ndf
!!                        ib=(inen-1)*ndf+jdof
!!                        block(jdof,idof,inode) =
!!     $                      block(jdof,idof,inode) + selm(ib,ia,ie)
!!                    enddo
!!                enddo
!!            enddo

!.............diagonal assembly
            do inen=1,nen
                inode=ien(inen,ie)
                do idof=1,ndf
                    ia=(inen-1)*ndf+idof
                    bp1(idof,inode) = bp1(idof,inode) + p_rot(ia)
                    wg1(idof,inode) = wg1(idof,inode) + selm(ia,ia)
                enddo
            enddo

!...........section for data verification
            icode = 0

            if(icode.eq.1) then
!                write(100+my_pe,*) ((selm(i,j,ie),i=1,nee),j=1,nee)
                write(105+my_pe,*) (p_rot(i),i=1,nee)
            else if(icode.eq.2) then
                do inen=1,nen
                    inode=ien(inen,ie)
!!                    write(100+my_pe,*) 
!!     $                  ((block(i,j,inode),i=1,ndf),j=1,ndf)
                    write(105+my_pe,*)
     $                  (bp1(idof,inode),idof=1,ndf)
                enddo
            else if(icode.eq.3) then
                if(ie.eq.nec) then
                    do inode=1,nnl
!!                        write(90+my_pe,*)
!!     $                   ((block(i,j,inode),i=1,ndf),j=1,ndf)
                        write(120+my_pe,*) (bp1(idof,inode),idof=1,ndf)
                    enddo
                endif
            endif

!.......ELEMENT LOOP ENDS...............................!
        enddo   

        return
        endsubroutine
