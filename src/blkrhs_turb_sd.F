	subroutine blkrhs_turb_sd(ien,xg1,dg1,dg0,pgt,wgt,
     &                    dg1_t,dg0_t,ynode,
     $			  umag_trans,he_trans,
     &                    dist_trans,vort_trans,
     $			  ntrans)

              
c**************************************************************************
c	BLOCK FOR 3D INCOMPRESSIBLE SEMI-DISCRETE GAL/LS METHOD
c       TO COMPUTE RESIDUALS FOR MATRIX FREE GMRES TECHNIQUE. 
c	LOOPS HAVE BEEN UNROLLED TO ACHIEVE SPEED UP
c	WRITTEN BY SANJAY MITTAL December 31, 2012
c****************************************************************************
	implicit none
	include "global.h"
	include "turbpara.h"
C...EXTERNAL VARIABLES
	integer ien(nen,nec)
	real*8  xg1(nsd,nnl),dg1(ndf,nnl) 
        real*8  wgt(ndf_t,nnl)
	real*8  dg0(ndf,nnl)
        real*8  dg1_t(ndf_t,nnl),dg0_t(ndf_t,nnl)
	real*8  pgt(ndf_t,nnl)

	integer ntrans
        real*8  umag_trans(nnl),he_trans(nnl)
	real*8  dist_trans(nnl),vort_trans(nnl)
        real*8  ynode(nnl)

C...INTERNAL VARIABLES
	integer mgau2 
	parameter(mgau2=21)

	real*8 x1(nsd,nen)

	real*8 lturb(nen)
	real*8 vis,dif_turb

	real*8 sh(0:3,nen),det
	real*8 xs(nsd,nsd), gaus(nsd,mgau2), weis(mgau2)

	real*8 do(ndf,nen), d1(ndf,nen)

        real*8 do_t(ndf_t,nen),d1_t(ndf_t,nen)
        real*8 p_t(ndf_t,nen),w_t(ndf_t,nen)

        real*8 distl_t(nen),vortl_t(nen),hel_t(nen),umagl_t(nen)

	real*8 yt(nen)
        real*8 distance
	real*8 vort_trip,he_trip,umag_trip,dist_trip

        real*8 nutt,nutx,nuty,nutz,nut
        real*8 nut_n,nutx_n,nuty_n,nutz_n
        real*8 nut_np1,nutx_np1,nuty_np1,nutz_np1
	real*8 nutxx_np1,nutyy_np1,nutzz_np1,nutxy_np1
	real*8 nutyz_np1,nutxz_np1
	real*8 nutxx_n,nutyy_n,nutzz_n,nutxy_n
	real*8 nutyz_n,nutxz_n
	real*8 nutxx,nutyy,nutzz,nutxy,nutyz,nutxz

	real*8 grad_nut,gradnu_sq

        real*8 u(nsd),un(nsd),unp1(nsd),umag,umagn

	real*8 rho_np1, rhou_np1,rhov_np1,rhow_np1,rhoe_np1
	real*8 rho_n  , rhou_n  ,rhov_n  ,rhow_n  ,rhoe_n
	real*8          rhou    ,rhov    , rhow   , rhoe

	real*8 c_sound,umagc

	real*8 rhoux_np1,rhouy_np1,rhouz_np1
	real*8 rhovx_np1,rhovy_np1,rhovz_np1
	real*8 rhowx_np1,rhowy_np1,rhowz_np1
	real*8 rhoux_n  ,rhouy_n  ,rhouz_n  
	real*8 rhovx_n  ,rhovy_n  ,rhovz_n
	real*8 rhowx_n  ,rhowy_n  ,rhowz_n
	real*8 rhoux    ,rhouy    ,rhouz
	real*8 rhovx    ,rhovy    ,rhovz
	real*8 rhowx    ,rhowy    ,rhowz

	real*8 uy_np1,uz_np1,vx_np1,vz_np1,wx_np1,wy_np1
	real*8 uy_n  ,uz_n  ,vx_n  ,vz_n  ,wx_n  ,wy_n
	real*8 uy    ,uz    ,vx    ,vz    ,wx    ,wy

	real*8 umag0,umagx,umagy,umagz,grad_umag
        real*8 rx,ry,rz,sx,sy,sz
        real*8 he_min,he_max,he_rgn,he_ugn,he_rgnt
        real*8 dis1,dis2,dis3,dis4,dis5,dis6
	real*8 temp1, temp2,temp3,temp4,temp5,temp6
	real*8 temp7, temp8,temp9,temp10,temp11,temp12
	real*8 temp13, temp14,temp15,temp16,temp17,temp18

        real*8  h1,h2,h3,hel,l1,w1,area,vol
	real*8 tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8
	real*8 sgrad,fv1,fv2,sbar,udif

	real*8 gt,fac1,fac2,t1,ft1,ft2,r,g,fw

	real*8 p_function,q_function,r_function
	real*8 tau_supg
	real*8 eff0
	real*8 b0,b1,b2,b3,b4,b5,b6,ub
	real*8 dif,dif_mol
	real*8 Rcon,Cv
c......ARRAYS FOR SHAPE ROUTINES
        real*8   xr(nen,nen),cf(nen,nen)
        real*8   cf11,cf21,cf31,cf41,cf12,cf22,cf32,cf42
        real*8   cf13,cf23,cf33,cf43,cf14,cf24,cf34,cf44
        real*8   xr11,xr21,xr31,xr41,xr12,xr22,xr32,xr42
        real*8   xr13,xr23,xr33,xr43,xr14,xr24,xr34,xr44
        integer na,nb,ig,inl,knl,jnl,k,isd,jsd
        integer imat,iee,inode,level,ie,idf,i,j
!**************************************************************************

	Rcon = 1.0d0/gamma
	Cv= Rcon/(gamma-1.0d0)

	call lgau3t(iquad,nquad3,gaus,weis)   !Quadrature Points

        dif = viscosity       ! free stream dynamic viscosity (mu)
        dif_mol = dif

!	write(6,*)"ndf ndf_t", ndf,ndf_t
!	write(6,*)"dif", dif
!	write(6,*) "nquad", nquad3

C...CLEAR BUFFER ARRAY
	call fclear(pgt,nnl*ndf_t)
	call fclear(wgt,nnl*ndf_t)

!	turb_kappa = 0.15
!	dif = 1.0/reynolds

c...3D TERMS

	do ie = 1,nec

C...CLEAR LOCAL ARRAYS
		do inl = 1,nen
			do idf = 1,ndf_t
				p_t(idf,inl) = 0.0
				w_t(idf,inl) = 0.0
			enddo
		enddo

C...LOCALIZE  ARRAYS
               	do inl = 1,nen 
			k = ien(inl,ie)
       	               		do idf = 1,ndf
                       		do(idf,inl) = dg0(idf,k)
                               	d1(idf,inl) = dg1(idf,k)
                       	enddo

                       	do isd = 1,nsd
                               	x1(isd,inl)  = xg1(isd,k)
                       	enddo

c..........LOCALIZE EDDY VISCOSITY
				do idf = 1,ndf_t
                       		do_t(idf,inl) = dg0_t(idf,k)
				d1_t(idf,inl) = dg1_t(idf,k)
!				write(6,*) do_t(idf,inl),d1_t(idf,inl)
                       		enddo

C...........LOCALIZE DISTANCE FROM WALL
			yt(inl) = ynode(k)

			
C...........LOCALIZE TRANSITION TERMS
!	check murali 
                        distl_t(inl)= dist_trans(k)
                        vortl_t(inl)= vort_trans(k)
                          hel_t(inl)= he_trans(k)
                        umagl_t(inl)= umag_trans(k) 
!	check murali
                enddo

		do ig=1,nquad3
                call sf36o(gaus(1,ig),x1,det,sh,1,xs)
                det=dabs(det)

!			if(nen.eq.8) then
!				include "sf3d8n.h"
!			elseif(nen.eq.6)then
!				include "sf3d6n.h"
!			else
!				include "sf3d4n.h"
!			endif

!	Distance of integration point from the wall
	  distance = sh(0,1)*yt(1) + sh(0,2)*yt(2) + sh(0,3)*yt(3)
     $             + sh(0,4)*yt(4) + sh(0,5)*yt(5) + sh(0,6)*yt(6)

!	 write(6,*) "distance", distance

!.......Distance of integration point from trip point

	  if(trip) then
	  dist_trip = sh(0,1)*distl_t(1) + sh(0,2)*distl_t(2)
     &              + sh(0,3)*distl_t(3) + sh(0,4)*distl_t(4) 
     &              + sh(0,5)*distl_t(5) + sh(0,6)*distl_t(6)

!.......Vorticity at trip
	  vort_trip = sh(0,1)*vortl_t(1) + sh(0,2)*vortl_t(2)
     &              + sh(0,3)*vortl_t(3) + sh(0,4)*vortl_t(4)
     &              + sh(0,5)*vortl_t(5) + sh(0,6)*vortl_t(6)

!.......Velocity magnitude at trip
	  umag_trip = sh(0,1)*umagl_t(1) + sh(0,2)*umagl_t(2)
     &              + sh(0,3)*umagl_t(3) + sh(0,4)*umagl_t(4)
     &              + sh(0,5)*umagl_t(5) + sh(0,6)*umagl_t(6)


!.......Element length (along the wall) at trip
	  he_trip = sh(0,1)*hel_t(1) + sh(0,2)*hel_t(2)
     &            + sh(0,3)*hel_t(3) + sh(0,4)*hel_t(4)
     &            + sh(0,5)*hel_t(5) + sh(0,6)*hel_t(6)
	   end if


	  nut_np1 =  sh(0,1)*d1_t(1,1) + sh(0,2)*d1_t(1,2)
     &             + sh(0,3)*d1_t(1,3) + sh(0,4)*d1_t(1,4)
     &             + sh(0,5)*d1_t(1,5) + sh(0,6)*d1_t(1,6)

	  nut_n =    sh(0,1)*do_t(1,1) + sh(0,2)*do_t(1,2)
     &             + sh(0,3)*do_t(1,3) + sh(0,4)*do_t(1,4)
     &             + sh(0,5)*do_t(1,5) + sh(0,6)*do_t(1,6)
	nut        = alpha*nut_np1 + oma*nut_n

	  nutx_np1 = sh(1,1)*d1_t(1,1) + sh(1,2)*d1_t(1,2)
     &             + sh(1,3)*d1_t(1,3) + sh(1,4)*d1_t(1,4)
     &             + sh(1,5)*d1_t(1,5) + sh(1,6)*d1_t(1,6)

	  nuty_np1 = sh(2,1)*d1_t(1,1) + sh(2,2)*d1_t(1,2)
     &             + sh(2,3)*d1_t(1,3) + sh(2,4)*d1_t(1,4)
     &             + sh(2,5)*d1_t(1,5) + sh(2,6)*d1_t(1,6)

	  nutz_np1 = sh(3,1)*d1_t(1,1) + sh(3,2)*d1_t(1,2)
     &             + sh(3,3)*d1_t(1,3) + sh(3,4)*d1_t(1,4)
     &             + sh(3,5)*d1_t(1,5) + sh(3,6)*d1_t(1,6)

	  nutx_n   = sh(1,1)*do_t(1,1) + sh(1,2)*do_t(1,2)
     &             + sh(1,3)*do_t(1,3) + sh(1,4)*do_t(1,4)
     &             + sh(1,5)*do_t(1,5) + sh(1,6)*do_t(1,6)

	  nuty_n   = sh(2,1)*do_t(1,1)  + sh(2,2)*do_t(1,2)
     &             + sh(2,3)*do_t(1,3)  + sh(2,4)*do_t(1,4)
     &             + sh(2,5)*do_t(1,5)  + sh(2,6)*do_t(1,6)

	  nutz_n   = sh(3,1)*do_t(1,1)  + sh(3,2)*do_t(1,2)
     &             + sh(3,3)*do_t(1,3)  + sh(3,4)*do_t(1,4)
     &             + sh(3,5)*do_t(1,5)  + sh(3,6)*do_t(1,6)


	  nutx  = alpha * nutx_np1 + oma * nutx_n 
	  nuty  = alpha * nuty_np1 + oma * nuty_n
	  nutz  = alpha * nutz_np1 + oma * nutz_n

	  grad_nut = dsqrt(nutx_n*nutx_n+nuty_n*nuty_n+nutz_n*nutz_n)

	  rho_np1 = sh(0,1)*d1(1,1) + sh(0,2)*d1(1,2)
     &            + sh(0,3)*d1(1,3) + sh(0,4)*d1(1,4)
     &            + sh(0,5)*d1(1,5) + sh(0,6)*d1(1,6)


	  rhou_np1= sh(0,1)*d1(2,1) + sh(0,2)*d1(2,2)
     &            + sh(0,3)*d1(2,3) + sh(0,4)*d1(2,4)
     &            + sh(0,5)*d1(2,5) + sh(0,6)*d1(2,6)


	  rhov_np1= sh(0,1)*d1(3,1) + sh(0,2)*d1(3,2)
     &            + sh(0,3)*d1(3,3) + sh(0,4)*d1(3,4)
     &            + sh(0,5)*d1(3,5) + sh(0,6)*d1(3,6)


	  rhow_np1= sh(0,1)*d1(4,1) + sh(0,2)*d1(4,2)
     &            + sh(0,3)*d1(4,3) + sh(0,4)*d1(4,4)
     &            + sh(0,5)*d1(4,5) + sh(0,6)*d1(4,6)


	  rhoe_np1= sh(0,1)*d1(5,1) + sh(0,2)*d1(5,2)
     &            + sh(0,3)*d1(5,3) + sh(0,4)*d1(5,4)
     &            + sh(0,5)*d1(5,5) + sh(0,6)*d1(5,6)


	  rho_n   = sh(0,1)*do(1,1) + sh(0,2)*do(1,2)
     &            + sh(0,3)*do(1,3) + sh(0,4)*do(1,4)
     &            + sh(0,5)*do(1,5) + sh(0,6)*do(1,6)

	  rhou_n  = sh(0,1)*do(2,1) + sh(0,2)*do(2,2)
     &            + sh(0,3)*do(2,3) + sh(0,4)*do(2,4)
     &            + sh(0,5)*do(2,5) + sh(0,6)*do(2,6)

	  rhov_n  = sh(0,1)*do(3,1) + sh(0,2)*do(3,2)
     &            + sh(0,3)*do(3,3) + sh(0,4)*do(3,4)
     &            + sh(0,5)*do(3,5) + sh(0,6)*do(3,6)

	  rhow_n  = sh(0,1)*do(4,1) + sh(0,2)*do(4,2)
     &            + sh(0,3)*do(4,3) + sh(0,4)*do(4,4)
     &            + sh(0,5)*do(4,5) + sh(0,6)*do(4,6)


	  rhoe_n  = sh(0,1)*do(5,1) + sh(0,2)*do(5,2)
     &            + sh(0,3)*do(5,3) + sh(0,4)*do(5,4)
     &		  + sh(0,5)*do(5,5) + sh(0,6)*do(5,6)

	  unp1(1) = rhou_np1/rho_np1
	  unp1(2) = rhov_np1/rho_np1
	  unp1(3) = rhow_np1/rho_np1
	  un(1)   = rhou_n/rho_n
	  un(2)   = rhov_n/rho_n
	  un(3)   = rhow_n/rho_n

	  u(1) = alpha * unp1(1) + oma * un(1)
	  u(2) = alpha * unp1(2) + oma * un(2)
	  u(3) = alpha * unp1(3) + oma * un(3)

! 	write(6,*) "u  v   w ", u(1), u(2), u(3)
	  rho  = alpha * rho_np1 + oma * rho_n
	  rhou = alpha * rhou_np1+ oma * rhou_n
	  rhov = alpha * rhov_np1+ oma * rhov_n
	  rhow = alpha * rhow_np1+ oma * rhow_n
	  rhoe = alpha * rhoe_np1+ oma * rhoe_n

	  umag  = dsqrt(u(1)*u(1) + u(2)*u(2) + u(3)*u(3))        !used in he_ugn
	  umagn = dsqrt(un(1)*un(1) + un(2)*un(2) + un(3)*un(3))   !used in tau_supg
	  c_sound = sqrt((rhoe - 0.50d0*
     $	  (rhou*rhou+rhov*rhov+rhow*rhow)/rho)/rho/Cv)
	  umagc = umagn + c_sound
!.......NOTE: 'umagc' may be used in place of 'umagn' in tau_supg

!	  rhox_np1 = sh(1,1)*d1(1,1) + sh(1,2)*d1(1,2)
!     &             + sh(1,3)*d1(1,3) 
!	  rhoy_np1 = sh(2,1)*d1(1,1) + sh(2,2)*d1(1,2)
!     &             + sh(2,3)*d1(1,3) 
	  rhoux_np1= sh(1,1)*d1(2,1) + sh(1,2)*d1(2,2)
     &             + sh(1,3)*d1(2,3) + sh(1,4)*d1(2,4) 
     &             + sh(1,5)*d1(2,5) + sh(1,6)*d1(2,6)
 
!	  rhouy_np1= sh(2,1)*d1(2,1) + sh(2,2)*d1(2,2)
!     &             + sh(2,3)*d1(2,3) 
!	  rhovx_np1= sh(1,1)*d1(3,1) + sh(1,2)*d1(3,2)
!     &             + sh(1,3)*d1(3,3) 
	  rhouy_np1= sh(2,1)*d1(2,1) + sh(2,2)*d1(2,2)
     &             + sh(2,3)*d1(2,3) + sh(2,4)*d1(2,4) 
     &             + sh(2,5)*d1(2,5) + sh(2,6)*d1(2,6)
 
	  rhouz_np1= sh(3,1)*d1(2,1) + sh(3,2)*d1(2,2)
     &             + sh(3,3)*d1(2,3) + sh(3,4)*d1(2,4) 
     &             + sh(3,5)*d1(2,5) + sh(3,6)*d1(2,6)

	  rhovx_np1= sh(1,1)*d1(3,1) + sh(1,2)*d1(3,2)
     &             + sh(1,3)*d1(3,3) + sh(1,4)*d1(3,4)
     &             + sh(1,5)*d1(3,5) + sh(1,6)*d1(3,6)

	  rhovy_np1= sh(2,1)*d1(3,1) + sh(2,2)*d1(3,2)
     &             + sh(2,3)*d1(3,3) + sh(2,4)*d1(3,4)
     &             + sh(2,5)*d1(3,5) + sh(2,6)*d1(3,6)

	  rhovz_np1= sh(3,1)*d1(3,1) + sh(3,2)*d1(3,2)
     &             + sh(3,3)*d1(3,3) + sh(3,4)*d1(3,4)
     &             + sh(3,5)*d1(3,5) + sh(3,6)*d1(3,6)

	  rhowx_np1= sh(1,1)*d1(4,1) + sh(1,2)*d1(4,2)
     &             + sh(1,3)*d1(4,3) + sh(1,4)*d1(4,4)
     &             + sh(1,5)*d1(4,5) + sh(1,6)*d1(4,6)

	  rhowy_np1= sh(2,1)*d1(4,1) + sh(2,2)*d1(4,2)
     &             + sh(2,3)*d1(4,3) + sh(2,4)*d1(4,4)
     &             + sh(2,5)*d1(4,5) + sh(2,6)*d1(4,6)

	  rhowz_np1= sh(3,1)*d1(4,1) + sh(3,2)*d1(4,2)
     &             + sh(3,3)*d1(4,3) + sh(3,4)*d1(4,4)
     &             + sh(3,5)*d1(4,5) + sh(3,6)*d1(4,6)

	  rhoux_n  = sh(1,1)*do(2,1) + sh(1,2)*do(2,2)
     &             + sh(1,3)*do(2,3) + sh(1,4)*do(2,4) 
     &             + sh(1,5)*do(2,5) + sh(1,6)*do(2,6)

	  rhouy_n  = sh(2,1)*do(2,1) + sh(2,2)*do(2,2)
     &             + sh(2,3)*do(2,3) + sh(2,4)*do(2,4) 
     &             + sh(2,5)*do(2,5) + sh(2,6)*do(2,6)
  
	  rhouz_n  = sh(3,1)*do(2,1) + sh(3,2)*do(2,2)
     &             + sh(3,3)*do(2,3) + sh(3,4)*do(2,4) 
     &             + sh(3,5)*do(2,5) + sh(3,6)*do(2,6) 

	  rhovx_n = sh(1,1)*do(3,1)  + sh(1,2)*do(3,2)
     &             + sh(1,3)*do(3,3) + sh(1,4)*do(3,4)
     &             + sh(1,5)*do(3,5) + sh(1,6)*do(3,6)

	  rhovy_n  = sh(2,1)*do(3,1) + sh(2,2)*do(3,2)
     &             + sh(2,3)*do(3,3) + sh(2,4)*do(3,4)
     &             + sh(2,5)*do(3,5) + sh(2,6)*do(3,6)

	  rhovz_n  = sh(3,1)*do(3,1) + sh(3,2)*do(3,2)
     &             + sh(3,3)*do(3,3) + sh(3,4)*do(3,4)
     &             + sh(3,5)*do(3,5) + sh(3,6)*do(3,6)

	  rhowx_n  = sh(1,1)*do(4,1) + sh(1,2)*do(4,2)
     &             + sh(1,3)*do(4,3) + sh(1,4)*do(4,4)
     &             + sh(1,5)*do(4,5) + sh(1,6)*do(4,6)

	  rhowy_n  = sh(2,1)*do(4,1) + sh(2,2)*do(4,2)
     &             + sh(2,3)*do(4,3) + sh(2,4)*do(4,4)
     &             + sh(2,5)*do(4,5) + sh(2,6)*do(4,6)

	  rhowz_n  = sh(3,1)*do(4,1) + sh(3,2)*do(4,2)
     &             + sh(3,3)*do(4,3) + sh(3,4)*do(4,4)
     &             + sh(3,5)*do(4,5) + sh(3,6)*do(4,6)

	  rhoux = alpha * rhoux_np1+ oma * rhoux_n
	  rhouy = alpha * rhouy_np1+ oma * rhouy_n
	  rhouz = alpha * rhouz_np1+ oma * rhouz_n

	  rhovx = alpha * rhovx_np1+ oma * rhovx_n
	  rhovy = alpha * rhovy_np1+ oma * rhovy_n
	  rhovz = alpha * rhovz_np1+ oma * rhovz_n

	  rhowx = alpha * rhowx_np1+ oma * rhowx_n
	  rhowy = alpha * rhowy_np1+ oma * rhowy_n
	  rhowz = alpha * rhowz_np1+ oma * rhowz_n
!	  rhoux_n  = sh(1,1)*do(2,1) + sh(1,2)*do(2,2)
!     &             + sh(1,3)*do(2,3) + sh(1,4)*do(2,4)
!     &             + sh(1,5)*do(2,5) + sh(1,6)*do(2,6)
!
!	  rhouz_n  = sh(3,1)*do(2,1) + sh(3,2)*do(2,2)
!     &             + sh(3,3)*do(2,3) + sh(3,4)*do(2,4)
!     &             + sh(3,5)*do(2,5) + sh(3,6)*do(2,6)
!
!	  rhovx_n  = sh(2,1)*do(3,1) + sh(2,2)*do(3,2)
!     &             + sh(2,3)*do(3,3) + sh(2,4)*do(3,4)
!     &             + sh(2,5)*do(3,5) + sh(2,6)*do(3,6)
!
!	  rhovz_n  = sh(3,1)*do(3,1) + sh(3,2)*do(3,2)
!     &             + sh(3,3)*do(3,3) + sh(3,4)*do(3,4)
!     &             + sh(3,5)*do(3,5) + sh(3,6)*do(3,6)

!	  rhoux = alpha * rhoux_np1+ oma * rhouy_n
!	  rhovy = alpha * rhovy_np1+ oma * rhovy_n

	  uy_np1=sh(2,1)*d1(2,1)/d1(1,1)+sh(2,2)*d1(2,2)/d1(1,2)
     &          +sh(2,3)*d1(2,3)/d1(1,3)+sh(2,4)*d1(2,4)/d1(1,4)
     &          +sh(2,5)*d1(2,5)/d1(1,5)+sh(2,6)*d1(2,6)/d1(1,6)

	  uz_np1=sh(3,1)*d1(2,1)/d1(1,1)+sh(3,2)*d1(2,2)/d1(1,2)
     &          +sh(3,3)*d1(2,3)/d1(1,3)+sh(3,4)*d1(2,4)/d1(1,4)
     &          +sh(3,5)*d1(2,5)/d1(1,5)+sh(3,6)*d1(2,6)/d1(1,6)

	  vx_np1=sh(1,1)*d1(3,1)/d1(1,1)+sh(1,2)*d1(3,2)/d1(1,2)
     &          +sh(1,3)*d1(3,3)/d1(1,3)+sh(1,4)*d1(3,4)/d1(1,4)
     &          +sh(1,5)*d1(3,5)/d1(1,5)+sh(1,6)*d1(3,6)/d1(1,6)

	  vz_np1=sh(3,1)*d1(3,1)/d1(1,1)+sh(3,2)*d1(3,2)/d1(1,2)
     &          +sh(3,3)*d1(3,3)/d1(1,3)+sh(3,4)*d1(3,4)/d1(1,4)
     &          +sh(3,5)*d1(3,5)/d1(1,5)+sh(3,6)*d1(3,6)/d1(1,6)

	  wx_np1=sh(1,1)*d1(4,1)/d1(1,1)+sh(1,2)*d1(4,2)/d1(1,2)
     &          +sh(1,3)*d1(4,3)/d1(1,3)+sh(1,4)*d1(4,4)/d1(1,4)
     &          +sh(1,5)*d1(4,5)/d1(1,5)+sh(1,6)*d1(4,6)/d1(1,6)

	  wy_np1=sh(2,1)*d1(4,1)/d1(1,1)+sh(2,2)*d1(4,2)/d1(1,2)
     &          +sh(2,3)*d1(4,3)/d1(1,3)+sh(2,4)*d1(4,4)/d1(1,4)
     &          +sh(2,5)*d1(4,5)/d1(1,5)+sh(2,6)*d1(4,6)/d1(1,6)

	  uy_n  =sh(2,1)*do(2,1)/do(1,1)+sh(2,2)*do(2,2)/do(1,2)
     &          +sh(2,3)*do(2,3)/do(1,3)+sh(2,4)*do(2,4)/do(1,4)
     &          +sh(2,5)*do(2,5)/do(1,5)+sh(2,6)*do(2,6)/do(1,6)

	  uz_n  =sh(3,1)*do(2,1)/do(1,1)+sh(3,2)*do(2,2)/do(1,2)
     &          +sh(3,3)*do(2,3)/do(1,3)+sh(3,4)*do(2,4)/do(1,4)
     &          +sh(3,5)*do(2,5)/do(1,5)+sh(3,6)*do(2,6)/do(1,6)

	  vx_n  =sh(1,1)*do(3,1)/do(1,1)+sh(1,2)*do(3,2)/do(1,2)
     &          +sh(1,3)*do(3,3)/do(1,3)+sh(1,4)*do(3,4)/do(1,4)
     &          +sh(1,5)*do(3,5)/do(1,5)+sh(1,6)*do(3,6)/do(1,6)

	  vz_n  =sh(3,1)*do(3,1)/do(1,1)+sh(3,2)*do(3,2)/do(1,2)
     &          +sh(3,3)*do(3,3)/do(1,3)+sh(3,4)*do(3,4)/do(1,4)
     &          +sh(3,5)*do(3,5)/do(1,5)+sh(3,6)*do(3,6)/do(1,6)

	  wx_n  =sh(1,1)*do(4,1)/do(1,1)+sh(1,2)*do(4,2)/do(1,2)
     &          +sh(1,3)*do(4,3)/do(1,3)+sh(1,4)*do(4,4)/do(1,4)
     &          +sh(1,5)*do(4,5)/do(1,5)+sh(1,6)*do(4,6)/do(1,6)

	  wy_n  =sh(2,1)*do(4,1)/do(1,1)+sh(2,2)*do(4,2)/do(1,2)
     &          +sh(2,3)*do(4,3)/do(1,3)+sh(2,4)*do(4,4)/do(1,4)
     &          +sh(2,5)*do(4,5)/do(1,5)+sh(2,6)*do(4,6)/do(1,6)

	  uy = alpha * uy_np1 + oma * uy_n
	  uz = alpha * uz_np1 + oma * uz_n
	  vx = alpha * vx_np1 + oma * vx_n
	  vz = alpha * vz_np1 + oma * vz_n
	  wx = alpha * wx_np1 + oma * wx_n
	  wy = alpha * wy_np1 + oma * wy_n


	 goto 101


	      if(flag_he.eq.1)then
                if(umag.gt.1e-8)then
                   h1 = u(1)/umag  !cos(theta), theta=angle of vel vector with hor.
                   h2 = u(2)/umag  !sin(theta) 
		   h3 = u(3)/umag

                 hel =   abs(h1*sh(1,1)+h2*sh(2,1)+h3*sh(3,1))
     $                 + abs(h1*sh(1,2)+h2*sh(2,2)+h3*sh(3,2))
     $                 + abs(h1*sh(1,3)+h2*sh(2,3)+h3*sh(3,3))
     $                 + abs(h1*sh(1,4)+h2*sh(2,4)+h3*sh(3,4))
     $                 + abs(h1*sh(1,5)+h2*sh(2,5)+h3*sh(3,5))
     $                 + abs(h1*sh(1,6)+h2*sh(2,6)+h3*sh(3,6))

                   hel = 0.5d0/hel
                else
                   hel=min(
     &             sqrt((x1(1,2)-x1(1,1))**2+(x1(2,2)-x1(2,1))**2),
     $             sqrt((x1(1,3)-x1(1,1))**2+(x1(2,3)-x1(2,1))**2),
     $             sqrt((x1(1,3)-x1(1,2))**2+(x1(2,3)-x1(2,2))**2))
!     $		    abs(x1(3,4)-x1(3,1)))
                endif
		he_ugn = hel
		he_rgn = hel

              else if(flag_he.eq.2) then
!.......MY Suggestion: Length based on AVERAGE of lengths of sides
                  hel=(
     &           sqrt((x1(1,2)-x1(1,1))**2+(x1(2,2)-x1(2,1))**2)+
     $           sqrt((x1(1,3)-x1(1,1))**2+(x1(2,3)-x1(2,1))**2)+
     $           sqrt((x1(1,3)-x1(1,2))**2+(x1(2,3)-x1(2,2))**2)+
     $		 abs(x1(3,4)-x1(3,1))) /4.0
		he_ugn = hel
		he_rgn = hel

              else if(flag_he.eq.3) then
!.......Element Length(hel) based on Length of Edge of Triangle
!.......Based on MINIMUM length
                 hel=dmin1(
     &           sqrt((x1(1,2)-x1(1,1))**2+(x1(2,2)-x1(2,1))**2),
     $           sqrt((x1(1,3)-x1(1,1))**2+(x1(2,3)-x1(2,1))**2),
     $           sqrt((x1(1,3)-x1(1,2))**2+(x1(2,3)-x1(2,2))**2))
!     $		 sqrt((x1(1,4)-x1(1,1))**2+(x1(3,4)-x1(3,1))**2))

		he_ugn = hel
		he_rgn = hel

		else if(flag_he.eq.4) then
		l1 = 0.5*((x1(1,2)*x1(2,3)-x1(1,3)*x1(2,2))
     $                  - (x1(1,1)*x1(2,3)-x1(1,3)*x1(2,1))
     $                  + (x1(1,1)*x1(2,2)-x1(1,2)*x1(2,1)))
		w1 = x1(3,4)-x1(3,1)
		vol = abs(l1*w1)
		tmp1 = (x1(1,1)-x1(1,2))*(x1(1,1)-x1(1,2)) 
     $		      +(x1(2,1)-x1(2,2))*(x1(2,1)-x1(2,2))
		tmp2 = (x1(3,1)-x1(3,4))*(x1(3,1)-x1(3,4))
		tmp3 = (x1(1,3)-x1(1,2))*(x1(1,3)-x1(1,2))   
     $                +(x1(2,3)-x1(2,2))*(x1(2,3)-x1(2,2))
		tmp4 = (x1(3,1)-x1(3,4))*(x1(3,1)-x1(3,4))
		tmp5 = (x1(1,3)-x1(1,1))*(x1(1,3)-x1(1,1))
     $                +(x1(2,3)-x1(2,1))*(x1(2,3)-x1(2,1))
		tmp6 = (x1(3,1)-x1(3,4))*(x1(3,1)-x1(3,4))
!		area = 2.0*abs(l1)   + sqrt(tmp1)*sqrt(tmp2)
!     $	       +sqrt(tmp3)*sqrt(tmp4)+ sqrt(tmp5)*sqrt(tmp6)

		area = 2.0*abs(l1)
     $		      +abs((x1(1,1)-x1(1,2))*(x1(3,1)-x1(3,4)))
     $                +abs((x1(2,3)-x1(2,2))*(x1(3,2)-x1(3,5)))
     $                +abs((x1(2,3)-x1(2,1))*(x1(3,3)-x1(3,6)))
		hel = vol/area
		he_ugn = hel
		he_rgn = hel

                else if(flag_he.eq.5) then
!.......Based on MAXIMUM length
!.......NOTE: For Fine MESH max length may be used
                 hel=dmax1(
     &            sqrt((x1(1,2)-x1(1,1))**2+(x1(2,2)-x1(2,1))**2),
     $            sqrt((x1(1,3)-x1(1,1))**2+(x1(2,3)-x1(2,1))**2),
     $            sqrt((x1(1,3)-x1(1,2))**2+(x1(2,3)-x1(2,2))**2))
!     $		   abs(x1(3,4)-x1(3,1)))
		he_ugn = hel
		he_rgn = hel

		else

	 temp1=sqrt(d1(1,1)*d1(1,1)+d1(2,1)*d1(2,1)+d1(3,1)*d1(3,1))
         temp2=sqrt(d1(1,2)*d1(1,2)+d1(2,2)*d1(2,2)+d1(3,2)*d1(3,2))
         temp3=sqrt(d1(1,3)*d1(1,3)+d1(2,3)*d1(2,3)+d1(3,3)*d1(3,3))
         temp4=sqrt(d1(1,4)*d1(1,4)+d1(2,4)*d1(2,4)+d1(3,4)*d1(3,4))
	 temp5=sqrt(d1(1,5)*d1(1,5)+d1(2,5)*d1(2,5)+d1(3,5)*d1(3,5))
         temp6=sqrt(d1(1,6)*d1(1,6)+d1(2,6)*d1(2,6)+d1(3,6)*d1(3,6))

	umag0=sh(0,1)*temp1+sh(0,2)*temp2+sh(0,3)*temp3+sh(0,4)*temp4
     &       +sh(0,5)*temp5+sh(0,6)*temp6
 
	umagx=sh(1,1)*temp1+sh(1,2)*temp2+sh(1,3)*temp3+sh(1,4)*temp4
     &       +sh(1,5)*temp5+sh(1,6)*temp6

	umagy=sh(2,1)*temp1+sh(2,2)*temp2+sh(2,3)*temp3+sh(2,4)*temp4
     &       +sh(2,5)*temp5+sh(2,6)*temp6

	umagz=sh(3,1)*temp1+sh(3,2)*temp2+sh(3,3)*temp3+sh(3,4)*temp4
     &       +sh(3,5)*temp5+sh(3,6)*temp6

       grad_umag=dsqrt(umagx*umagx + umagy*umagy + umagz*umagz)


	temp1=x1(1,1)-x1(1,2)
        temp2=x1(2,1)-x1(2,2)
        temp3=x1(3,1)-x1(3,2)
        dis1 = sqrt(temp1*temp1+temp2*temp2+temp3*temp3)
        temp1=x1(1,1)-x1(1,3)
        temp2=x1(2,1)-x1(2,3)
        temp3=x1(3,1)-x1(3,3)
        dis2 = sqrt(temp1*temp1+temp2*temp2+temp3*temp3)
        temp1=x1(1,2)-x1(1,3)
        temp2=x1(2,2)-x1(2,3)
        temp3=x1(3,2)-x1(3,3)
        dis3 = sqrt(temp1*temp1+temp2*temp2+temp3*temp3)

        he_max=0.2*max(dis1,dis2,dis3)                                                                          
C...GEOMETRIC ELEMENT LENGTH
!	fix it 
	he_min=min(dis1,dis2,dis3)
!	fix it 

	if (grad_umag .gt. 1.0e-9) then
             rx = nutx_n/grad_umag
             ry = nuty_n/grad_umag
             rz = nuty_n/grad_umag
                                                                                    
             he_rgn = dabs(rx*sh(1,1)+ry*sh(2,1)+rz*sh(3,1))
     &              + dabs(rx*sh(1,2)+ry*sh(2,2)+rz*sh(3,2))
     &              + dabs(rx*sh(1,3)+ry*sh(2,3)+rz*sh(3,3))
     &              + dabs(rx*sh(1,4)+ry*sh(2,4)+rz*sh(3,4))
     &              + dabs(rx*sh(1,5)+ry*sh(2,5)+rz*sh(3,5))
     &              + dabs(rx*sh(1,6)+ry*sh(2,6)+rz*sh(3,6))
                                                                                    
             he_rgn = 2.0/he_rgn
                                                                                    
            else
             he_rgn = he_min !he_max
        end if

	if (umag .gt. 1.0d-6) then
            sx=un(1)/umag
            sy=un(2)/umag
            sz=un(3)/umag
                                                                                    
            he_ugn = dabs(sx*sh(1,1)+sy*sh(2,1)+sz*sh(3,1))
     &             + dabs(sx*sh(1,2)+sy*sh(2,2)+sz*sh(3,2))
     &             + dabs(sx*sh(1,3)+sy*sh(2,3)+sz*sh(3,3))
     &             + dabs(sx*sh(1,4)+sy*sh(2,4)+sz*sh(3,4))
     &             + dabs(sx*sh(1,5)+sy*sh(2,5)+sz*sh(3,5))
     &             + dabs(sx*sh(1,6)+sy*sh(2,6)+sz*sh(3,6))
                                                                          
	   he_ugn = 2.0/he_ugn
            else
            he_ugn = he_min !he_max
       end if

	if (grad_nut .gt. 1.0e-9) then
            rx = nutx_n/grad_nut
            ry = nuty_n/grad_nut
            rz = nutz_n/grad_nut
                                                                                    
            he_rgnt= dabs(rx*sh(1,1)+ry*sh(2,1)+rz*sh(3,1))
     &             + dabs(rx*sh(1,2)+ry*sh(2,2)+rz*sh(3,2))
     &             + dabs(rx*sh(1,3)+ry*sh(2,3)+rz*sh(3,3))
     &             + dabs(rx*sh(1,4)+ry*sh(2,4)+rz*sh(3,4))
     &             + dabs(rx*sh(1,5)+ry*sh(2,5)+rz*sh(3,5))
     &             + dabs(rx*sh(1,6)+ry*sh(2,6)+rz*sh(3,6))
                                                                                    
            he_rgnt= 2.0/he_rgnt
                                                                                    
            else
            he_rgnt= he_min !he_max
        end if
	end if


101	continue

!	write(6,*) umag, umagc, he_ugn, he_rgn
        he_min=dmin1(
     &           sqrt((x1(1,2)-x1(1,1))**2+(x1(2,2)-x1(2,1))**2),
     $           sqrt((x1(1,3)-x1(1,1))**2+(x1(2,3)-x1(2,1))**2),
     $           sqrt((x1(1,3)-x1(1,2))**2+(x1(2,3)-x1(2,2))**2))
!     $          sqrt((x1(1,4)-x1(1,1))**2+(x1(3,4)-x1(3,1))**2))

                he_ugn = he_min

        he_max=dmax1(
     &            sqrt((x1(1,2)-x1(1,1))**2+(x1(2,2)-x1(2,1))**2),
     $            sqrt((x1(1,3)-x1(1,1))**2+(x1(2,3)-x1(2,1))**2),
     $            sqrt((x1(1,3)-x1(1,2))**2+(x1(2,3)-x1(2,2))**2))

                he_rgn = he_max


       if (steady) then
		nutt = 0.0d0
	else 
		nutt= (rho_np1*nut_np1 - rho_n*nut_n)/dt
	end if
!!.......gradnu_sq = nu,i*nu,i
	gradnu_sq = nutx*nutx+nuty*nuty+nutz*nutz

	sgrad =dsqrt((uy-vx)*(uy-vx)+(wx-uz)*(wx-uz)+(vz-wy)*(vz-wy))

	chi = nut/(dif_mol/rho)              
   
	  fv1 = chi**3.0d0/(chi**3.0d0+cv1**3.0d0)
	  fv2 = 1.0d0-chi/(1.0d0+chi*fv1)
	  sbar= sgrad + nut/r_kappa/r_kappa/distance/distance*fv2
!.......pass values of he_trip, umag_trip, vort_trip, dist_trip (at trip)
!.......ft1 and ft2 functions control the transition
	  umag_trip = 0.0d0   !trip point is located on the wall      

	  udif = dabs(umag-umag_trip)
!	  gt   = 0.0d0
!	  fac1 = 0.0d0
!	  fac2 = 0.0d0
	 if(trip) then
	  gt   = dmin1(0.1d0,udif/vort_trip/he_trip)
	  fac1 = distance*distance + gt*gt*dist_trip*dist_trip
	  fac2 = ct2* vort_trip*vort_trip/udif/udif*fac1     

	  t1 = dexp(-fac2)
	  	if (t1 .lt. 1.0e-20) then
	    		ft1 = 0.0000
	  	else
	    		ft1 = ct1*gt*t1
	  	end if

	  	ft2 = ct3*dexp(-ct4*chi*chi)

	 else 
	  ft1=0.0 ! for fully developed turbulence
	  ft2=0.0 ! for fully developed turbulence
	 end if

	  r = nut/sbar/r_kappa/r_kappa/distance/distance

!	  if(r.ge.100.0d0) r=100.0d0
	  if(r.ge.10.0)r=10.0d0
	  g = r + cw2*(r**6.0d0-r)
	  fw= g*((1.0d0+cw3**6.0d0)/
     $        (g**6.0d0+cw3**6.0d0))**(1.0d0/6.0)

	  p_function = rho*cb1*(1.0d0-ft2)*sgrad
	  r_function = rho*ft1*udif*udif      

!	  q_function = rho*(-cw1*fw !  +
!     $   (cb1/r_kappa/r_kappa*(1.0-ft2)*fv2))*(1.0/distance**2.0)
!	  q_function = rho*(-cw1*fw) !  +
!     $	  (cb1/r_kappa/r_kappa*(1.0-ft2)*fv2))*(1.0/distance**2.0)
	  q_function = rho*(-cw1*fw + (cb1/r_kappa/r_kappa*ft2)
     $       +(cb1/r_kappa/r_kappa*(1.0-ft2)*fv2))*(1.0/distance**2.0)

!	write(6,*)sgrad, r_function
!	  write(6,*) q_function, distance, rho, fw, fv2

!		dif_turb = nut*fv1
!		vis = (dif_mol/rho) + dif_turb

!.......SUPG Stabilisation Term
	  tau_supg = 1.0/dsqrt((2.0d0*umagn/he_ugn)**2
!     $             + (12.0d0*(dif_mol/rho)/he_rgnt/he_rgnt)**2)
     $             + (12.0d0*(dif_mol/rho)/he_rgn/he_rgn)**2)
!     $             + (12.0d0*vis/he_rgn/he_rgn)**2)



!	write(6,*) umagn,he_ugn, he_rgn,dif_mol,rho,tau_supg

!	fix it
!.......Galerkin stuff
!         det  = abs(det)
!         eff0 = wq3(ig)*det
	  eff0 = weis(ig)*det
!	fix it 


	  do na =1,nen !na loop
	    b0 = sh(0,na)
	    b1 = sh(1,na)
	    b2 = sh(2,na)
	    b3 = sh(3,na)
!	    b4 = sh(4,na)
!	    b5 = sh(5,na)
!	    b6 = sh(6,na)
	    if(stokes) then
	      ub =0.0d0
	    else
	      ub = u(1)*b1 + u(2)*b2 + u(3)*b3
	    end if

!.......Gravity terms
!	fix it 
!	    p(nb) = p(nb)-                 
!	write(6,*) iunsteady, dt

	    p_t(1,na) = p_t(1,na)-
     &            (b0 + ub * tau_supg)
     &           *(nutt + (rhou*nutx+nut*rhoux+rhov*nuty+rhovy*nut
     &                  +  rhow*nutz+nut*rhowz) 
     &           - p_function*nut - q_function*nut*nut
     &           - r_function - rho*cb2/sigma*gradnu_sq)*eff0
     &           - rho*(nut+dif_mol/rho)/sigma*
     &             (b1*nutx + b2*nuty +b3*nutz )*eff0
!       added (stabilization term)
!       same as incompressible  code
!     $           + rho*(nut+dif_mol/rho)/sigma*ub*tau_supg*
!     $             (nutxx + nutyy+nutzz)*eff0
!       added

	    w_t(1,na)=w_t(1,na)
     &           + (b0 + ub * tau_supg)
     &           * (rho*b0*iunsteady/dt
     &           + ((rhou*b1+rhoux*b0+rhov*b2+rhovy*b0
     &              +rhow*b3+rhowz*b0)
     &           - p_function*b0 - 2.0*q_function*nut*b0
     &           - 2.0*rho*cb2/sigma*
     &             (b1*nutx+b2*nuty+b3*nutz))*alpha)*eff0
     &           +(dif_mol)/sigma*(b1*b1+b2*b2+b3*b3)*alpha*eff0
     &           + rho*b0/sigma *(b1*nutx+b2*nuty+b3*nutz)*alpha*eff0
 !       added(stabilization terms)
!       stabilisation term
!     $           - (dif_mol)/sigma*(b4+b5+b6)*ub*tau_supg*alpha*eff0
!       stabilisation term
!       second order terms
!     $           - rho*b0/sigma *(nutxx+nutyy+nutzz)*
!     $			ub*tau_supg*alpha*eff0
!       second order terms
!       added

c       write(80,180)p_t(1,na),w_t(1,na)
c       write(6,180)p_t(1,na),w_t(1,na)
c       stop
180             format(2e25.15)


	  enddo !na loop
	enddo ! quadrature loop

C...ASSEMBLE INTO GLOBAL VECTOR (STEP 1)
		do inl = 1,nen
			inode = ien(inl,ie)
			do idf = 1,ndf_t
			pgt(idf,inode) = pgt(idf,inode) + p_t(idf,inl)
			wgt(idf,inode) = wgt(idf,inode) + w_t(idf,inl)
			enddo
		enddo


        enddo ! element loop

!        close(81)
        return
        end

