        subroutine blkrhs_turb_sd(ien,xg1,dg1,dg0,pgt,wgt,
     &                    dg1_t,dg0_t,ynode,
!     $                   umag_trans,
!     &                    he_trans,dist_trans,vort_trans,
!     &                    inode_trans,
     $                    ntrans)
c**************************************************************************
c	BLOCK FOR 3D INCOMPRESSIBLE SEMI-DISCRETE GAL/LS METHOD
c       TO COMPUTE RESIDUALS FOR MATRIX FREE GMRES TECHNIQUE. 
c	LOOPS HAVE BEEN UNROLLED TO ACHIEVE SPEED UP
c	WRITTEN BY SANJAY MITTAL AUGUST 31, 2004
c****************************************************************************
	implicit none
	include "global.h"
	include "turbpara.h"
C...EXTERNAL VARIABLES
	integer ien(nen,nec)
	real*8  xg1(nsd,nnl),dg1(ndf,nnl) 
        real*8  wgt(ndf_t,nnl)
	real*8  dg0(ndf,nnl)
        real*8  dg1_t(ndf_t,nnl),dg0_t(ndf_t,nnl),deps_t(ndf_t,nnl)
	real*8  pgt(ndf_t,nnl)
        integer  inode_trans(nnl),ntrans
        real*8  umag_trans(nnl),vort_trans(nnl)
        real*8  he_trans(nnl),dist_trans(nnl)
        real*8  ynode(nnl)

C...INTERNAL VARIABLES
	real*8 do(ndf,nen)
	real*8 x1(nen,nen),d1(ndf,nen)

	real*8 lturb(nen)
	real*8 vis

	real*8 sh(0:9,nen),det

        real*8 do_t(1,nen),d1_t(1,nen),dell_t(1,nen)
        real*8 p_t(1,nen) !,w_t(1,nen) 
        real*8 dist_t(nen),vort_t(nen),he_t(nen),umag_t(nen)
	real*8 yt(nen)
        real*8 distance
	real*8 vort_trip,he_trip,umag_trip,dist_trip

        real*8 nutt,nutx,nuty,nutz,nut
        real*8 nut_n,nutx_n,nuty_n,nutz_n
        real*8 nut_np1,nutx_np1,nuty_np1,nutz_np1
	real*8 grad_nut,gradnu_sq

!	added
	real*8 delnut_np1,delnut,delnutub
	real*8 delnutx_np1, delnutx,delnuty_np1,delnuty
	real*8 delnutz_np1, delnutz
	real*8 delnutxx_np1, delnutxx,delnutyy_np1,delnutyy
	real*8 delnutzz_np1, delnutzz
!	added

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

	real*8 tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8
	real*8 sgrad,fv1,fv2,sbar,umag_trip,udif

	real*8 gt,fac1,fac2,t1,ft1,ft2,r,g,fw

	real*8 p_function,q_function,r_function
	real*8 tau_supg
	real*8 eff0
	real*8 b0,b1,b2,b3,ub
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
        dif = viscosity       ! free stream dynamic viscosity (mu)
!        dif_mol = dif

	Rcon = 1.0d0/gamma
	Cv= Rcon/(gamma-1.0d0)


C...CLEAR BUFFER ARRAY
	call fclear(pgt,nnl*ndf_t)
!	call fclear(wgt,nnl*ndf_t)

!	turb_kappa = 0.15
!	dif = 1.0/reynolds

c...3D TERMS

	do ie = 1,nec

C...CLEAR LOCAL ARRAYS
		do inl = 1,nen
			do idf = 1,ndf_t
				p_t(idf,inl) = 0.0
				w_t(idf,inl) = 0.0
!			enddo
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
				dell_t(idf,inl)= deps_t(idf,k)
                       		do_t(idf,inl) = dg0_t(idf,k)
				d1_t(idf,inl) = dg1_t(idf,k)
                       		enddo

C...........LOCALIZE DISTANCE FROM WALL
			yt(inl) = ynode(k)

			
C...........LOCALIZE TRANSITION TERMS 
                        dist_t(inl)= dist_trans(k)
                        vort_t(inl)= vort_trans(k)
                          he_t(inl)= he_trans(k)
                        umag_t(inl)= umag_trans(k) 

                enddo

		
		do ig=1,nquad3
		
			if(nen.eq.8) then
				include "sf3d8n.h"
			elseif(nen.eq.6)then
				include "sf3d6n.h"
			else
				include "sf3d4n.h"
			endif

!	Distance of integration point from the wall

	  distance = sh(0,1)*yt(1) + sh(0,2)*yt(2) + sh(0,3)*yt(3)
     $             + sh(0,4)*yt(4) + sh(0,5)*yt(5) + sh(0,6)*yt(6)
     $             + sh(0,7)*yt(7) + sh(0,8)*yt(8) + sh(0,9)*yt(9)
     $             + sh(0,10)*yt(10) + sh(0,11)*yt(11) + sh(0,12)*yt(12)
     $             + sh(0,13)*yt(13) + sh(0,14)*yt(14) + sh(0,15)*yt(15)
     $             + sh(0,16)*yt(16) + sh(0,17)*yt(17) + sh(0,18)*yt(18)

!.......Distance of integration point from trip point
	  dist_trip = sh(0,1)*dist_t(1) + sh(0,2)*dist_t(2)
     &              + sh(0,3)*dist_t(3) + sh(0,4)*dist_t(4) 
     &              + sh(0,5)*dist_t(5) + sh(0,6)*dist_t(6)
     &              + sh(0,7)*dist_t(7) + sh(0,8)*dist_t(8) 
     &              + sh(0,9)*dist_t(9) + sh(0,10)*dist_t(10) 
     &              + sh(0,11)*dist_t(11) + sh(0,12)*dist_t(12)
     &              + sh(0,13)*dist_t(13) + sh(0,14)*dist_t(14) 
     &              + sh(0,15)*dist_t(15) + sh(0,16)*dist_t(16)
     &              + sh(0,17)*dist_t(17) + sh(0,18)*dist_t(18) 

!.......Vorticity at trip
	  vort_trip = sh(0,1)*vort_t(1) + sh(0,2)*vort_t(2)
     &              + sh(0,3)*vort_t(3) + sh(0,4)*vort_t(4)
     &              + sh(0,5)*vort_t(5) + sh(0,6)*vort_t(6)
     &              + sh(0,7)*vort_t(7) + sh(0,8)*vort_t(8)
     &              + sh(0,9)*vort_t(9) + sh(0,10)*vort_t(10)
     &              + sh(0,11)*vort_t(11) + sh(0,12)*vort_t(12)
     &              + sh(0,13)*vort_t(13) + sh(0,14)*vort_t(14)
     &              + sh(0,15)*vort_t(15) + sh(0,16)*vort_t(16)
     &              + sh(0,17)*vort_t(17) + sh(0,18)*vort_t(18)

!.......Velocity magnitude at trip
	  umag_trip = sh(0,1)*umag_t(1) + sh(0,2)*umag_t(2)
     &              + sh(0,3)*umag_t(3) + sh(0,4)*umag_t(4)
     &              + sh(0,5)*umag_t(5) + sh(0,6)*umag_t(6)
     &              + sh(0,7)*umag_t(7) + sh(0,8)*umag_t(8)
     &              + sh(0,9)*umag_t(9) + sh(0,10)*umag_t(10)
     &              + sh(0,11)*umag_t(11) + sh(0,12)*umag_t(12)
     &              + sh(0,13)*umag_t(13) + sh(0,14)*umag_t(14)
     &              + sh(0,15)*umag_t(15) + sh(0,16)*umag_t(16)
     &              + sh(0,17)*umag_t(17) + sh(0,18)*umag_t(18)


!.......Element length (along the wall) at trip
	  he_trip = sh(0,1)*he_t(1) + sh(0,2)*he_t(2)
     &            + sh(0,3)*he_t(3) + sh(0,4)*he_t(4)
     &            + sh(0,5)*he_t(5) + sh(0,6)*he_t(6)
     &            + sh(0,7)*he_t(7) + sh(0,8)*he_t(8)
     &            + sh(0,9)*he_t(9) + sh(0,10)*he_t(10)
     &            + sh(0,11)*he_t(11) + sh(0,12)*he_t(12)
     &            + sh(0,13)*he_t(13) + sh(0,14)*he_t(14)
     &            + sh(0,15)*he_t(15) + sh(0,16)*he_t(16)
     &            + sh(0,17)*he_t(18) + sh(0,18)*he_t(18)

	  nut_np1 =  sh(0,1)*d1_t(1,1) + sh(0,2)*d1_t(1,2)
     &             + sh(0,3)*d1_t(1,3) + sh(0,4)*d1_t(1,4)
     &             + sh(0,5)*d1_t(1,5) + sh(0,6)*d1_t(1,6)
     &             + sh(0,7)*d1_t(1,7) + sh(0,8)*d1_t(1,8)
     &             + sh(0,9)*d1_t(1,9) + sh(0,10)*d1_t(1,10)
     &             + sh(0,11)*d1_t(1,11) + sh(0,12)*d1_t(1,12)
     &             + sh(0,13)*d1_t(1,13) + sh(0,14)*d1_t(1,14)
     &             + sh(0,15)*d1_t(1,15) + sh(0,16)*d1_t(1,16)
     &             + sh(0,17)*d1_t(1,17) + sh(0,18)*d1_t(1,18)

	  nut_n =    sh(0,1)*do_t(1,1) + sh(0,2)*do_t(1,2)
     &             + sh(0,3)*do_t(1,3) + sh(0,4)*do_t(1,4)
     &             + sh(0,5)*do_t(1,5) + sh(0,6)*do_t(1,6)
     &             + sh(0,7)*do_t(1,7) + sh(0,8)*do_t(1,8)
     &             + sh(0,9)*do_t(1,9) + sh(0,10)*do_t(1,10)
     &             + sh(0,11)*do_t(1,11) + sh(0,12)*do_t(1,12)
     &             + sh(0,13)*do_t(1,13) + sh(0,14)*do_t(1,14)
     &             + sh(0,15)*do_t(1,15) + sh(0,16)*do_t(1,16)
     &             + sh(0,17)*do_t(1,17) + sh(0,18)*do_t(1,18)
	nut        = alpha*nut_np1 + oma*nut_n

	delnut_np1 =  sh(0,1)*dell_t(1,1) + sh(0,2)*dell_t(1,2)
     &              + sh(0,3)*dell_t(1,3) + sh(0,4)*dell_t(1,4)
     &              + sh(0,5)*dell_t(1,5) + sh(0,6)*dell_t(1,6)
     &              + sh(0,7)*dell_t(1,7) + sh(0,8)*dell_t(1,8)
     &              + sh(0,9)*dell_t(1,9) + sh(0,10)*dell_t(1,10)
     &              + sh(0,11)*dell_t(1,11) + sh(0,12)*dell_t(1,12)
     &              + sh(0,13)*dell_t(1,13) + sh(0,14)*dell_t(1,14)
     &              + sh(0,15)*dell_t(1,15) + sh(0,16)*dell_t(1,16)
     &              + sh(0,17)*dell_t(1,17) + sh(0,18)*dell_t(1,18)
	delnut= alpha * delnut_np1

	  nutx_np1 = sh(1,1)*d1_t(1,1) + sh(1,2)*d1_t(1,2)
     &             + sh(1,3)*d1_t(1,3) + sh(1,4)*d1_t(1,4)
     &             + sh(1,5)*d1_t(1,5) + sh(1,6)*d1_t(1,6)
     &             + sh(1,7)*d1_t(1,7) + sh(1,8)*d1_t(1,8)
     &             + sh(1,9)*d1_t(1,9) + sh(1,10)*d1_t(1,10)
     &             + sh(1,11)*d1_t(1,11) + sh(1,12)*d1_t(1,12)
     &             + sh(1,13)*d1_t(1,13) + sh(1,14)*d1_t(1,14)
     &             + sh(1,15)*d1_t(1,15) + sh(1,16)*d1_t(1,16)
     &             + sh(1,17)*d1_t(1,17) + sh(1,18)*d1_t(1,18)


	  nuty_np1 = sh(2,1)*d1_t(1,1) + sh(2,2)*d1_t(1,2)
     &             + sh(2,3)*d1_t(1,3) + sh(2,4)*d1_t(1,4)
     &             + sh(2,5)*d1_t(1,5) + sh(2,6)*d1_t(1,6)
     &             + sh(2,7)*d1_t(1,7) + sh(2,8)*d1_t(1,8)
     &             + sh(2,9)*d1_t(1,9) + sh(2,10)*d1_t(1,10)
     &             + sh(2,11)*d1_t(1,11) + sh(2,12)*d1_t(1,12)
     &             + sh(2,13)*d1_t(1,13) + sh(2,14)*d1_t(1,14)
     &             + sh(2,15)*d1_t(1,15) + sh(2,16)*d1_t(1,16)
     &             + sh(2,17)*d1_t(1,17) + sh(2,18)*d1_t(1,18)


	  nutz_np1 = sh(3,1)*d1_t(1,1) + sh(3,2)*d1_t(1,2)
     &             + sh(3,3)*d1_t(1,3) + sh(3,4)*d1_t(1,4)
     &             + sh(3,5)*d1_t(1,5) + sh(3,6)*d1_t(1,6)
     &             + sh(3,7)*d1_t(1,7) + sh(3,8)*d1_t(1,8)
     &             + sh(3,9)*d1_t(1,9) + sh(3,10)*d1_t(1,10)
     &             + sh(3,11)*d1_t(1,11) + sh(3,12)*d1_t(1,12)
     &             + sh(3,13)*d1_t(1,13) + sh(3,14)*d1_t(1,14)
     &             + sh(3,15)*d1_t(1,15) + sh(3,16)*d1_t(1,16)
     &             + sh(3,17)*d1_t(1,17) + sh(3,18)*d1_t(1,18)

	delnutx_np1= sh(1,1)*dell_t(1,1) + sh(1,2)*dell_t(1,2)
     &             + sh(1,3)*dell_t(1,3) + sh(1,4)*dell_t(1,4)
     &             + sh(1,5)*dell_t(1,5) + sh(1,6)*dell_t(1,6)
     &             + sh(1,7)*dell_t(1,7) + sh(1,8)*dell_t(1,8)
     &             + sh(1,9)*dell_t(1,9) + sh(1,10)*dell_t(1,10)
     &             + sh(1,11)*dell_t(1,11) + sh(1,12)*dell_t(1,12)
     &             + sh(1,13)*dell_t(1,13) + sh(1,14)*dell_t(1,14)
     &             + sh(1,15)*dell_t(1,15) + sh(1,16)*dell_t(1,16)
     &             + sh(1,17)*dell_t(1,17) + sh(1,18)*dell_t(1,18)


	delnuty_np1= sh(2,1)*dell_t(1,1) + sh(2,2)*dell_t(1,2)
     &             + sh(2,3)*dell_t(1,3) + sh(2,4)*dell_t(1,4)
     &             + sh(2,5)*dell_t(1,5) + sh(2,6)*dell_t(1,6)
     &             + sh(2,7)*dell_t(1,7) + sh(2,8)*dell_t(1,8)
     &             + sh(2,9)*dell_t(1,9) + sh(2,10)*dell_t(1,10)
     &             + sh(2,11)*dell_t(1,11) + sh(2,12)*dell_t(1,12)
     &             + sh(2,13)*dell_t(1,13) + sh(2,14)*dell_t(1,14)
     &             + sh(2,15)*dell_t(1,15) + sh(2,16)*dell_t(1,16)
     &             + sh(2,17)*dell_t(1,17) + sh(2,18)*dell_t(1,18)


	delnutz_np1= sh(3,1)*dell_t(1,1) + sh(3,2)*dell_t(1,2)
     &             + sh(3,3)*dell_t(1,3) + sh(3,4)*dell_t(1,4)
     &             + sh(3,5)*dell_t(1,5) + sh(3,6)*dell_t(1,6)
     &             + sh(3,7)*dell_t(1,7) + sh(3,8)*dell_t(1,8)
     &             + sh(3,9)*dell_t(1,9) + sh(3,10)*dell_t(1,10)
     &             + sh(3,11)*dell_t(1,11) + sh(3,12)*dell_t(1,12)
     &             + sh(3,13)*dell_t(1,13) + sh(3,14)*dell_t(1,14)
     &             + sh(3,15)*dell_t(1,15) + sh(3,16)*dell_t(1,16)
     &             + sh(3,17)*dell_t(1,17) + sh(3,18)*dell_t(1,18)

	delnutxx_np1= sh(4,1)*dell_t(1,1)+ sh(4,2)*dell_t(1,2)
     &             + sh(4,3)*dell_t(1,3) + sh(4,4)*dell_t(1,4)
     &             + sh(4,5)*dell_t(1,5) + sh(4,6)*dell_t(1,6)
     &             + sh(4,7)*dell_t(1,7) + sh(4,8)*dell_t(1,8)
     &             + sh(4,9)*dell_t(1,9) + sh(4,10)*dell_t(1,10)
     &             + sh(4,11)*dell_t(1,11) + sh(4,12)*dell_t(1,12)
     &             + sh(4,13)*dell_t(1,13) + sh(4,14)*dell_t(1,14)
     &             + sh(4,15)*dell_t(1,15) + sh(4,16)*dell_t(1,16)
     &             + sh(4,17)*dell_t(1,17) + sh(4,18)*dell_t(1,18)

	delnutyy_np1= sh(5,1)*dell_t(1,1)+ sh(5,2)*dell_t(1,2)
     &             + sh(5,3)*dell_t(1,3) + sh(5,4)*dell_t(1,4)
     &             + sh(5,5)*dell_t(1,5) + sh(5,6)*dell_t(1,6)
     &             + sh(5,7)*dell_t(1,7) + sh(5,8)*dell_t(1,8)
     &             + sh(5,9)*dell_t(1,9) + sh(5,10)*dell_t(1,10)
     &             + sh(5,11)*dell_t(1,11) + sh(5,12)*dell_t(1,12)
     &             + sh(5,13)*dell_t(1,13) + sh(5,14)*dell_t(1,14)
     &             + sh(5,15)*dell_t(1,15) + sh(5,16)*dell_t(1,16)
     &             + sh(5,17)*dell_t(1,17) + sh(5,18)*dell_t(1,18)

	delnutzz_np1= sh(6,1)*dell_t(1,1)+ sh(6,2)*dell_t(1,2)
     &             + sh(6,3)*dell_t(1,3) + sh(6,4)*dell_t(1,4)
     &             + sh(6,5)*dell_t(1,5) + sh(6,6)*dell_t(1,6)
     &             + sh(6,7)*dell_t(1,7) + sh(6,8)*dell_t(1,8)
     &             + sh(6,9)*dell_t(1,9) + sh(6,10)*dell_t(1,10)
     &             + sh(6,11)*dell_t(1,11) + sh(6,12)*dell_t(1,12)
     &             + sh(6,13)*dell_t(1,13) + sh(6,14)*dell_t(1,14)
     &             + sh(6,15)*dell_t(1,15) + sh(6,16)*dell_t(1,16)
     &             + sh(6,17)*dell_t(1,17) + sh(6,18)*dell_t(1,18)

	delnutx = alpha * delnutx_np1
	delnuty = alpha * delnuty_np1
	delnutz = alpha * delnutz_np1
	delnutxx = alpha * delnutxx_np1
	delnutyy = alpha * delnutyy_np1
	delnutzz = alpha * delnutzz_np1


	  nutx_n   = sh(1,1)*do_t(1,1) + sh(1,2)*do_t(1,2)
     &             + sh(1,3)*do_t(1,3) + sh(1,4)*do_t(1,4)
     &             + sh(1,5)*do_t(1,5) + sh(1,6)*do_t(1,6)
     &             + sh(1,7)*do_t(1,7) + sh(1,8)*do_t(1,8)
     &             + sh(1,9)*do_t(1,9) + sh(1,10)*do_t(1,10)
     &             + sh(1,11)*do_t(1,11) + sh(1,12)*do_t(1,12)
     &             + sh(1,13)*do_t(1,13) + sh(1,14)*do_t(1,14)
     &             + sh(1,15)*do_t(1,15) + sh(1,16)*do_t(1,16)
     &             + sh(1,17)*do_t(1,17) + sh(1,18)*do_t(1,18)


	  nuty_n   = sh(2,1)*do_t(1,1)  + sh(2,2)*do_t(1,2)
     &             + sh(2,3)*do_t(1,3)  + sh(2,4)*do_t(1,4)
     &             + sh(2,5)*do_t(1,5)  + sh(2,6)*do_t(1,6)
     &             + sh(2,7)*do_t(1,7)  + sh(2,8)*do_t(1,8)
     &             + sh(2,9)*do_t(1,9)  + sh(2,10)*do_t(1,10)
     &             + sh(2,11)*do_t(1,11)  + sh(2,12)*do_t(1,12)
     &             + sh(2,13)*do_t(1,13)  + sh(2,14)*do_t(1,14)
     &             + sh(2,15)*do_t(1,15)  + sh(2,16)*do_t(1,16)
     &             + sh(2,17)*do_t(1,17)  + sh(2,18)*do_t(1,18)


	  nutz_n   = sh(3,1)*do_t(1,1)  + sh(3,2)*do_t(1,2)
     &             + sh(3,3)*do_t(1,3)  + sh(3,4)*do_t(1,4)
     &             + sh(3,5)*do_t(1,5)  + sh(3,6)*do_t(1,6)
     &             + sh(3,7)*do_t(1,7)  + sh(3,8)*do_t(1,8)
     &             + sh(3,9)*do_t(1,9)  + sh(3,10)*do_t(1,10)
     &             + sh(3,11)*do_t(1,11)  + sh(3,12)*do_t(1,12)
     &             + sh(3,13)*do_t(1,13)  + sh(3,14)*do_t(1,14)
     &             + sh(3,15)*do_t(1,15)  + sh(3,16)*do_t(1,16)
     &             + sh(3,17)*do_t(1,17)  + sh(3,18)*do_t(1,18)



	  nutx  = alpha * nutx_np1 + oma * nutx_n 
	  nuty  = alpha * nuty_np1 + oma * nuty_n
	  nutz  = alpha * nutz_np1 + oma * nutz_n

	  grad_nut = dsqrt(nutx_n*nutx_n+nuty_n*nuty_n+nutz_n*nutz_n)


	  rho_np1 = sh(0,1)*d1(1,1) + sh(0,2)*d1(1,2)
     &            + sh(0,3)*d1(1,3) + sh(0,4)*d1(1,4)
     &            + sh(0,5)*d1(1,5) + sh(0,6)*d1(1,6)
     &            + sh(0,7)*d1(1,7) + sh(0,8)*d1(1,8)
     &            + sh(0,9)*d1(1,9) + sh(0,10)*d1(1,10)
     &            + sh(0,11)*d1(1,11) + sh(0,12)*d1(1,12)
     &            + sh(0,13)*d1(1,13) + sh(0,14)*d1(1,14)
     &            + sh(0,15)*d1(1,15) + sh(0,16)*d1(1,16)
     &            + sh(0,17)*d1(1,17) + sh(0,18)*d1(1,18)


	  rhou_np1= sh(0,1)*d1(2,1) + sh(0,2)*d1(2,2)
     &            + sh(0,3)*d1(2,3) + sh(0,4)*d1(2,4)
     &            + sh(0,5)*d1(2,5) + sh(0,6)*d1(2,6)
     &            + sh(0,7)*d1(2,7) + sh(0,8)*d1(2,8)
     &            + sh(0,9)*d1(2,9) + sh(0,10)*d1(2,10)
     &            + sh(0,11)*d1(2,11) + sh(0,12)*d1(2,12)
     &            + sh(0,13)*d1(2,13) + sh(0,14)*d1(2,14)
     &            + sh(0,15)*d1(2,15) + sh(0,16)*d1(2,16)
     &            + sh(0,17)*d1(2,17) + sh(0,18)*d1(2,18)


	  rhov_np1= sh(0,1)*d1(3,1) + sh(0,2)*d1(3,2)
     &            + sh(0,3)*d1(3,3) + sh(0,4)*d1(3,4)
     &            + sh(0,5)*d1(3,5) + sh(0,6)*d1(3,6)
     &            + sh(0,7)*d1(3,7) + sh(0,8)*d1(3,8)
     &            + sh(0,9)*d1(3,9) + sh(0,10)*d1(3,10)
     &            + sh(0,11)*d1(3,11) + sh(0,12)*d1(3,12)
     &            + sh(0,13)*d1(3,13) + sh(0,14)*d1(3,14)
     &            + sh(0,15)*d1(3,15) + sh(0,16)*d1(3,16)
     &            + sh(0,17)*d1(3,17) + sh(0,18)*d1(3,18)


	  rhow_np1= sh(0,1)*d1(4,1) + sh(0,2)*d1(4,2)
     &            + sh(0,3)*d1(4,3) + sh(0,4)*d1(4,4)
     &            + sh(0,5)*d1(4,5) + sh(0,6)*d1(4,6)
     &            + sh(0,7)*d1(4,7) + sh(0,8)*d1(4,8)
     &            + sh(0,9)*d1(4,9) + sh(0,10)*d1(4,10)
     &            + sh(0,11)*d1(4,11) + sh(0,12)*d1(4,12)
     &            + sh(0,13)*d1(4,13) + sh(0,14)*d1(4,14)
     &            + sh(0,15)*d1(4,15) + sh(0,16)*d1(4,16)
     &            + sh(0,17)*d1(4,17) + sh(0,18)*d1(4,18)


	  rhoe_np1= sh(0,1)*d1(5,1) + sh(0,2)*d1(5,2)
     &            + sh(0,3)*d1(5,3) + sh(0,4)*d1(5,4)
     &            + sh(0,5)*d1(5,5) + sh(0,6)*d1(5,6)
     &            + sh(0,7)*d1(5,7) + sh(0,8)*d1(5,8)
     &            + sh(0,9)*d1(5,9) + sh(0,10)*d1(5,10)
     &            + sh(0,11)*d1(5,11) + sh(0,12)*d1(5,12)
     &            + sh(0,13)*d1(5,13) + sh(0,14)*d1(5,14)
     &            + sh(0,15)*d1(5,15) + sh(0,16)*d1(5,16)
     &            + sh(0,17)*d1(5,17) + sh(0,18)*d1(5,18)


	  rho_n   = sh(0,1)*do(1,1) + sh(0,2)*do(1,2)
     &            + sh(0,3)*do(1,3) + sh(0,4)*do(1,4)
     &            + sh(0,5)*do(1,5) + sh(0,6)*do(1,6)
     &            + sh(0,7)*do(1,7) + sh(0,8)*do(1,8)
     &            + sh(0,9)*do(1,9) + sh(0,10)*do(1,10)
     &            + sh(0,11)*do(1,11) + sh(0,12)*do(1,12)
     &            + sh(0,13)*do(1,13) + sh(0,14)*do(1,14)
     &            + sh(0,15)*do(1,15) + sh(0,16)*do(1,16)
     &            + sh(0,17)*do(1,17) + sh(0,18)*do(1,18)



	  rhou_n  = sh(0,1)*do(2,1) + sh(0,2)*do(2,2)
     &            + sh(0,3)*do(2,3) + sh(0,4)*do(2,4)
     &            + sh(0,5)*do(2,5) + sh(0,6)*do(2,6)
     &            + sh(0,7)*do(2,7) + sh(0,8)*do(2,8)
     &            + sh(0,9)*do(2,9) + sh(0,10)*do(2,10)
     &            + sh(0,11)*do(2,11) + sh(0,12)*do(2,12)
     &            + sh(0,13)*do(2,13) + sh(0,14)*do(2,14)
     &            + sh(0,15)*do(2,15) + sh(0,16)*do(2,16)
     &            + sh(0,17)*do(2,17) + sh(0,18)*do(2,18)

	  rhov_n  = sh(0,1)*do(3,1) + sh(0,2)*do(3,2)
     &            + sh(0,3)*do(3,3) + sh(0,4)*do(3,4)
     &            + sh(0,5)*do(3,5) + sh(0,6)*do(3,6)
     &            + sh(0,7)*do(3,7) + sh(0,8)*do(3,8)
     &            + sh(0,9)*do(3,9) + sh(0,10)*do(3,10)
     &            + sh(0,11)*do(3,11) + sh(0,12)*do(3,12)
     &            + sh(0,13)*do(3,13) + sh(0,14)*do(3,14)
     &            + sh(0,15)*do(3,15) + sh(0,16)*do(3,16)
     &            + sh(0,17)*do(3,17) + sh(0,18)*do(3,18)

	  rhow_n  = sh(0,1)*do(4,1) + sh(0,2)*do(4,2)
     &            + sh(0,3)*do(4,3) + sh(0,4)*do(4,4)
     &            + sh(0,5)*do(4,5) + sh(0,6)*do(4,6)
     &            + sh(0,7)*do(4,7) + sh(0,8)*do(4,8)
     &            + sh(0,9)*do(4,9) + sh(0,10)*do(4,10)
     &            + sh(0,11)*do(4,11) + sh(0,12)*do(4,12)
     &            + sh(0,13)*do(4,13) + sh(0,14)*do(4,14)
     &            + sh(0,15)*do(4,15) + sh(0,16)*do(4,16)
     &            + sh(0,17)*do(4,17) + sh(0,18)*do(4,18)


	  rhoe_n  = sh(0,1)*do(5,1) + sh(0,2)*do(5,2)
     &            + sh(0,3)*do(5,3) + sh(0,4)*do(5,4)
     &		  + sh(0,5)*do(5,5) + sh(0,6)*do(5,6)
     &            + sh(0,7)*do(5,7) + sh(0,8)*do(5,8)
     &            + sh(0,9)*do(5,9) + sh(0,10)*do(5,10)
     &            + sh(0,11)*do(5,11) + sh(0,12)*do(5,12)
     &            + sh(0,13)*do(5,13) + sh(0,14)*do(5,14)
     &            + sh(0,15)*do(5,15) + sh(0,16)*do(5,16)
     &            + sh(0,17)*do(5,17) + sh(0,18)*do(5,18)

	  unp1(1) = rhou_np1/rho_np1
	  unp1(2) = rhov_np1/rho_np1
	  unp1(3) = rhow_np1/rho_np1
	  un(1)   = rhou_n/rho_n
	  un(2)   = rhov_n/rho_n
	  un(3)   = rhow_n/rho_n

	  u(1) = alpha * unp1(1) + oma * un(1)
	  u(2) = alpha * unp1(2) + oma * un(2)
	  u(3) = alpha * unp1(3) + oma * un(3)

	  rho  = alpha * rho_np1 + oma * rho_n
	  rhou = alpha * rhou_np1+ oma * rhou_n
	  rhov = alpha * rhov_np1+ oma * rhov_n
	  rhow = alpha * rhow_np1+ oma * rhow_n
	  rhoe = alpha * rhoe_np1+ oma * rhoe_n

	  umag  = dsqrt(u(1)*u(1) +  u(2)*u(2)+u(3)*u(3))        !used in he_ugn
	  umagn = dsqrt(un(1)*un(1) + un(2)*un(2)+un(3)*un(3))   !used in tau_supg
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
     &             + sh(1,7)*d1(2,7) + sh(1,8)*d1(2,8) 
     &             + sh(1,9)*d1(2,9) + sh(1,10)*d1(2,10)
     &             + sh(1,11)*d1(2,11) + sh(1,12)*d1(2,12) 
     &             + sh(1,13)*d1(2,13) + sh(1,14)*d1(2,14)
     &             + sh(1,15)*d1(2,15) + sh(1,16)*d1(2,16) 
     &             + sh(1,17)*d1(2,17) + sh(1,18)*d1(2,18)
 
!	  rhouy_np1= sh(2,1)*d1(2,1) + sh(2,2)*d1(2,2)
!     &             + sh(2,3)*d1(2,3) 
!	  rhovx_np1= sh(1,1)*d1(3,1) + sh(1,2)*d1(3,2)
!     &             + sh(1,3)*d1(3,3) 
	  rhouy_np1= sh(2,1)*d1(2,1) + sh(2,2)*d1(2,2)
     &             + sh(2,3)*d1(2,3) + sh(2,4)*d1(2,4) 
     &             + sh(2,5)*d1(2,5) + sh(2,6)*d1(2,6)
     &             + sh(2,7)*d1(2,7) + sh(2,8)*d1(2,8) 
     &             + sh(2,9)*d1(2,9) + sh(2,10)*d1(2,10)
     &             + sh(2,11)*d1(2,11) + sh(2,12)*d1(2,12) 
     &             + sh(2,13)*d1(2,13) + sh(2,14)*d1(2,14)
     &             + sh(2,15)*d1(2,15) + sh(2,16)*d1(2,16) 
     &             + sh(2,17)*d1(2,17) + sh(2,18)*d1(2,18)
 

	  rhouz_np1= sh(3,1)*d1(2,1) + sh(3,2)*d1(2,2)
     &             + sh(3,3)*d1(2,3) + sh(3,4)*d1(2,4) 
     &             + sh(3,5)*d1(2,5) + sh(3,6)*d1(2,6)
     &             + sh(3,7)*d1(2,7) + sh(3,8)*d1(2,8) 
     &             + sh(3,9)*d1(2,9) + sh(3,10)*d1(2,10)
     &             + sh(3,11)*d1(2,11) + sh(3,12)*d1(2,12) 
     &             + sh(3,13)*d1(2,13) + sh(3,14)*d1(2,14)
     &             + sh(3,15)*d1(2,15) + sh(3,16)*d1(2,16) 
     &             + sh(3,17)*d1(2,17) + sh(3,18)*d1(2,18)

	  rhovx_np1= sh(1,1)*d1(3,1) + sh(1,2)*d1(3,2)
     &             + sh(1,3)*d1(3,3) + sh(1,4)*d1(3,4)
     &             + sh(1,5)*d1(3,5) + sh(1,6)*d1(3,6)
     &             + sh(1,7)*d1(3,7) + sh(1,8)*d1(3,8) 
     &             + sh(1,9)*d1(3,9) + sh(1,10)*d1(3,10)
     &             + sh(1,11)*d1(3,11) + sh(1,12)*d1(3,12) 
     &             + sh(1,13)*d1(3,13) + sh(1,14)*d1(3,14)
     &             + sh(1,15)*d1(3,15) + sh(1,16)*d1(3,16) 
     &             + sh(1,17)*d1(3,17) + sh(1,18)*d1(3,18)

	  rhovy_np1= sh(2,1)*d1(3,1) + sh(2,2)*d1(3,2)
     &             + sh(2,3)*d1(3,3) + sh(2,4)*d1(3,4)
     &             + sh(2,5)*d1(3,5) + sh(2,6)*d1(3,6)
     &             + sh(2,7)*d1(3,7) + sh(2,8)*d1(3,8) 
     &             + sh(2,9)*d1(3,9) + sh(2,10)*d1(3,10)
     &             + sh(2,11)*d1(3,11) + sh(2,12)*d1(3,12) 
     &             + sh(2,13)*d1(3,13) + sh(2,14)*d1(3,14)
     &             + sh(2,15)*d1(3,15) + sh(2,16)*d1(3,16) 
     &             + sh(2,17)*d1(3,17) + sh(2,18)*d1(3,18)


	  rhovz_np1= sh(3,1)*d1(3,1) + sh(3,2)*d1(3,2)
     &             + sh(3,3)*d1(3,3) + sh(3,4)*d1(3,4)
     &             + sh(3,5)*d1(3,5) + sh(3,6)*d1(3,6)
     &             + sh(3,7)*d1(3,7) + sh(3,8)*d1(3,8) 
     &             + sh(3,9)*d1(3,9) + sh(3,10)*d1(3,10)
     &             + sh(3,11)*d1(3,11) + sh(3,12)*d1(3,12) 
     &             + sh(3,13)*d1(3,13) + sh(3,14)*d1(3,14)
     &             + sh(3,15)*d1(3,15) + sh(3,16)*d1(3,16) 
     &             + sh(3,17)*d1(3,17) + sh(3,18)*d1(3,18)

	  rhowx_np1= sh(1,1)*d1(4,1) + sh(1,2)*d1(4,2)
     &             + sh(1,3)*d1(4,3) + sh(1,4)*d1(4,4)
     &             + sh(1,5)*d1(4,5) + sh(1,6)*d1(4,6)
     &             + sh(1,7)*d1(4,7) + sh(1,8)*d1(4,8) 
     &             + sh(1,9)*d1(4,9) + sh(1,10)*d1(4,10)
     &             + sh(1,11)*d1(4,11) + sh(1,12)*d1(4,12) 
     &             + sh(1,13)*d1(4,13) + sh(1,14)*d1(4,14)
     &             + sh(1,15)*d1(4,15) + sh(1,16)*d1(4,16) 
     &             + sh(1,17)*d1(4,17) + sh(1,18)*d1(4,18)

	  rhowy_np1= sh(2,1)*d1(4,1) + sh(2,2)*d1(4,2)
     &             + sh(2,3)*d1(4,3) + sh(2,4)*d1(4,4)
     &             + sh(2,5)*d1(4,5) + sh(2,6)*d1(4,6)
     &             + sh(2,7)*d1(4,7) + sh(2,8)*d1(4,8) 
     &             + sh(2,9)*d1(4,9) + sh(2,10)*d1(4,10)
     &             + sh(2,11)*d1(4,11) + sh(2,12)*d1(4,12) 
     &             + sh(2,13)*d1(4,13) + sh(2,14)*d1(4,14)
     &             + sh(2,15)*d1(4,15) + sh(2,16)*d1(4,16) 
     &             + sh(2,17)*d1(4,17) + sh(2,18)*d1(4,18)

	  rhowz_np1= sh(3,1)*d1(4,1) + sh(3,2)*d1(4,2)
     &             + sh(3,3)*d1(4,3) + sh(3,4)*d1(4,4)
     &             + sh(3,5)*d1(4,5) + sh(3,6)*d1(4,6)
     &             + sh(3,7)*d1(4,7) + sh(3,8)*d1(4,8) 
     &             + sh(3,9)*d1(4,9) + sh(3,10)*d1(4,10)
     &             + sh(3,11)*d1(4,11) + sh(3,12)*d1(4,12) 
     &             + sh(3,13)*d1(4,13) + sh(3,14)*d1(4,14)
     &             + sh(3,15)*d1(4,15) + sh(3,16)*d1(4,16) 
     &             + sh(3,17)*d1(4,17) + sh(3,18)*d1(4,18)

	  rhoux_n  = sh(1,1)*do(2,1) + sh(1,2)*do(2,2)
     &             + sh(1,3)*do(2,3) + sh(1,4)*do(2,4) 
     &             + sh(1,5)*do(2,5) + sh(1,6)*do(2,6)
     &             + sh(1,7)*do(2,7) + sh(1,8)*do(2,8) 
     &             + sh(1,9)*do(2,9) + sh(1,10)*do(2,10)
     &             + sh(1,11)*do(2,11) + sh(1,12)*do(2,12) 
     &             + sh(1,13)*do(2,13) + sh(1,14)*do(2,14)
     &             + sh(1,15)*do(2,15) + sh(1,16)*do(2,16) 
     &             + sh(1,17)*do(2,17) + sh(1,18)*do(2,18)


	  rhouy_n  = sh(2,1)*do(2,1) + sh(2,2)*do(2,2)
     &             + sh(2,3)*do(2,3) + sh(2,4)*do(2,4) 
     &             + sh(2,5)*do(2,5) + sh(2,6)*do(2,6)
     &             + sh(2,7)*do(2,7) + sh(2,8)*do(2,8) 
     &             + sh(2,9)*do(2,9) + sh(2,10)*do(2,10)
     &             + sh(2,11)*do(2,11) + sh(2,12)*do(2,12) 
     &             + sh(2,13)*do(2,13) + sh(2,14)*do(2,14)
     &             + sh(2,15)*do(2,15) + sh(2,16)*do(2,16) 
     &             + sh(2,17)*do(2,17) + sh(2,18)*do(2,18)
  
	  rhouz_n  = sh(3,1)*do(2,1) + sh(3,2)*do(2,2)
     &             + sh(3,3)*do(2,3) + sh(3,4)*do(2,4) 
     &             + sh(3,5)*do(2,5) + sh(3,6)*do(2,6) 
     &             + sh(3,7)*do(2,7) + sh(3,8)*do(2,8) 
     &             + sh(3,9)*do(2,9) + sh(3,10)*do(2,10)
     &             + sh(3,11)*do(2,11) + sh(3,12)*do(2,12) 
     &             + sh(3,13)*do(2,13) + sh(3,14)*do(2,14)
     &             + sh(3,15)*do(2,15) + sh(3,16)*do(2,16) 
     &             + sh(3,17)*do(2,17) + sh(3,18)*do(2,18)


	  rhovx_n = sh(1,1)*do(3,1)  + sh(1,2)*do(3,2)
     &             + sh(1,3)*do(3,3) + sh(1,4)*do(3,4)
     &             + sh(1,5)*do(3,5) + sh(1,6)*do(3,6)
     &             + sh(1,7)*do(3,7) + sh(1,8)*do(3,8) 
     &             + sh(1,9)*do(3,9) + sh(1,10)*do(3,10)
     &             + sh(1,11)*do(3,11) + sh(1,12)*do(3,12) 
     &             + sh(1,13)*do(3,13) + sh(1,14)*do(3,14)
     &             + sh(1,15)*do(3,15) + sh(1,16)*do(3,16) 
     &             + sh(1,17)*do(3,17) + sh(1,18)*do(3,18)


	  rhovy_n  = sh(2,1)*do(3,1) + sh(2,2)*do(3,2)
     &             + sh(2,3)*do(3,3) + sh(2,4)*do(3,4)
     &             + sh(2,5)*do(3,5) + sh(2,6)*do(3,6)
     &             + sh(2,7)*do(3,7) + sh(2,8)*do(3,8) 
     &             + sh(2,9)*do(3,9) + sh(2,10)*do(3,10)
     &             + sh(2,11)*do(3,11) + sh(2,12)*do(3,12) 
     &             + sh(2,13)*do(3,13) + sh(2,14)*do(3,14)
     &             + sh(2,15)*do(3,15) + sh(2,16)*do(3,16) 
     &             + sh(2,17)*do(3,17) + sh(2,18)*do(3,18)

	  rhovz_n  = sh(3,1)*do(3,1) + sh(3,2)*do(3,2)
     &             + sh(3,3)*do(3,3) + sh(3,4)*do(3,4)
     &             + sh(3,5)*do(3,5) + sh(3,6)*do(3,6)
     &             + sh(3,7)*do(3,7) + sh(3,8)*do(3,8)
     &             + sh(3,9)*do(3,9) + sh(3,10)*do(3,10)
     &             + sh(3,11)*do(3,11) + sh(3,12)*do(3,12) 
     &             + sh(3,13)*do(3,13) + sh(3,14)*do(3,14)
     &             + sh(3,15)*do(3,15) + sh(3,16)*do(3,16) 
     &             + sh(3,17)*do(3,17) + sh(3,18)*do(3,18)

	  rhowx_n  = sh(1,1)*do(4,1) + sh(1,2)*do(4,2)
     &             + sh(1,3)*do(4,3) + sh(1,4)*do(4,4)
     &             + sh(1,5)*do(4,5) + sh(1,6)*do(4,6)
     &             + sh(1,7)*do(4,7) + sh(1,8)*do(4,8) 
     &             + sh(1,9)*do(4,9) + sh(1,10)*do(4,10)
     &             + sh(1,11)*do(4,11) + sh(1,12)*do(4,12) 
     &             + sh(1,13)*do(4,13) + sh(1,14)*do(4,14)
     &             + sh(1,15)*do(4,15) + sh(1,16)*do(4,16) 
     &             + sh(1,17)*do(4,17) + sh(1,18)*do(4,18)


	  rhowy_n  = sh(2,1)*do(4,1) + sh(2,2)*do(4,2)
     &             + sh(2,3)*do(4,3) + sh(2,4)*do(4,4)
     &             + sh(2,5)*do(4,5) + sh(2,6)*do(4,6)
     &             + sh(2,7)*do(4,7) + sh(2,8)*do(4,8) 
     &             + sh(2,9)*do(4,9) + sh(2,10)*do(4,10)
     &             + sh(2,11)*do(4,11) + sh(2,12)*do(4,12) 
     &             + sh(2,13)*do(4,13) + sh(2,14)*do(4,14)
     &             + sh(2,15)*do(4,15) + sh(2,16)*do(4,16) 
     &             + sh(2,17)*do(4,17) + sh(2,18)*do(4,18)

	  rhowz_n  = sh(3,1)*do(4,1) + sh(3,2)*do(4,2)
     &             + sh(3,3)*do(4,3) + sh(3,4)*do(4,4)
     &             + sh(3,5)*do(4,5) + sh(3,6)*do(4,6)
     &             + sh(3,7)*do(4,7) + sh(3,8)*do(4,8) 
     &             + sh(3,9)*do(4,9) + sh(3,10)*do(4,10)
     &             + sh(3,11)*do(4,11) + sh(3,12)*do(4,12)
     &             + sh(3,13)*do(4,13) + sh(3,14)*do(4,14)
     &             + sh(3,15)*do(4,15) + sh(3,16)*do(4,16)
     &             + sh(3,17)*do(4,17) + sh(3,18)*do(4,18)

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
     &          +sh(2,7)*d1(2,7)/d1(1,7)+sh(2,8)*d1(2,8)/d1(1,8)
     &          +sh(2,9)*d1(2,9)/d1(1,9)+sh(2,10)*d1(2,10)/d1(1,10)
     &          +sh(2,11)*d1(2,11)/d1(1,11)+sh(2,12)*d1(2,12)/d1(1,12)
     &          +sh(2,13)*d1(2,13)/d1(1,13)+sh(2,14)*d1(2,14)/d1(1,14)
     &          +sh(2,15)*d1(2,15)/d1(1,15)+sh(2,16)*d1(2,16)/d1(1,16)
     &          +sh(2,17)*d1(2,17)/d1(1,17)+sh(2,18)*d1(2,18)/d1(1,18)

	  uz_np1=sh(3,1)*d1(2,1)/d1(1,1)+sh(3,2)*d1(2,2)/d1(1,2)
     &          +sh(3,3)*d1(2,3)/d1(1,3)+sh(3,4)*d1(2,4)/d1(1,4)
     &          +sh(3,5)*d1(2,5)/d1(1,5)+sh(3,6)*d1(2,6)/d1(1,6)
     &          +sh(3,7)*d1(2,7)/d1(1,7)+sh(3,8)*d1(2,8)/d1(1,8)
     &          +sh(3,9)*d1(2,9)/d1(1,9)+sh(3,10)*d1(2,10)/d1(1,10)
     &          +sh(3,11)*d1(2,11)/d1(1,11)+sh(3,12)*d1(2,12)/d1(1,12)
     &          +sh(3,13)*d1(2,13)/d1(1,13)+sh(3,14)*d1(2,14)/d1(1,14)
     &          +sh(3,15)*d1(2,15)/d1(1,15)+sh(3,16)*d1(2,16)/d1(1,16)
     &          +sh(3,17)*d1(2,17)/d1(1,17)+sh(3,18)*d1(2,18)/d1(1,18)


	  vx_np1=sh(1,1)*d1(3,1)/d1(1,1)+sh(1,2)*d1(3,2)/d1(1,2)
     &          +sh(1,3)*d1(3,3)/d1(1,3)+sh(1,4)*d1(3,4)/d1(1,4)
     &          +sh(1,5)*d1(3,5)/d1(1,5)+sh(1,6)*d1(3,6)/d1(1,6)
     &          +sh(1,7)*d1(3,7)/d1(1,7)+sh(1,8)*d1(3,8)/d1(1,8)
     &          +sh(1,9)*d1(3,9)/d1(1,9)+sh(1,10)*d1(3,10)/d1(1,10)
     &          +sh(1,11)*d1(3,11)/d1(1,11)+sh(1,12)*d1(3,12)/d1(1,12)
     &          +sh(1,13)*d1(3,13)/d1(1,13)+sh(1,14)*d1(3,14)/d1(1,14)
     &          +sh(1,15)*d1(3,15)/d1(1,15)+sh(1,16)*d1(3,16)/d1(1,16)
     &          +sh(1,17)*d1(3,17)/d1(1,17)+sh(1,18)*d1(3,18)/d1(1,18)


	  vz_np1=sh(3,1)*d1(3,1)/d1(1,1)+sh(3,2)*d1(3,2)/d1(1,2)
     &          +sh(3,3)*d1(3,3)/d1(1,3)+sh(3,4)*d1(3,4)/d1(1,4)
     &          +sh(3,5)*d1(3,5)/d1(1,5)+sh(3,6)*d1(3,6)/d1(1,6)
     &          +sh(3,7)*d1(3,7)/d1(1,7)+sh(3,8)*d1(3,8)/d1(1,8)
     &          +sh(3,9)*d1(3,9)/d1(1,9)+sh(3,10)*d1(3,10)/d1(1,10)
     &          +sh(3,11)*d1(3,11)/d1(1,11)+sh(3,12)*d1(3,12)/d1(1,12)
     &          +sh(3,13)*d1(3,13)/d1(1,13)+sh(3,14)*d1(3,14)/d1(1,14)
     &          +sh(3,15)*d1(3,15)/d1(1,15)+sh(3,16)*d1(3,16)/d1(1,16)
     &          +sh(3,17)*d1(3,17)/d1(1,17)+sh(3,18)*d1(3,18)/d1(1,18)

	  wx_np1=sh(1,1)*d1(4,1)/d1(1,1)+sh(1,2)*d1(4,2)/d1(1,2)
     &          +sh(1,3)*d1(4,3)/d1(1,3)+sh(1,4)*d1(4,4)/d1(1,4)
     &          +sh(1,5)*d1(4,5)/d1(1,5)+sh(1,6)*d1(4,6)/d1(1,6)
     &          +sh(1,7)*d1(4,7)/d1(1,7)+sh(1,8)*d1(4,8)/d1(1,8)
     &          +sh(1,9)*d1(4,9)/d1(1,9)+sh(1,10)*d1(4,10)/d1(1,10)
     &          +sh(1,11)*d1(4,11)/d1(1,11)+sh(1,12)*d1(4,12)/d1(1,12)
     &          +sh(1,13)*d1(4,13)/d1(1,13)+sh(1,14)*d1(4,14)/d1(1,14)
     &          +sh(1,15)*d1(4,15)/d1(1,15)+sh(1,16)*d1(4,16)/d1(1,16)
     &          +sh(1,17)*d1(4,17)/d1(1,17)+sh(1,18)*d1(4,18)/d1(1,18)

	  wy_np1=sh(2,1)*d1(4,1)/d1(1,1)+sh(2,2)*d1(4,2)/d1(1,2)
     &          +sh(2,3)*d1(4,3)/d1(1,3)+sh(2,4)*d1(4,4)/d1(1,4)
     &          +sh(2,5)*d1(4,5)/d1(1,5)+sh(2,6)*d1(4,6)/d1(1,6)
     &          +sh(2,7)*d1(4,7)/d1(1,7)+sh(2,8)*d1(4,8)/d1(1,8)
     &          +sh(2,9)*d1(4,9)/d1(1,9)+sh(2,10)*d1(4,10)/d1(1,10)
     &          +sh(2,11)*d1(4,11)/d1(1,11)+sh(2,12)*d1(4,12)/d1(1,12)
     &          +sh(2,13)*d1(4,13)/d1(1,13)+sh(2,14)*d1(4,14)/d1(1,14)
     &          +sh(2,15)*d1(4,15)/d1(1,15)+sh(2,16)*d1(4,16)/d1(1,16)
     &          +sh(2,17)*d1(4,17)/d1(1,17)+sh(2,18)*d1(4,18)/d1(1,18)

	  uy_n  =sh(2,1)*do(2,1)/do(1,1)+sh(2,2)*do(2,2)/do(1,2)
     &          +sh(2,3)*do(2,3)/do(1,3)+sh(2,4)*do(2,4)/do(1,4)
     &          +sh(2,5)*do(2,5)/do(1,5)+sh(2,6)*do(2,6)/do(1,6)
     &          +sh(2,7)*do(2,7)/do(1,7)+sh(2,8)*do(2,8)/do(1,8)
     &          +sh(2,9)*do(2,9)/do(1,9)+sh(2,10)*do(2,10)/do(1,10)
     &          +sh(2,11)*do(2,11)/do(1,11)+sh(2,12)*do(2,12)/do(1,12)
     &          +sh(2,13)*do(2,13)/do(1,13)+sh(2,14)*do(2,14)/do(1,14)
     &          +sh(2,15)*do(2,15)/do(1,15)+sh(2,16)*do(2,16)/do(1,16)
     &          +sh(2,17)*do(2,17)/do(1,17)+sh(2,18)*do(2,18)/do(1,18)

	  uz_n  =sh(3,1)*do(2,1)/do(1,1)+sh(3,2)*do(2,2)/do(1,2)
     &          +sh(3,3)*do(2,3)/do(1,3)+sh(3,4)*do(2,4)/do(1,4)
     &          +sh(3,5)*do(2,5)/do(1,5)+sh(3,6)*do(2,6)/do(1,6)
     &          +sh(3,7)*do(2,7)/do(1,7)+sh(3,8)*do(2,8)/do(1,8)
     &          +sh(3,9)*do(2,9)/do(1,9)+sh(3,10)*do(2,10)/do(1,10)
     &          +sh(3,11)*do(2,11)/do(1,11)+sh(3,12)*do(2,12)/do(1,12)
     &          +sh(3,13)*do(2,13)/do(1,13)+sh(3,14)*do(2,14)/do(1,14)
     &          +sh(3,15)*do(2,15)/do(1,15)+sh(3,16)*do(2,16)/do(1,16)
     &          +sh(3,17)*do(2,17)/do(1,17)+sh(3,18)*do(2,18)/do(1,18)


	  vx_n  =sh(1,1)*do(3,1)/do(1,1)+sh(1,2)*do(3,2)/do(1,2)
     &          +sh(1,3)*do(3,3)/do(1,3)+sh(1,4)*do(3,4)/do(1,4)
     &          +sh(1,5)*do(3,5)/do(1,5)+sh(1,5)*do(3,6)/do(1,6)
     &          +sh(1,7)*do(3,7)/do(1,7)+sh(1,8)*do(3,8)/do(1,8)
     &          +sh(1,9)*do(3,9)/do(1,9)+sh(1,10)*do(3,10)/do(1,10)
     &          +sh(1,11)*do(3,11)/do(1,11)+sh(1,12)*do(3,12)/do(1,12)
     &          +sh(1,13)*do(3,13)/do(1,13)+sh(1,14)*do(3,14)/do(1,14)
     &          +sh(1,15)*do(3,15)/do(1,15)+sh(1,16)*do(3,16)/do(1,16)
     &          +sh(1,17)*do(3,17)/do(1,17)+sh(1,18)*do(3,18)/do(1,18)

	  vz_n  =sh(3,1)*do(3,1)/do(1,1)+sh(3,2)*do(3,2)/do(1,2)
     &          +sh(3,3)*do(3,3)/do(1,3)+sh(3,4)*do(3,4)/do(1,4)
     &          +sh(3,5)*do(3,5)/do(1,5)+sh(3,5)*do(3,6)/do(1,6)
     &          +sh(3,7)*do(3,7)/do(1,7)+sh(3,8)*do(3,8)/do(1,8)
     &          +sh(3,9)*do(3,9)/do(1,9)+sh(3,10)*do(3,10)/do(1,10)
     &          +sh(3,11)*do(3,11)/do(1,11)+sh(3,12)*do(3,12)/do(1,12)
     &          +sh(3,13)*do(3,13)/do(1,13)+sh(3,14)*do(3,14)/do(1,14)
     &          +sh(3,15)*do(3,15)/do(1,15)+sh(3,16)*do(3,16)/do(1,16)
     &          +sh(3,17)*do(3,17)/do(1,17)+sh(3,18)*do(3,18)/do(1,18)


	  wx_n  =sh(1,1)*do(4,1)/do(1,1)+sh(1,2)*do(4,2)/do(1,2)
     &          +sh(1,3)*do(4,3)/do(1,3)+sh(1,4)*do(4,4)/do(1,4)
     &          +sh(1,5)*do(4,5)/do(1,5)+sh(1,5)*do(4,6)/do(1,6)
     &          +sh(1,7)*do(4,7)/do(1,7)+sh(1,8)*do(4,8)/do(1,8)
     &          +sh(1,9)*do(4,9)/do(1,9)+sh(1,10)*do(4,10)/do(1,10)
     &          +sh(1,11)*do(4,11)/do(1,11)+sh(1,12)*do(4,12)/do(1,12)
     &          +sh(1,13)*do(4,13)/do(1,13)+sh(1,14)*do(4,14)/do(1,14)
     &          +sh(1,15)*do(4,15)/do(1,15)+sh(1,16)*do(4,16)/do(1,16)
     &          +sh(1,17)*do(4,17)/do(1,17)+sh(1,18)*do(4,18)/do(1,18)

	  wy_n  =sh(2,1)*do(4,1)/do(1,1)+sh(2,2)*do(4,2)/do(1,2)
     &          +sh(2,3)*do(4,3)/do(1,3)+sh(2,4)*do(4,4)/do(1,4)
     &          +sh(2,5)*do(4,5)/do(1,5)+sh(2,6)*do(4,6)/do(1,6)
     &          +sh(2,7)*do(4,7)/do(1,7)+sh(2,8)*do(4,8)/do(1,8)
     &          +sh(2,9)*do(4,9)/do(1,9)+sh(2,10)*do(4,10)/do(1,10)
     &          +sh(2,11)*do(4,11)/do(1,11)+sh(2,12)*do(4,12)/do(1,12)
     &          +sh(2,13)*do(4,13)/do(1,13)+sh(2,14)*do(4,14)/do(1,14)
     &          +sh(2,15)*do(4,15)/do(1,15)+sh(2,16)*do(4,16)/do(1,16)
     &          +sh(2,17)*do(4,17)/do(1,17)+sh(2,18)*do(4,18)/do(1,18)


	  uy = alpha * uy_np1 + oma * uy_n
	  uz = alpha * uz_np1 + oma * uz_n
	  vx = alpha * vx_np1 + oma * vx_n
	  vz = alpha * vz_np1 + oma * vz_n
	  wx = alpha * wx_np1 + oma * wx_n
	  wy = alpha * wy_np1 + oma * wy_n

	 temp1=sqrt(d1(1,1)*d1(1,1)+d1(2,1)*d1(2,1)+d1(3,1)*d1(3,1))
         temp2=sqrt(d1(1,2)*d1(1,2)+d1(2,2)*d1(2,2)+d1(3,2)*d1(3,2))
         temp3=sqrt(d1(1,3)*d1(1,3)+d1(2,3)*d1(2,3)+d1(3,3)*d1(3,3))
         temp4=sqrt(d1(1,4)*d1(1,4)+d1(2,4)*d1(2,4)+d1(3,4)*d1(3,4))
	 temp5=sqrt(d1(1,5)*d1(1,5)+d1(2,5)*d1(2,5)+d1(3,5)*d1(3,5))
         temp6=sqrt(d1(1,6)*d1(1,6)+d1(2,6)*d1(2,6)+d1(3,6)*d1(3,6))
	 temp7=sqrt(d1(1,7)*d1(1,7)+d1(2,7)*d1(2,7)+d1(3,7)*d1(3,7))
         temp8=sqrt(d1(1,8)*d1(1,8)+d1(2,8)*d1(2,8)+d1(3,8)*d1(3,8))
         temp9=sqrt(d1(1,9)*d1(1,9)+d1(2,9)*d1(2,9)+d1(3,9)*d1(3,9))
         temp10=sqrt(d1(1,10)*d1(1,10)+d1(2,10)*d1(2,10)+d1(3,10)*d1(3,10))
	 temp11=sqrt(d1(1,11)*d1(1,11)+d1(2,11)*d1(2,11)+d1(3,11)*d1(3,11))
         temp12=sqrt(d1(1,12)*d1(1,12)+d1(2,12)*d1(2,12)+d1(3,12)*d1(3,12))
	 temp13=sqrt(d1(1,13)*d1(1,13)+d1(2,13)*d1(2,13)+d1(3,13)*d1(3,13))
         temp14=sqrt(d1(1,14)*d1(1,14)+d1(2,14)*d1(2,14)+d1(3,14)*d1(3,14))
         temp15=sqrt(d1(1,15)*d1(1,15)+d1(2,15)*d1(2,15)+d1(3,15)*d1(3,15))
         temp16=sqrt(d1(1,16)*d1(1,16)+d1(2,16)*d1(2,16)+d1(3,16)*d1(3,16))
	 temp17=sqrt(d1(1,17)*d1(1,17)+d1(2,17)*d1(2,17)+d1(3,17)*d1(3,17))
         temp18=sqrt(d1(1,18)*d1(1,18)+d1(2,18)*d1(2,18)+d1(3,18)*d1(3,18))


	umag0=sh(0,1)*temp1+sh(0,2)*temp2+sh(0,3)*temp3+sh(0,4)*temp4
     &       +sh(0,5)*temp5+sh(0,6)*temp6+sh(0,7)*temp7+sh(0,8)*temp8
     &       +sh(0,9)*temp9+sh(0,10)*temp10+sh(0,11)*temp11+sh(0,12)*temp12
     &       +sh(0,13)*temp13+sh(0,14)*temp14+sh(0,15)*temp15+sh(0,16)*temp16
     &       +sh(0,17)*temp17+sh(0,18)*temp18
 
	umagx=sh(1,1)*temp1+sh(1,2)*temp2+sh(1,3)*temp3+sh(1,4)*temp4
     &       +sh(1,5)*temp5+sh(1,6)*temp6+sh(1,7)*temp7+sh(1,8)*temp8
     &       +sh(1,9)*temp9+sh(1,10)*temp10+sh(1,11)*temp11+sh(1,12)*temp12
     &       +sh(1,13)*temp13+sh(1,14)*temp14+sh(1,15)*temp15+sh(1,16)*temp16
     &       +sh(1,17)*temp17+sh(1,18)*temp18

	umagy=sh(2,1)*temp1+sh(2,2)*temp2+sh(2,3)*temp3+sh(2,4)*temp4
     &       +sh(2,5)*temp5+sh(2,6)*temp6+sh(2,7)*temp7+sh(2,8)*temp8
     &       +sh(2,9)*temp9+sh(2,10)*temp10+sh(2,11)*temp11+sh(2,12)*temp12
     &       +sh(2,13)*temp13+sh(2,14)*temp14+sh(2,15)*temp15+sh(2,16)*temp16
     &       +sh(2,17)*temp17+sh(2,18)*temp18

	umagz=sh(3,1)*temp1+sh(3,2)*temp2+sh(3,3)*temp3+sh(3,4)*temp4
     &       +sh(3,5)*temp5+sh(3,6)*temp6+sh(3,7)*temp7+sh(3,8)*temp8
     &       +sh(3,9)*temp9+sh(3,10)*temp10+sh(3,11)*temp11+sh(3,12)*temp12
     &       +sh(3,13)*temp13+sh(3,14)*temp14+sh(3,15)*temp15+sh(3,16)*temp16
     &       +sh(3,17)*temp17+sh(3,18)*temp18

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
	he_min=0.5*min(dis1,dis2,dis3)
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
             he_rgn = he_min
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
     &             + dabs(sx*sh(1,7)+sy*sh(2,7)+sz*sh(3,7))
     &             + dabs(sx*sh(1,8)+sy*sh(2,8)+sz*sh(3,8))
     &             + dabs(sx*sh(1,9)+sy*sh(2,9)+sz*sh(3,9))
     &             + dabs(sx*sh(1,10)+sy*sh(2,10)+sz*sh(3,10))
     &             + dabs(sx*sh(1,11)+sy*sh(2,11)+sz*sh(3,11))
     &             + dabs(sx*sh(1,12)+sy*sh(2,12)+sz*sh(3,12))
     &             + dabs(sx*sh(1,13)+sy*sh(2,13)+sz*sh(3,13))
     &             + dabs(sx*sh(1,14)+sy*sh(2,14)+sz*sh(3,14))
     &             + dabs(sx*sh(1,15)+sy*sh(2,15)+sz*sh(3,15))
     &             + dabs(sx*sh(1,16)+sy*sh(2,16)+sz*sh(3,16))
     &             + dabs(sx*sh(1,17)+sy*sh(2,17)+sz*sh(3,17))
     &             + dabs(sx*sh(1,18)+sy*sh(2,18)+sz*sh(3,18))
                                                                                   
                                                                                    
            he_ugn = 2.0/he_ugn
            else
            he_ugn = he_min
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
     &             + dabs(rx*sh(1,7)+ry*sh(2,7)+rz*sh(3,7))
     &             + dabs(rx*sh(1,8)+ry*sh(2,8)+rz*sh(3,8))
     &             + dabs(rx*sh(1,9)+ry*sh(2,9)+rz*sh(3,9))
     &             + dabs(rx*sh(1,10)+ry*sh(2,10)+rz*sh(3,10))
     &             + dabs(rx*sh(1,11)+ry*sh(2,11)+rz*sh(3,11))
     &             + dabs(rx*sh(1,12)+ry*sh(2,12)+rz*sh(3,12))
     &             + dabs(rx*sh(1,13)+ry*sh(2,13)+rz*sh(3,13))
     &             + dabs(rx*sh(1,14)+ry*sh(2,14)+rz*sh(3,14))
     &             + dabs(rx*sh(1,15)+ry*sh(2,15)+rz*sh(3,15))
     &             + dabs(rx*sh(1,16)+ry*sh(2,16)+rz*sh(3,16))
     &             + dabs(rx*sh(1,17)+ry*sh(2,17)+rz*sh(3,17))
     &             + dabs(rx*sh(1,18)+ry*sh(2,18)+rz*sh(3,18))
                                                                                    
            he_rgnt= 2.0/he_rgnt
                                                                                    
            else
            he_rgnt= he_min
       end if

       if (steady) then
		nutt = 0.0d0
		delnutt= 0.0d0
	else 
		nutt= (rho_np1*nut_np1 - rho_n*nut_n)/dt
		delnutt= rho_np1*delnut_np1/dt
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
	  gt   = 0.0d0
	  fac1 = 0.0d0
	  fac2 = 0.0d0
!	  gt   = dmin1(0.1d0,udif/vort_trip/he_trip)
!	  fac1 = distance*distance + gt*gt*dist_trip*dist_trip
!	  fac2 = ct2* vort_trip*vort_trip/udif/udif*fac1     

!	  t1 = dexp(-fac2)
!	  if (t1 .lt. 1.0e-20) then
!	    ft1 = 0.0000
!	  else
!	    ft1 = ct1*gt !*t1
!	  end if
!
!	  ft2 = ct3*dexp(-ct4*chi*chi)

	  ft1=0.0 ! for fully developed turbulence
	  ft2=0.0 ! for fully developed turbulence

	  r = nut/sbar/r_kappa/r_kappa/distance/distance

!	  if(r.ge.100.0d0) r=100.0d0
	  if(r.ge.1.0e+30)r=10.0
	  g = r + cw2*(r**6.0d0-r)
	  fw= g*((1.0d0+cw3**6.0d0)/
     $        (g**6.0d0+cw3**6.0d0))**(1.0d0/6.0d0)
	  p_function = rho*cb1*(1.0d0-ft2)*sgrad
	  r_function = rho*ft1*udif*udif      
	  q_function = rho*(-cw1*fw + (cb1/r_kappa/r_kappa*ft2)
     $       +(cb1/r_kappa/r_kappa*(1.0-ft2)*fv2))*(1.0/distance**2.0)    

!.......SUPG Stabilisation Term
	  tau_supg = 1.0/dsqrt((2.0d0*umagn/he_ugn)**2
     $             + (12.0d0*(dif_mol/rho)/he_rgnt/he_rgnt)**2)

!	fix it
!.......Galerkin stuff
         det  = abs(det)
         eff0 = wq3(ig)*det
	delnutub = delnutt
     $         + (rhou*delnutx+delnut*rhoux+rhov*delnuty+rhovy*delnut)
!	  eff0 = weis2(k)*det2
!	fix it 


	  do na =1,nen !na loop
	    b0 = sh(0,na)
	    b1 = sh(1,na)
	    b2 = sh(2,na)
	    b3 = sh(3,na)
	    if(stokes) then
	      ub =0.0d0
	    else
	      ub = u(1)*b1 + u(2)*b2 + u(3)*b3
	    end if

!.......Gravity terms
!	fix it 
	    p_t(1,na) = p_t(1,na)
     &            + (b0 + ub * tau_supg)
     &           *(delnutub 
     &       - p_function*delnut - 2.0*q_function*nut*delnut
     &       - 2.0*rho*cb2/sigma*
     &            (nutx*delnutx+nuty*delnuty+nutz*delnutz))*eff0
     &      + rho*(nut+dif_mol/rho)/sigma*
     $			(b1*delnutx+b2*delnuty+b3*delnutz)*eff0
     &      + rho*delnut/sigma*(b1*nutx+b2*nuty+b3*nutz)*eff0
!       added
     $      - rho*(nut+dif_mol/rho)/sigma*(ub*tau_supg)*
     $           (delnutxx + delnutyy+delnutzz)*eff0
!!     $      + rho*delnut/sigma*(nutxx +nutyy)*ub*tau_supg *eff0
!       added

c       write(80,180)p_t(1,na)
c       write(6,180)p_t(1,na)
c       stop
180             format(2e25.15)


	  enddo !na loop
	enddo ! quadrature loop

C...ASSEMBLE INTO GLOBAL VECTOR (STEP 1)
		do inl = 1,nen
			inode = ien(inl,ie)
			do idf = 1,ndf_t
			pgt(idf,inode) = pgt(idf,inode) + p_t(idf,inl)
			enddo
		enddo


        enddo ! element loop

        close(81)
        return
        end

