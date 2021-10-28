	subroutine trans_data(ien,xl,dl,dlo,rng,xn,reynolds,inode_trans,
     $         umagn_trans,vortn_trans,hen_trans,distn_trans,ntrans,
     $	      bdy,flag,nsd,ndf,nnl,nnc,nn,ne,nen,nec,neface,my_pe,
     $	        	  numprocs,alpha,oma,maxnnc,heigt,npart)
	implicit none
c	include "global.h"
	include "fmalloc.h"
	
	integer flag,nsd,ndf,nnl,nnc,nn,nen,nec,neface,my_pe,numprocs
	integer ne
	integer bdy,maxnnc
	real*8 reynolds,oma,alpha
	integer ien(nen,nec),rng(neface,nec),npart(20)
	real*8  xl(nsd,nnl),xn(nsd,nnc),dl(ndf,nnl),dlo(ndf,nnl)
	real*8  wfl(nnl),hwall(nnc)

	integer inode_trans(nnc),ntrans
	integer	ielem_trans(nnc)
	real*8  umagn_trans(nnc),hen_trans(nnc)
	real*8  vortn_trans(nnc),distn_trans(nnc)


	integer hnode(1)
	pointer (hnodeptr,hnode)
	data    hnodeptr /0/
	save    hnodeptr


	real*8  sdiv(1)
	pointer (sdivptr,sdiv)
	data    sdivptr /0/
	save    sdivptr

	real*8 xntrans(nsd,nnc),vtrans(1,nnc),utrans(1,nnc)
	real*8 htrans(1,nnc),heigt(nec)
	pointer (xntransptr,xntrans),(vtransptr,vtrans)
	pointer (utransptr,utrans),(htransptr,htrans)
	real*8 xx,yy,zz,ymn,ymx,rr

	real*8 wfn(1),vort(1),wfn_t(1,nnc)
	real*8 ufn(1),umagt(1),ufn_t(1,nnc)
	real*8 hfn(1),ht(1),hfn_t(1,nnc)
	pointer (wfnptr,wfn),(vortptr,vort),(wfn_tptr,wfn_t)
	pointer	(ufnptr,ufn),(umagtptr,umagt),(ufn_tptr,ufn_t)
	pointer (hfnptr,hfn),(htptr,ht),(hfn_tptr,hfn_t)

	character*255  turbfile,transfile
	character*255  elemfile
	real* 8 det,wt
        real* 8 sh(0:3,nen),sql(0:3,nen),xql(3),x(3,nen),d(4,nen),dq(3,3)
        real* 8 xr11,xr12,xr13,xr21,xr22,xr23,xr31,xr32,xr33
        real* 8 cf11,cf12,cf13,cf21,cf22,cf23,cf31,cf32,cf33
        real* 8 tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8
	real*8 xm1,ym1,zm1,xp1,yp1,zp1
	real*8 u0,v0,w0
	real*8 ux,uy,uz,vx,vy,vz,wx,wy,wz,vis,tmp
	real*8 txx,txy,txz,tyy,tyz,tzz,Aplus,yplus,uwall,turb_kappa
	real*8 hk,sx,sy,sz
	integer i,j,iface,ie,je,inl,inode,isd,idf,iloc

	integer ntrans1,knode,local_trans(ntrans)
	pointer (localtransptr,local_trans)

C...SETUP COMMUNICATION, FORM HWALL AND BELEM
	if(flag.eq.1) then
                transfile = "trans.hnode"
		elemfile = "trans.helem"
                call frif(transfile,inode_trans,1,nn,nec,npart)
		call frif(elemfile,ielem_trans,1,nn,nec,npart)


                call turbsetup(inode_trans,nnc,maxnnc,my_pe,numprocs)

	

C...FORM HWALL
	xntransptr = malloc(fsize*nsd*nnc)

	call equal(xn,xntrans,nsd*nnc)

       	call turbsend(xntrans,3,nnc,my_pe,numprocs)

	do i = 1,nnc
		xx =  xn(1,i)- xntrans(1,i)
		yy =  xn(2,i)- xntrans(2,i)
		zz =  xn(3,i)- xntrans(3,i)
		distn_trans(i) = sqrt(xx*xx+yy*yy+zz*zz)
	enddo

	call free(xntransptr)
	call error("FORMED HWALL",-999,.false.)
	return
	endif

	if(nen.eq.4) then
                xql(1)    = 0.25
                xql(2)    = 0.25
                xql(3)    = 0.25
                wt        = 0.16666667
                sql(0,1) = xql(1)
                sql(0,2) = xql(2)
                sql(0,3) = xql(3)
                sql(0,4) = 1 - xql(1) - xql(2) - xql(3)
c               local first derivatives
                        sql(1,1) = +1.0
                        sql(1,2) =  0.0
                        sql(1,3) =  0.0
                        sql(1,4) = -1.0
                                                                                                                   
                        sql(2,1) =  0.0
                        sql(2,2) = +1.0
                        sql(2,3) =  0.0
                        sql(2,4) = -1.0
                                                                                                                   
                        sql(3,1) =  0.0
                        sql(3,2) =  0.0
                        sql(3,3) = +1.0
                        sql(3,4) = -1.0

        elseif (nen.eq.6)then
			xql(1) = + 0.333333333333333
      			xql(2) = + 0.333333333333333
      			xql(3) = + 0.0
	          
                           wt  = 1.0e0

		sql(0,1) = (1.0 - xql(1) - xql(2))
     $                                  * (1.0 - xql(3)) / 2.0
                sql(0,2) = xql(1) * (1.0 - xql(3)) / 2.0
                sql(0,3) = xql(2) * (1.0 - xql(3)) / 2.0
                sql(0,4) = (1.0 - xql(1) - xql(2))
     $                                  * (1.0 + xql(3)) / 2.0
                sql(0,5) = xql(1) * (1.0 + xql(3)) / 2.0
                sql(0,6) = xql(2) * (1.0 + xql(3)) / 2.0

c...............Local shape function derivatives

		sql(1,1) = -1.0 * (1.0 - xql(3)) / 2.0
                sql(1,2) = +1.0 * (1.0 - xql(3)) / 2.0
                sql(1,3) = +0.0
                sql(1,4) = -1.0 * (1.0 + xql(3)) / 2.0
                sql(1,5) = +1.0 * (1.0 + xql(3)) / 2.0
                sql(1,6) = +0.0
                                                                                                                            
                                                                                                                            
                sql(2,1) = -1.0 * (1.0 - xql(3)) / 2.0
                sql(2,2) = +0.0
                sql(2,3) = +1.0 * (1.0 - xql(3)) / 2.0
                sql(2,4) = -1.0 * (1.0 + xql(3)) / 2.0
                sql(2,5) = +0.0
                sql(2,6) = +1.0 * (1.0 + xql(3)) / 2.0

	        sql(3,1) = -1.0 * (1.0 - xql(1)-xql(2))/2.0
                sql(3,2) = -1.0 * xql(1) / 2.0
                sql(3,3) = -1.0 * xql(2) / 2.0
                sql(3,4) = +1.0 * (1.0 - xql(1)-xql(2))/2.0
                sql(3,5) = +1.0 * xql(1) / 2.0
                sql(3,6) = +1.0 * xql(2) / 2.0	

	else

                xql(1)    = 0.0
                xql(2)    = 0.0
                xql(3)    = 0.0
                wt       = 8.0
		xm1      = 1-xql(1)
		xp1      = 1+xql(1)
		ym1      = 1-xql(2)
		yp1      = 1+xql(2)
		zm1      = 1-xql(3)
		zp1      = 1+xql(3)

                sql(0,1) = xm1*ym1*zm1/8 
                sql(0,2) = xp1*ym1*zm1/8 
                sql(0,3) = xp1*yp1*zm1/8 
                sql(0,4) = xm1*yp1*zm1/8 
                sql(0,5) = xm1*ym1*zp1/8 
                sql(0,6) = xp1*ym1*zp1/8 
                sql(0,7) = xp1*yp1*zp1/8 
                sql(0,8) = xm1*yp1*zp1/8 

                sql(1,1) = -ym1*zm1/8 
                sql(1,2) = +ym1*zm1/8 
                sql(1,3) = +yp1*zm1/8 
                sql(1,4) = -yp1*zm1/8 
                sql(1,5) = -ym1*zp1/8 
                sql(1,6) = +ym1*zp1/8 
                sql(1,7) = +yp1*zp1/8 
                sql(1,8) = -yp1*zp1/8 

                sql(2,1) = -xm1*zm1/8 
                sql(2,2) = -xp1*zm1/8 
                sql(2,3) = +xp1*zm1/8 
                sql(2,4) = +xm1*zm1/8 
                sql(2,5) = -xm1*zp1/8 
                sql(2,6) = -xp1*zp1/8 
                sql(2,7) = +xp1*zp1/8 
                sql(2,8) = +xm1*zp1/8 

                sql(3,1) = -xm1*ym1/8 
                sql(3,2) = -xp1*ym1/8 
                sql(3,3) = -xp1*yp1/8 
                sql(3,4) = -xm1*yp1/8 
                sql(3,5) = +xm1*ym1/8 
                sql(3,6) = +xp1*ym1/8 
                sql(3,7) = +xp1*yp1/8 
                sql(3,8) = +xm1*yp1/8 


        endif



	vortptr  = malloc(fsize*nnl)
	umagtptr = malloc(fsize*nnl)
	htptr    = malloc(fsize*nnl)
	call fclear(vort,nnl)
	call fclear(umagt,nnl)
	call fclear (ht,nnl)





C....FIND STRESS AT THE CENTROID OF ELEMENTS  ON THE SURFACE
	do ie = 1,nec
		do inl = 1,nen
			iloc = ien(inl,ie)
	        	do isd = 1,nsd
                                x(isd,inl)  = xl(isd,iloc)				
                        enddo

                        do idf = 1,ndf
                                d(idf,inl)  = alpha*dl(idf,iloc)
     $				            +  oma*dlo(idf,iloc)
                        enddo


		enddo
		if(nen.eq.8) then
			include "sh3d8nac.h"
		elseif(nen.eq.6)then
			include "sh3d6nac.h"
		else
			include "sh3d4nac.h"
		endif

		u0 = 0.0
		v0 = 0.0
		w0 = 0.0
		ux = 0.0
		uy = 0.0
		uz = 0.0
		vx = 0.0
		vy = 0.0
		vz = 0.0
		wx = 0.0
		wy = 0.0
		wz = 0.0

		do inl = 1,nen
			u0 = u0+sh(0,inl)*d(1,inl)
			v0 = v0+sh(0,inl)*d(2,inl)
			w0 = w0+sh(0,inl)*d(3,inl)
			ux =  ux + sh(1,inl)*d(1,inl)
			uy =  uy + sh(2,inl)*d(1,inl)
			uz =  uz + sh(3,inl)*d(1,inl)
			vx =  vx + sh(1,inl)*d(2,inl)
			vy =  vy + sh(2,inl)*d(2,inl)
			vz =  vz + sh(3,inl)*d(2,inl)
			wx =  wx + sh(1,inl)*d(3,inl)
			wy =  wy + sh(2,inl)*d(3,inl)
			wz =  wz + sh(3,inl)*d(3,inl)

		xx  = vz -wy
		yy  = wx -uz
		zz  = uy -vx
		tmp = sqrt(xx*xx+yy*yy+zz*zz)
		end do

		tmp1=sqrt(u0*u0+v0*v0+w0*w0)

		sx=u0/tmp1
		sy=v0/tmp1
		sz=w0/tmp1
		if(tmp1>1.0e-09)then
			tmp2=abs(sx*sh(1,1)+sy*sh(2,1)+sz*sh(3,1))
     &			   +abs(sx*sh(1,2)+sy*sh(2,2)+sz*sh(3,2))
     &			   +abs(sx*sh(1,3)+sy*sh(2,3)+sz*sh(3,3))
     &			   +abs(sx*sh(1,4)+sy*sh(2,4)+sz*sh(3,4))
     &                     +abs(sx*sh(1,5)+sy*sh(2,5)+sz*sh(3,5))
     &                     +abs(sx*sh(1,6)+sy*sh(2,6)+sz*sh(3,6))
			
			tmp2=2.0d0/tmp2
		else
			tmp2=heigt(ie)
		end if	




		do inl = 1,nen
			iloc = ien(inl,ie)
  	                if(tmp.gt.vort(iloc)) vort(iloc) =  tmp
c			if(tmp1.lt. umagt(iloc))umagt(iloc)=tmp1
			if(tmp2.gt.ht(iloc))ht(iloc)=tmp2
			u0=d(1,inl)
			v0=d(2,inl)
			w0=d(3,inl)
			tmp1=dsqrt(u0*u0+v0*v0+w0*w0)
			umagt(iloc)=tmp1
                enddo
	enddo

	
C...FIND  VORTICITY AT THE  CLOSEST TRANSITION NODE TO A GIVEN NODE
	
	wfnptr = malloc(fsize*nnc)
	wfn_tptr=malloc(fsize*nnc)
	
	
	ufnptr = malloc(fsize*nnc)
	ufn_tptr=malloc(fsize*nnc)

	hfnptr=malloc(fsize*nnc)
	hfn_tptr=malloc(fsize*nnc)

	call scatter(vort,wfn,1,1)
	call scatter(umagt,ufn,1,1)
	

	call scatter(ht,hfn,1,1)
	
	do i=1,nnc
		wfn_t(1,i)=wfn(i)
		ufn_t(1,i)=ufn(i)
		hfn_t(1,i)=hfn(i)
	end do


	vtransptr = malloc(fsize*nnc)
	utransptr =malloc(fsize*nnc)
	htransptr =malloc(fsize*nnc)
	call equal(wfn_t,vtrans,1*nnc)
	call equal(ufn_t,utrans,1*nnc)
	call equal(hfn_t,htrans,1*nnc)

       	call turbsend(vtrans,1,nnc,my_pe,numprocs)
	call turbsend(utrans,1,nnc,my_pe,numprocs)
	call turbsend(htrans,1,nnc,my_pe,numprocs)
	

	do i = 1,nnc
		vortn_trans(i)=vtrans(1,i)
		umagn_trans(i)=utrans(1,i)
		hen_trans(i)=htrans(1,i)
	enddo




	call free(vtransptr)
	call free(utransptr)
	call free(vortptr)
	call free(umagtptr)
	call free(wfn_tptr)
	call free(ufn_tptr)
	call free(htransptr)
	call free(htptr)
	call free(hfn_tptr)






	return
	end
	




	
