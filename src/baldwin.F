	subroutine turbsetup(hnode,nnc,maxnnc,my_pe,numprocs)

	implicit none
	include "mpf.h"
	include "fmalloc.h"

        integer nnc,maxnnc,my_pe,numprocs
	integer hnode(nnc)
     

        integer nrcv,nrh,rproc,my_rh
        integer ria(1),rja(1),rcv(1),rh(1)
        pointer (riaptr,ria),(rjaptr,rja),(rcvptr,rcv),(rhptr,rh)

        integer nsnd,nsh,sproc,my_sh
        integer sia(1),sja(1),snd(1),sh(1)
        pointer (siaptr,sia),(sjaptr,sja),(sndptr,snd),(shptr,sh)

	integer perm(1)
	pointer (permptr,perm)

	common/ctcomm/nrcv,nrh,rproc,my_rh,riaptr,rjaptr,rcvptr,rhptr,
     $		      nsnd,nsh,sproc,my_sh,siaptr,sjaptr,sndptr,shptr,  
     $		      permptr

        integer SSTAT( MPI_STATUS_SIZE,512),RSTAT(MPI_STATUS_SIZE,512)
        integer IRQ(512),ISQ(512),ierr


        
        integer home(1),loc(1)
        pointer (homeptr,home),(locptr,loc)

        integer iloc,ihome,inode,jnode,iproc,jproc
	integer isrc,ides,i,j


	call error("TURBSETUP: STAGE 0",-999,.false.)

C...FIND HOMES OF NODES 
	homeptr = malloc(isize*nnc)
	locptr  = malloc(isize*nnc)
	do inode =1,nnc
		jnode = hnode(inode)
		ihome  = (jnode-1)/maxnnc
                iloc   = jnode - ihome*maxnnc
		home(inode) = ihome+1
		loc(inode)  = iloc
	enddo

C...FORM RCV 
	rcvptr = malloc(isize*numprocs)
	do iproc  = 1,numprocs
		rcv(iproc) =  0
	enddo
	do inode =1,nnc
		iproc = home(inode)
		rcv(iproc) = rcv(iproc) + 1
	enddo	

        nrcv = 0
        nrh   = 0
	do iproc  = 1,numprocs
		nrcv = nrcv + rcv(iproc)
		if(rcv(iproc).gt.0) nrh = nrh +1
	enddo

C...FORM RH
	if(nrh.gt.0) then

		rhptr = malloc(isize*nrh)
        	ihome   = 0
		do iproc  = 1,numprocs
			if(rcv(iproc).gt.0) then
				 ihome = ihome +1
				rh(ihome) = iproc
			endif
		enddo

C...FORM RIA
		riaptr = malloc(isize*(nrh+1))
        	ihome  = 1
		ria(ihome) = 1
		do iproc  = 1,numprocs
			if(rcv(iproc).ne.0) then
				ihome = ihome + 1
			        ria(ihome) = ria(ihome-1) + rcv(iproc)
			endif
		enddo

C...FORM RJA AND PERM
		rjaptr  = malloc(isize*nrcv)
		permptr = malloc(isize*nrcv)
		do inode =1,nnc
			iproc = home(inode)
			ihome = 1
			do while(rh(ihome).ne.iproc)
				ihome = ihome + 1
			enddo
			iloc  = ria(ihome)
			perm(inode) = iloc
			rja(iloc) = loc(inode)
			ria(ihome) = ria(ihome) + 1
		enddo
	

C...REFORM RIA
        	ihome  = 1
		ria(ihome) = 1
		do iproc  = 1,numprocs
			if(rcv(iproc).ne.0) then
				ihome = ihome + 1
			        ria(ihome) = ria(ihome-1) + rcv(iproc)
			endif
		enddo
	endif

	call free(homeptr)
	call free(locptr)
	
	
C...FIND OUT AMOUNTS OF DATA TO BE TRANSFERRED BETWEEN PROCESSORS
	call error("TURBSETUP: STAGE 1",-999,.false.)



C...FROM SND
	sndptr = malloc(isize*numprocs)

#if defined(UNI)
#else
	do iproc = 1,numprocs
		isrc = iproc-1
		call MPI_IRECV(snd(iproc),1,MPI_INTEGER,isrc,100,
     $                         MPI_COMM_WORLD,IRQ(iproc),ierr)
	enddo

	do iproc = 1,numprocs
		ides = iproc -1
		call MPI_ISEND(rcv(iproc),1,MPI_INTEGER,ides,100,
     $		               MPI_COMM_WORLD,ISQ(iproc),ierr)
	enddo

	call MPI_WAITALL(numprocs,ISQ,SSTAT,ierr)	
	call MPI_WAITALL(numprocs,IRQ,RSTAT,ierr)	
#endif

	call error("TURBSETUP: STAGE 2",-999,.false.)

        nsnd = 0
        nsh   = 0
	do iproc  = 1,numprocs
		nsnd = nsnd + snd(iproc)
		if(snd(iproc).gt.0) nsh = nsh +1
	enddo

C...FORM SH
	if(nsh.gt.0)  then
        	shptr = malloc(isize*nsh)
                ihome   = 0
		do iproc  = 1,numprocs
			if(snd(iproc).gt.0) then
		        	ihome = ihome +1
				sh(ihome) = iproc
			endif
		enddo

C...FORM SIA
		siaptr = malloc(isize*(nsh+1))
        	ihome  = 1
		sia(ihome) = 1
		do iproc  = 1,numprocs
			if(snd(iproc).gt.0) then
				ihome = ihome + 1
			        sia(ihome) = sia(ihome-1) + snd(iproc)
			endif
		enddo

C...RENUMBER SH WITH 0 BASED INDEX
		my_sh = -1
		do ihome = 1,nsh
			sh(ihome) = sh(ihome) -1
			if(sh(ihome).eq.my_pe)my_sh = ihome
		enddo
	endif

C...RENUMBER RH  WITH 0 BASED INDEX
	my_rh = -1
	do ihome = 1,nrh
		rh(ihome) = rh(ihome) -1
		if(rh(ihome).eq.my_pe)my_rh = ihome
	enddo

	if(my_rh*my_sh.lt.0) write(6,*) 'HOME ERROR ON ', my_pe

	call free(sndptr)
	call free(rcvptr)

C...REFORM SND AND RCV
	if(nsh.gt.0) then 
		sndptr = malloc(nsh*isize)
		do ihome = 1,nsh
			snd(ihome) = sia(ihome+1) - sia(ihome)
		enddo
	endif

	if(nrh.gt.0) then
		 rcvptr = malloc(nrh*isize)
		do ihome = 1,nrh
			rcv(ihome) = ria(ihome+1) - ria(ihome)
		enddo
	endif

	
C...RECV SJA
	if(nsnd.gt.0) sjaptr = malloc(isize*nsnd)
	rproc = 0
#if defined(UNI)
#else
	do iproc = 1,nsh
		iloc = sia(iproc)
		if(sh(iproc).ne.my_pe) then
		rproc = rproc + 1
		call MPI_IRECV(sja(iloc),snd(iproc),MPI_INTEGER,sh(iproc),100,
     $                         MPI_COMM_WORLD,IRQ(rproc),ierr)
		endif
	enddo
#endif

C...SEND RJA
	sproc = 0
#if defined(UNI)
#else
	do iproc = 1,nrh
		iloc = ria(iproc)
		if(rh(iproc).ne.my_pe) then
			sproc = sproc + 1
		        call MPI_ISEND(rja(iloc),rcv(iproc),MPI_INTEGER,
     $		         rh(iproc),100,MPI_COMM_WORLD,ISQ(sproc),ierr)
		else
	if(rcv(iproc).ne.snd(my_sh)) write(6,*) 'ERROR',my_pe,iproc,my_rh,my_sh
			j    = sia(my_sh)
			do i = 1,rcv(iproc)
				sja(j) = rja(iloc)
				j = j+1
				iloc = iloc + 1
			enddo
		endif
	enddo

	if(rproc.gt.0)call MPI_WAITALL(rproc,IRQ,RSTAT,ierr)	
	if(sproc.gt.0)call MPI_WAITALL(sproc,ISQ,SSTAT,ierr)	
#endif
	call error("TURBSETUP: STAGE 3",-999,.false.)

	return
	end

	subroutine turbsend(d,nx,nnc,my_pe,numprocs)

	implicit none
	include "mpf.h"
	include "fmalloc.h"

        integer nx,nnc,my_pe,numprocs
	real*8 d(nx,nnc)
     

        integer nrcv,nrh,rproc,my_rh
        integer ria(1),rja(1),rcv(1),rh(1)
        pointer (riaptr,ria),(rjaptr,rja),(rcvptr,rcv),(rhptr,rh)

        integer nsnd,nsh,sproc,my_sh
        integer sia(1),sja(1),snd(1),sh(1)
        pointer (siaptr,sia),(sjaptr,sja),(sndptr,snd),(shptr,sh)

	integer perm(1)
	pointer (permptr,perm)

	common/ctcomm/nrcv,nrh,rproc,my_rh,riaptr,rjaptr,rcvptr,rhptr,
     $		      nsnd,nsh,sproc,my_sh,siaptr,sjaptr,sndptr,shptr,
     $		      permptr  

        integer SSTAT( MPI_STATUS_SIZE,512),RSTAT(MPI_STATUS_SIZE,512)
        integer IRQ(512),ISQ(512),ierr

	real*8 sbuff(nx,1),rbuff(nx,1)
	pointer (sbuffptr,sbuff),(rbuffptr,rbuff)
	integer iloc,jloc,iproc,ist,ind,i,j

	if(nsnd.gt.0)sbuffptr  = malloc(fsize*nx*nsnd)
	rbuffptr  = malloc(fsize*nx*nrcv)

C...FILL IN  SEND BUFFER
	do iproc = 1,nsh
		ist = sia(iproc)
		ind = sia(iproc+1) -1
		do iloc = ist,ind
			jloc = sja(iloc)
			do i = 1,nx
				sbuff(i,iloc) = d(i,jloc)
			enddo
		enddo
	enddo
		
        
C...RECV IN RBUFF
#if defined(UNI)
#else
	rproc = 0
	do iproc = 1,nrh
		iloc = ria(iproc)
		if(rh(iproc).ne.my_pe) then
		rproc = rproc + 1
		call MPI_IRECV(rbuff(1,iloc),nx*rcv(iproc),
     $                         MPI_DOUBLE_PRECISION,rh(iproc),100,
     $			       MPI_COMM_WORLD,IRQ(rproc),ierr)
		endif
	enddo

C...SEND  SBUFF
	sproc = 0
	do iproc = 1,nsh
		iloc = sia(iproc)
		if(sh(iproc).ne.my_pe) then
			sproc = sproc + 1
		        call MPI_ISEND(sbuff(1,iloc),nx*snd(iproc),
     $			               MPI_DOUBLE_PRECISION,sh(iproc),100,
     $				       MPI_COMM_WORLD,ISQ(sproc),ierr)
		endif
	enddo


C...ON PN TRANSFER
	if(my_rh.ne.-1.and.my_sh.ne.-1) then
		iloc    = sia(my_sh)
		jloc    = ria(my_rh)
		do i = 1,rcv(my_rh)
			do j = 1,nx
				rbuff(j,jloc) = sbuff(j,iloc)
			enddo
			iloc = iloc + 1
			jloc = jloc + 1
		enddo
	endif

	if(rproc.gt.0)call MPI_WAITALL(rproc,IRQ,RSTAT,ierr)	
	if(sproc.gt.0)call MPI_WAITALL(sproc,ISQ,SSTAT,ierr)	

C...FILL OUT  RECV BUFFER
	do iloc = 1,nrcv
		jloc = perm(iloc)
		do i = 1,nx
			d(i,iloc) = rbuff(i,jloc)
		enddo
	enddo
#endif
		

	if(nsnd.gt.0)call free(sbuffptr)
	call free(rbuffptr)
	return
	end



	subroutine form_wall_func(ien,xl,dl,dlo,rng,xn,wfl,hwall,reynolds,
     $			  bdy,flag,nsd,ndf,nnl,nnc,nn,nen,nec,neface,my_pe,
     $			  numprocs,alpha,oma,maxnnc)
	implicit none
C	include "global.h"
	include "fmalloc.h"
	
	integer flag,nsd,ndf,nnl,nnc,nn,nen,nec,neface,my_pe,numprocs
	integer bdy,maxnnc
	real*8 reynolds,oma,alpha
	integer ien(nen,nec),rng(neface,nec)
	real*8  xl(nsd,nnl),xn(nsd,nnc),dl(ndf,nnl),dlo(ndf,nnl)
	real*8  wfl(nnl),hwall(nnc)


	integer hnode(1)
	pointer (hnodeptr,hnode)
	data    hnodeptr /0/
	save    hnodeptr


	real*8  sdiv(1)
	pointer (sdivptr,sdiv)
	data    sdivptr /0/
	save    sdivptr

	real*8 xnw(nsd,nnc)
	pointer (xnwptr,xnw)
	real*8 xx,yy,zz,ymn,ymx,rr

	real*8 wfn(1),vort(1)
	pointer (wfnptr,wfn),(vortptr,vort)

	character*255  turbfile
	real* 8 det,wt
        real* 8 sh(0:3,4),sql(0:3,4),xql(3),x(3,4),d(4,4),dq(3,3)
        real* 8 xr11,xr12,xr13,xr21,xr22,xr23,xr31,xr32,xr33
        real* 8 cf11,cf12,cf13,cf21,cf22,cf23,cf31,cf32,cf33
        real* 8 tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8
	real*8 xm1,ym1,zm1,xp1,yp1,zp1
	real*8 ux,uy,uz,vx,vy,vz,wx,wy,wz,vis,tmp
	real*8 txx,txy,txz,tyy,tyz,tzz,Aplus,yplus,uwall,turb_kappa

	integer i,iface,ie,je,inl,inode,isd,idf,iloc

C...SETUP COMMUNICATION, FORM HWALL AND BELEM
	if(flag.eq.1) then
C...FORM SDIV
                sdivptr = malloc(fsize*nnc)
		do i = 1,nnl
			wfl(i) = 1.0
		enddo
		call scatter(wfl,sdiv,1,1)
C...READ HNODE 
                hnodeptr = malloc(isize*nnc)
		turbfile = "mesh.hnode"
                call frif(turbfile,hnode,1,nn)
                call turbsetup(hnode,nnc,maxnnc,my_pe,numprocs)
	



	
C...FORM HWALL
	xnwptr = malloc(fsize*nsd*nnc)
	call equal(xn,xnw,nsd*nnc)
       	call turbsend(xnw,3,nnc,my_pe,numprocs)
	do i = 1,nnc
		xx =  xn(1,i)- xnw(1,i)
		yy =  xn(2,i)- xnw(2,i)
		zz =  xn(3,i)- xnw(3,i)
		hwall(i) = sqrt(xx*xx+yy*yy+zz*zz)
	enddo
	

	call free(xnwptr)
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

	vortptr = malloc(fsize*nnl)
	call fclear(vort,nnl)

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
		else
			include "sh3d4nac.h"
		endif

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
			ux =  ux + sh(1,inl)*d(1,inl)
			uy =  uy + sh(2,inl)*d(1,inl)
			uz =  uz + sh(3,inl)*d(1,inl)
			vx =  vx + sh(1,inl)*d(2,inl)
			vy =  vy + sh(2,inl)*d(2,inl)
			vz =  vz + sh(3,inl)*d(2,inl)
			wx =  wx + sh(1,inl)*d(3,inl)
			wy =  wy + sh(2,inl)*d(3,inl)
			wz =  wz + sh(3,inl)*d(3,inl)
		enddo


		xx  = vz -wy
		yy  = wx -uz
		zz  = uy -vx


		tmp = sqrt(xx*xx+yy*yy+zz*zz)

		do inl = 1,nen
			iloc = ien(inl,ie)
  	               if(tmp.gt.vort(iloc)) vort(iloc) =  tmp
                enddo
	enddo

C...FIND AVERAGE  VORTICITY AT THE NODE
	wfnptr = malloc(fsize*nnc)
	call scatter(vort,wfn,1,1)
	do i = 1,nnc
		wfn(i) = wfn(i)/sdiv(i)
	enddo
	call gather(vort,wfn,1)

	vis = 1.0/reynolds
	do i = 1,nnc
		wfn(i)   = sqrt(vis*wfn(i))
	enddo

 	call turbsend(wfn,1,nnc,my_pe,numprocs)

	turb_kappa = 0.41
	Aplus = 26.0
C	ymn = 10000.0
	do i = 1,nnc
		yplus = hwall(i)*wfn(i)/vis
C		if(yplus.lt.ymn) ymn = yplus
		wfn(i)  = 1.0 - exp(-yplus/Aplus) 
		wfn(i)  = wfn(i)*hwall(i)
	enddo
C	if(my_pe.eq.0) write(6,*) 'YPLUS', ymn
	call gather(wfl,wfn,1)
	call free(wfnptr)

	do i = 1,nnl
		tmp = turb_kappa*wfl(i)
		wfl(i) = tmp*tmp*vort(i)
	enddo

	call free(vortptr)
		
	return
	end
	




	
