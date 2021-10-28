	subroutine f2s_setup(maxnnc,numprocs,my_pe,nn_int)
	implicit none
	include "f2s.h"
	include "fmalloc.h"
	integer maxnnc,numprocs,my_pe,nn_int
	integer gnode,solid_node,fluid_home,fluid_node
	integer i,j,iloc,ihome,iproc
	integer nodeh(nn_f2s),nodel(nn_f2s),hsh(nn_f2s)
	!integer nodeh(1),nodel(1),hsh(1)
	pointer (nodehptr,nodeh),(nodelptr,nodel),(hshptr,hsh)

	character*255 file1,file2

	file1 = "surf.map"
	file2 = "surf.normal"

	nfs = 0
        nn_f2s = 2*nn_int

	if(my_pe.eq.0) then

C...READ IN NODEMAP AND NORMAL ARRAYS
		nodemapptr = malloc(nn_f2s*isize)
		normalptr  = malloc(nn_f2s*fsize)
!		call sfri(file1, nodemap,1,nn_f2s)
		call sfrd(file2,normal,1,nn_f2s)

C...PREPROCESS INFORMATION BETWEEN SURFACE AND VOLUME MESHES
        	ia_f2sptr = malloc((numprocs+1)*isize)
       	 	permptr   = malloc(nn_f2s*isize)
		nodehptr  = malloc(nn_f2s*isize)
		nodelptr  = malloc(nn_f2s*isize)
       	 	hshptr    = malloc(nn_f2s*isize)

C...FIND HOMES AND LOCAL NUMBERS (ON PN) FOR  FLUID NODES
		do solid_node = 1,nn_f2s
			gnode = nodemap(solid_node)
			fluid_home = (gnode -1)/maxnnc
			fluid_node = gnode - fluid_home*maxnnc
			nodeh(solid_node) = fluid_home
			nodel(solid_node) = fluid_node
		enddo


C...CONSTRUCT IA_F2S 
		do i = 1,numprocs + 1
			ia_f2s(i) = 0
		enddo

		do solid_node = 1,nn_f2s
			ihome = nodeh(solid_node) + 1
			ia_f2s(ihome) = ia_f2s(ihome) + 1
		enddo
		
		i = ia_f2s(1)
		ia_f2s(1) = 1

		do ihome = 2,numprocs + 1
			j = ia_f2s(ihome)
			ia_f2s(ihome) = ia_f2s(ihome-1) + i
			i = j
		enddo

C...CLEAR HSH, SORT NODES BY PROCESSORS AND CONSTRUCT PERMUTATION (PERM)
		do j = 1,nn_f2s
			hsh(j) = 0
		enddo	
		do solid_node = 1,nn_f2s
			ihome = nodeh(solid_node) + 1
			iloc  = ia_f2s(ihome)
			hsh(iloc) = nodel(solid_node)
			ia_f2s(ihome) = iloc + 1
			perm(iloc)    = solid_node
		enddo

C...RECONSTRUCT IA_F2S
		do ihome = 1,numprocs + 1
			ia_f2s(ihome) = 0
		enddo
		do solid_node = 1,nn_f2s
			ihome = nodeh(solid_node) + 1
			ia_f2s(ihome) = ia_f2s(ihome) + 1
		enddo
		
		i = ia_f2s(1)
		ia_f2s(1) = 1

		do ihome = 2,numprocs + 1
			j = ia_f2s(ihome)
			ia_f2s(ihome) = ia_f2s(ihome-1) + i
			i = j
		enddo

		call free(nodehptr)
		call free(nodelptr)

	endif

#if defined(UNI)
#else
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif


C...PN 0 COMMUNICATES NUMBER OF INTERFACE NODES LYING ON OTHER PNs
#if defined(UNI)
#else
	do iproc = 2,numprocs	
		if(my_pe.eq.0) then
			nfs = ia_f2s(iproc+1) - ia_f2s(iproc)
			call MPI_ISEND(nfs,1,MPI_INTEGER,iproc-1,100,
     $                          MPI_COMM_WORLD,request,ierr)
		endif
		if(my_pe.eq.iproc-1) then
                        call MPI_RECV(nfs,1,MPI_INTEGER,0,100,
     $                          MPI_COMM_WORLD,istat,ierr)
		endif
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	enddo
#endif

C...PN 0 COMMUNICATES NUMBER OF INTERFACE NODES TO ITSELF 
	if(my_pe.eq.0) nfs = ia_f2s(2) - ia_f2s(1)

C...ALLOCATE MEMORY FOR JA_F2S (ARRAY CONTAINING NODE NUMBERS ON INTERFACE)
	if(nfs.gt.0)  ja_f2sptr = malloc(nfs*isize)

C...FORM JA_F2S FOR OFF PN  (NODE NUMBERS ARE STORED IN HSH ARRAY)
#if defined(UNI)
#else
	do iproc = 2,numprocs	
		if(my_pe.eq.0) nfs = ia_f2s(iproc+1) - ia_f2s(iproc)
		if(my_pe.eq.0.and.nfs.gt.0) then
			iloc = ia_f2s(iproc)
			call MPI_ISEND(hsh(iloc),nfs,MPI_INTEGER,iproc-1
     $			,100,MPI_COMM_WORLD,request,ierr)
		endif
		if(my_pe.eq.iproc-1.and.nfs.gt.0) then
                        call MPI_RECV(ja_f2s,nfs,MPI_INTEGER,0,100,
     $                          MPI_COMM_WORLD,istat,ierr)
		endif
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	enddo
#endif


C...FORM JA_F2S FOR ON PN  (NODE NUMBERS ARE STORED IN HSH ARRAY)
	if(my_pe.eq.0)  then

		nfs = ia_f2s(2) - ia_f2s(1)

		do iloc = ia_f2s(1),ia_f2s(2) -1
			ja_f2s(iloc) = hsh(iloc)
		enddo

		call  free(hshptr)

	endif

	call error("F2S: SETUP COMM BETWEEN SURFACE AND VOLUME",-999,.false.)

	return
	end
		
	

C...SUBROUTINE TO EXCHANGE BOUNDARY CONDITIONS BETWEEN FLUID AND STRUCTURE
	subroutine f2s(flag,xf,df,dof,nsd,ndf,nn,my_pe,numprocs,
     $				alpha,dgc,dgcold,vgc,pg,tt,inf_pres)
	implicit none
	include "fmalloc.h"
	include "f2s.h"
	integer flag,ndf,nsd,nn,my_pe,numprocs
	real*8  alpha,beta,gamma,dtinv,opa,tt,inf_pres
	real*8 df(ndf,nn),xf(nsd,nn),dof(ndf,nn)
	real*8 dgc(nsd,1),vgc(nsd,1),pg(1)
	real*8 dgcold(nsd,1)
	integer nsd2f
	real*8 psl(1),pfl(1)
	real*8 vsl(2*nsd,1),vfl(2*nsd,1)
	!real*8 vsl(1,1),vfl(1,1)
	pointer(psptr,psl),(pfptr,pfl)
	pointer(vsptr,vsl),(vfptr,vfl)
        data psptr /0/, pfptr /0/, vsptr /0/, vfptr /0/
        save psptr,pfptr,vsptr,vfptr

	integer isd
	integer inode,jnode,knode,irq,isrc,ifs,iproc,iloc

C...HAVE TO SET THIS
	opa = 1.0 + alpha
	nsd2f = 2*nsd


C...TRANSFER PRESSURE  FROM VOLUME TO SURFACE 
	if (flag.eq.1) then

		if(my_pe.eq.0.and.psptr.eq.0) psptr = malloc(nn_f2s*fsize)

C...VOLUME -> VOLUME BUFFER (ON PN)
		if(nfs.gt.0) then
		       	if(pfptr.eq.0)  pfptr = malloc(nfs*fsize)
			do inode = 1,nfs
			jnode = ja_f2s(inode)
			pfl(inode) = opa*df(4,jnode) -alpha*dof(4,jnode)
			enddo

		endif

C...VOLUME BUFFER -> SURFACE BUFFER (ON PN)
		if(my_pe.eq.0.and.nfs.gt.0) then
			do inode = 1,nfs
				psl(inode) = pfl(inode) 
			enddo
		endif

C...VOLUME BUFFER -> SURFACE BUFFER (OFF PN)
#if defined(UNI)
#else
		do iproc = 2,numprocs
			if(my_pe.eq.iproc-1.and.nfs.gt.0) then
				call MPI_ISEND(pfl,nfs,MPI_DOUBLE_PRECISION,0,
     $					       100,MPI_COMM_WORLD,request,ierr)
			endif

			if(my_pe.eq.0) then
				ifs  = ia_f2s(iproc+1) - ia_f2s(iproc)
				if(ifs.gt.0) then
					iloc = ia_f2s(iproc)
					isrc = iproc -1
     				 	call MPI_RECV(psl(iloc),ifs,
     $					MPI_DOUBLE_PRECISION,
     $			               isrc,100,MPI_COMM_WORLD,istat,ierr)
				endif
			endif
		enddo

#endif

C...SURFACE BUFFER -> SURFACE (ON PN 0 USING PERM ARRAY)
C...ACCOUNT FOR DOUBLE VALUED NODES (USING NORMAL)
		if(my_pe.eq.0) then
			do inode = 1,nn_f2s/2
				pg(inode) = -inf_pres ! INTERNAL INFLATION
			enddo
			do inode = 1,nn_f2s
				jnode = (perm(inode)+1)/2
				pg(jnode) = pg(jnode) 
     $                                    + psl(inode)*normal(perm(inode))
			enddo
C			call free (psptr)
		endif

C		if(nfs.gt.0) call free(pfptr)

#if defined(UNI)
#else
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

	else

C...TRANSFER DISPLACEMENTS AND VELOCITY  FROM SURFACE TO VOLUME 
C...SURFACE BUFFER <- SURFACE  (ON PN 0 USING PERM)
		if(my_pe.eq.0) then
			if(vsptr.eq.0)vsptr = malloc(nsd2f*nn_f2s*fsize)
			do inode = 1,nn_f2s
				jnode = (perm(inode)+1)/2
				do isd = 1,nsd
					vsl(isd,inode) = vgc(isd,jnode) 
				enddo
				do isd =  1,nsd
				vsl(isd+ nsd,inode) = dgc(isd,jnode)-dgcold(isd,jnode)
				enddo
			enddo
		endif

		if(nfs.gt.0.and.vfptr.eq.0)vfptr = malloc(nsd2f*nfs*fsize)

C...VOLUME BUFFER <- SURFACE BUFFER (OFF PN) 
#if defined(UNI)
#else
		do iproc = 2,numprocs
		if(my_pe.eq.0) ifs=(ia_f2s(iproc+1)-ia_f2s(iproc))*nsd2f
			if(my_pe.eq.0.and.ifs.gt.0) then
			iloc = ia_f2s(iproc)
			isrc = iproc -1
     			call MPI_ISEND(vsl(1,iloc),ifs,MPI_DOUBLE_PRECISION,
     $				isrc,100,MPI_COMM_WORLD,request,ierr)
			endif
			if(my_pe.eq.iproc-1.and.nfs.gt.0) then
			call MPI_RECV(vfl,nfs*nsd2f,MPI_DOUBLE_PRECISION,0,
     $			       100,MPI_COMM_WORLD,istat,ierr)
			endif
		enddo
#endif


C...VOLUME BUFFER <- SURFACE BUFFER (ON PN) 
		if(my_pe.eq.0.and.nfs.gt.0) then
			do inode = 1,nfs
				do isd = 1,nsd2f
					vfl(isd,inode) = vsl(isd,inode) 
				enddo
			enddo
C			call free(vsptr)
		endif

C...VOLUME <- VOLUME BUFFER (ON PN)
		if(nfs.gt.0) then
			do  inode = 1,nfs
				jnode = ja_f2s(inode)
				do isd = 1,nsd
					df(isd,jnode) = vfl(isd,inode)
					!dof(isd,jnode) = vfl(isd,inode) ! added by vk.Apr09.02
				enddo
				do isd = 1,nsd
					xf(isd,jnode) = vfl(isd+ nsd,inode)
				enddo
					
			enddo
C			call free(vfptr)
		endif

#if defined(UNI)
#else
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
	endif
	

	return
	end
	


C
C THIS SUBROUTINE TRANSFERS SURFACE MESH INFORMATION TO THE VOLUME MESH
C     1. TO PREPARE THE CODE FOR CFD MESH MOVEMENT 
C           -TRANSFER SURFACE DISPLACEMENTS TO VOLUME
C           -TRANSFER THE AVERAGE SURFACE Z-DISPLACEMENT TO ALL 
C            OUTER BOUNDARY NODES
C     2. TRANSFERS SURFACE VELOCITIES TO VOLUME MESH

C...SUBROUTINE TO EXCHANGE BOUNDARY CONDITIONS BETWEEN FLUID AND STRUCTURE
	subroutine f2sFF(flag,xf,df,dof,nsd,ndf,nn,my_pe,numprocs,
     $				alpha,dgc,dgcold,vgc,pg,tt,inf_pres)
	implicit none
	include "fmalloc.h"
	include "f2s.h"
	integer flag,ndf,nsd,nn,my_pe,numprocs
	real*8  alpha,beta,gamma,dtinv,opa,tt,inf_pres
	real*8 df(ndf,nn),xf(nsd,nn),dof(ndf,nn)
	real*8 dgc(nsd,1),vgc(nsd,1),pg(1)
	real*8 dgcold(nsd,1)
	integer nsd2f
	real*8 psl(1),pfl(1)
	real*8 vsl(2*nsd,1),vfl(2*nsd,1)
	pointer(psptr,psl),(pfptr,pfl)
	pointer(vsptr,vsl),(vfptr,vfl)
        data psptr /0/, pfptr /0/, vsptr /0/, vfptr /0/
        save psptr,pfptr,vsptr,vfptr

	integer isd
	integer inode,jnode,knode,irq,isrc,ifs,iproc,iloc
	real*8 d_ave

C...HAVE TO SET THIS
	opa = 1.0 + alpha
	nsd2f = 2*nsd

	d_ave=0.0

C...TRANSFER DISPLACEMENTS AND VELOCITY  FROM SURFACE TO VOLUME 
C...SURFACE BUFFER <- SURFACE  (ON PN 0 USING PERM)
	if(my_pe.eq.0) then
		if(vsptr.eq.0)vsptr = malloc(nsd2f*nn_f2s*fsize)
		do inode = 1,nn_f2s
			jnode = (perm(inode)+1)/2
			do isd = 1,nsd
				vsl(isd,inode) = vgc(isd,jnode) 
			enddo
			do isd =  1,nsd
				vsl(isd+ nsd,inode) = dgc(isd,jnode)-dgcold(isd,jnode)
			enddo

			d_ave = d_ave + dgc(3,jnode)-dgcold(3,jnode)
		enddo

		d_ave = d_ave / nn_f2s
	endif

	if(nfs.gt.0.and.vfptr.eq.0)vfptr = malloc(nsd2f*nfs*fsize)

C...VOLUME BUFFER <- SURFACE BUFFER (OFF PN) 
#if defined(UNI)
#else
	do iproc = 2,numprocs
		if(my_pe.eq.0) ifs=(ia_f2s(iproc+1)-ia_f2s(iproc))*nsd2f
		if(my_pe.eq.0.and.ifs.gt.0) then
			iloc = ia_f2s(iproc)
			isrc = iproc -1
     			call MPI_ISEND(vsl(1,iloc),ifs,MPI_DOUBLE_PRECISION,
     $				isrc,100,MPI_COMM_WORLD,request,ierr)
C
	call MPI_ISEND(d_ave,1,MPI_DOUBLE_PRECISION,isrc,
     $		100,MPI_COMM_WORLD,request,ierr)
C
		endif
		if(my_pe.eq.iproc-1.and.nfs.gt.0) then
			call MPI_RECV(vfl,nfs*nsd2f,MPI_DOUBLE_PRECISION,0,
     $			       100,MPI_COMM_WORLD,istat,ierr)
C
	call MPI_RECV(d_ave,1,MPI_DOUBLE_PRECISION,0,
     $		100,MPI_COMM_WORLD,istat,ierr)
C
		endif
	enddo
#endif


C...VOLUME BUFFER <- SURFACE BUFFER (ON PN) 
C
	if(my_pe.eq.0.and.nfs.gt.0) then
		do inode = 1,nfs
			do isd = 1,nsd2f
				vfl(isd,inode) = vsl(isd,inode) 
			enddo
		enddo
	endif

C...VOLUME <- VOLUME BUFFER (ON PN)
C
C First, move all fluid nodes by d_ave
C
	do  inode = 1,nn
		xf(3,inode) = d_ave
	enddo
C
	if(nfs.gt.0) then
		do  inode = 1,nfs
			jnode = ja_f2s(inode)
			do isd = 1,nsd
				df(isd,jnode) = vfl(isd,inode)
				dof(isd,jnode) = vfl(isd,inode)  ! added by vk.Apr09.02
			enddo
C
C Second, move fluid surface nodes by sd displacement
C
			do isd = 1,nsd
				xf(isd,jnode) = vfl(isd+ nsd,inode)
			enddo

		enddo
C		call free(vfptr)
	endif

#if defined(UNI)
#else
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

	return
	end
	
