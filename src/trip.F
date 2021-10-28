        subroutine transition(xn,ien,dg1,dgo)

        implicit none
        include "global.h"

        integer i,j,idf,isd,inode,idof,inen,iet,ii
	integer isl,jj,ij,itemp,ia,ipr,iproc,pts,ntline
        real*8  xn(nsd,nnl),xl1(nsd,nen)
        integer ibase,i1,i2,nelemz,ne2d,iline,ipts
	integer ien(nen,nec)
	real*8 dg1(ndf,nnl),dgo(ndf,nnl)
	real*8 dlo(ndf,nen),dl1(ndf,nen)
	real*8 cv,R,gma,velc,uy,vx
	real*8 h1,h2,h3,h_min,h_max
!	real*8 umag_trip(100),he_trip(100),vort_trip(100)
	real*8 umag_trip,he_trip,vort_trip
!	real*8 un,unp1,u_temp
	real*8 unp1(ndf),un(ndf),u(ndf)
!	real*8 tmp1,tmp2,tmp3

        cv = 1.0d0/((1.4d0)*(0.4d0))
	gma= 1.40d0
        R    = 1.0d0/gma

!	write(6,*) "ntelem",ntelem
!	check

	nelemz = (nslices-1)
	ne2d   = ne/nelemz

	ntline = (ntelem/nelemz)
	ipts   = ntline*nslices


!	   do iproc =1,numprocs
!	   ipr = iproc-1

!	   ROOT PROCESSOR  STARTS HERE 
	   ipr = my_pe 

	   if(ipr.eq.0) then 
	   do i=1,ntelem
           if (ipr.eq.itrne_pe(i)) then
!	    write(6,*)ielem_t(i),itrne_pe(i),ielem_pe(i)
!	ODD NUMBER SLICE
	   do inen = 1,3
		inode = ien(inen,ielem_pe(i))
		do isd= 1,nsd
	   		xl1(isd,inen) = xn(isd,inode)
!			write(6,*) xl1(isd,inen)
	   	end do
!		write(6,*)xl1(1,inen),xl1(2,inen)
	   	do idof =1,ndf
	      	dl1(idof,inen) = dg1(idof,inode)
	      	dlo(idof,inen) = dgo(idof,inode)
!		write(6,*)dl1(idof,inen), dlo(idof,inen) 
	   	end do
	   end do 

	   iline = (i/nelemz) + 1
	   if(mod(i,nelemz).eq.0) iline = iline-1

	   itemp = (ielem_t(i)/ne2d)+ 1
	   if(mod(ielem_t(i),ne2d).eq.0) itemp = itemp - 1
	   isl = itemp
	   ia  = iline+isl-1

	   umag_t(ia)=0.0d0
!	   umag_trip = 0.0d0
!	   write(6,*) umag_t(iline,isl), iline,isl
!.......Element Length
	  h1 = dsqrt((xl1(1,3)-xl1(1,1))*(xl1(1,3)-xl1(1,1))
     &             + (xl1(2,3)-xl1(2,1))*(xl1(2,3)-xl1(2,1)))
	  h2 = dsqrt((xl1(1,3)-xl1(1,2))*(xl1(1,3)-xl1(1,2))
     &             + (xl1(2,3)-xl1(2,2))*(xl1(2,3)-xl1(2,2)))
	  h3 = dsqrt((xl1(1,1)-xl1(1,2))*(xl1(1,1)-xl1(1,2))
     &             + (xl1(2,1)-xl1(2,2))*(xl1(2,1)-xl1(2,2)))

	  h_min = dmin1(h1,h2,h3)
	  h_max = dmax1(h1,h2,h3)
        !.......he_trans()...Delta Xt = Grid Spacing Along the wall.
  	  he_t(ia)=he_t(ia) + sqrt(h_max**2-h_min**2)
!	  he_trip = sqrt(h_max**2-h_min**2)
!	  write(6,*)he_t(iline,isl),iline,isl
!.......magnitude of vorticity at Trip point
!.......Vorticity=(V,x-U,y).......V,x=0.0 and U,y.ne.0
           if(dl1(2,1)*dl1(2,2).eq.0.and.dl1(2,2)*dl1(2,3).eq.0
     $                             .and.dl1(2,3)*dl1(2,1).eq.0)then
              unp1(1)=(dl1(2,1)/dl1(1,1)+dl1(2,2)/dl1(1,2)
     $              +dl1(2,3)/dl1(1,3))
              un(1)  =(dlo(2,1)/dlo(1,1)+dlo(2,2)/dlo(1,2)
     $              +dlo(2,3)/dlo(1,3))
            else
              unp1(1)=(dl1(2,1)/dl1(1,1)+dl1(2,2)/dl1(1,2)
     $              +dl1(2,3)/dl1(1,3))/2.0
              un(1)  =(dlo(2,1)/dlo(1,1)+dlo(2,2)/dlo(1,2)
     $              +dlo(2,3)/dlo(1,3))/2.0
            endif

            u(1) = alpha*unp1(1)+oma*un(1)
            uy   = u(1)/h_min
            vx   = 0.0d0
            vort_t(ia)=vort_t(ia) + dabs(vx-uy)
!	    vort_trip = dabs(vx-uy)
!	    write(6,*) umag_t(ia), he_t(ia),vort_t(ia),ia

	end if 
	end do

!	RECEIVING DATA FROM WORKING PROCESSORS 

	   do ipr = 1,numprocs-1
	   do i=1,ntelem
           if (ipr.eq.itrne_pe(i)) then

	   iline = (i/nelemz) + 1
	   if(mod(i,nelemz).eq.0) iline = iline-1
	   itemp = (ielem_t(i)/ne2d)+ 1
	   if(mod(ielem_t(i),ne2d).eq.0) itemp = itemp - 1
	   isl = itemp
	   ia  = iline+isl-1

	call MPI_RECV(he_t(ia),1,MPI_DOUBLE_PRECISION,ipr,0,
     $	MPI_COMM_WORLD,istat,ierr)

!	call MPI_RECV(umag_t(ia),1,MPI_DOUBLE_PRECISION,ipr,0,
!     $	MPI_COMM_WORLD,istat,ierr)

	call MPI_RECV(vort_t(ia),1,MPI_DOUBLE_PRECISION,ipr,0,
     $	MPI_COMM_WORLD,istat,ierr)

	   if(i.eq.ntelem) then
	   iline = (i/nelemz) + 1
	   if(mod(i,nelemz).eq.0) iline = iline-1

	   itemp = (ielem_t(i)/ne2d)+ 1
	   isl = itemp
	   ia  = iline+isl-1

	call MPI_RECV(he_t(ia),1,MPI_DOUBLE_PRECISION,ipr,0,
     $	MPI_COMM_WORLD,istat,ierr)

!	call MPI_RECV(umag_t(ia),1,MPI_DOUBLE_PRECISION,ipr,0,
!     $	MPI_COMM_WORLD,istat,ierr)

	call MPI_RECV(vort_t(ia),1,MPI_DOUBLE_PRECISION,ipr,0,
     $	MPI_COMM_WORLD,istat,ierr)
	end if

	end if
	end do
	end do


	else 

!	   WORKING PROCESSORS LOOP STARTS HERE 
	   do i=1,ntelem
           if (ipr.eq.itrne_pe(i)) then
!	    write(6,*)ielem_t(i),itrne_pe(i),ielem_pe(i)
!	ODD NUMBER SLICE
	   do inen = 1,3
		inode = ien(inen,ielem_pe(i))
		do isd= 1,nsd
	   		xl1(isd,inen) = xn(isd,inode)
!			write(6,*) xl1(isd,inen)
	   	end do
!		write(6,*)xl1(1,inen),xl1(2,inen)
	   	do idof =1,ndf
	      	dl1(idof,inen) = dg1(idof,inode)
	      	dlo(idof,inen) = dgo(idof,inode)
!		write(6,*)dl1(idof,inen), dlo(idof,inen) 
	   	end do
	   end do 

	   iline = (i/nelemz) + 1
	   if(mod(i,nelemz).eq.0) iline = iline-1

	   itemp = (ielem_t(i)/ne2d)+ 1
	   if(mod(ielem_t(i),ne2d).eq.0) itemp = itemp - 1
	   isl = itemp
	   ia  = iline+isl-1

	   umag_t(ia)=0.0d0
!	   umag_trip = 0.0d0
!	   write(6,*) umag_t(iline,isl), iline,isl
!.......Element Length
	  h1 = dsqrt((xl1(1,3)-xl1(1,1))*(xl1(1,3)-xl1(1,1))
     &             + (xl1(2,3)-xl1(2,1))*(xl1(2,3)-xl1(2,1)))
	  h2 = dsqrt((xl1(1,3)-xl1(1,2))*(xl1(1,3)-xl1(1,2))
     &             + (xl1(2,3)-xl1(2,2))*(xl1(2,3)-xl1(2,2)))
	  h3 = dsqrt((xl1(1,1)-xl1(1,2))*(xl1(1,1)-xl1(1,2))
     &             + (xl1(2,1)-xl1(2,2))*(xl1(2,1)-xl1(2,2)))

	  h_min = dmin1(h1,h2,h3)
	  h_max = dmax1(h1,h2,h3)
        
!.......he_trans()...Delta Xt = Grid Spacing Along the wall.
  
	  he_t(ia)=he_t(ia) + sqrt(h_max**2-h_min**2)
!	  he_trip = sqrt(h_max**2-h_min**2)
!	  write(6,*)he_t(iline,isl),iline,isl
!.......magnitude of vorticity at Trip point
!.......Vorticity=(V,x-U,y).......V,x=0.0 and U,y.ne.0
           if(dl1(2,1)*dl1(2,2).eq.0.and.dl1(2,2)*dl1(2,3).eq.0
     $                             .and.dl1(2,3)*dl1(2,1).eq.0)then
              unp1(1)=(dl1(2,1)/dl1(1,1)+dl1(2,2)/dl1(1,2)
     $              +dl1(2,3)/dl1(1,3))
              un(1)  =(dlo(2,1)/dlo(1,1)+dlo(2,2)/dlo(1,2)
     $              +dlo(2,3)/dlo(1,3))
            else
              unp1(1)=(dl1(2,1)/dl1(1,1)+dl1(2,2)/dl1(1,2)
     $              +dl1(2,3)/dl1(1,3))/2.0
              un(1)  =(dlo(2,1)/dlo(1,1)+dlo(2,2)/dlo(1,2)
     $              +dlo(2,3)/dlo(1,3))/2.0
            endif

            u(1) = alpha*unp1(1)+oma*un(1)
            uy   = u(1)/h_min
            vx   = 0.0d0
            vort_t(ia)=vort_t(ia) + dabs(vx-uy)
!	    vort_trip = dabs(vx-uy)
!	    write(6,*) umag_t(ia), he_t(ia),vort_t(ia),ia

	call MPI_ISEND(he_t(ia),1,MPI_DOUBLE_PRECISION,0,0,
     $	MPI_COMM_WORLD,request,ierr)

!	call MPI_ISEND(umag_t(ia),1,MPI_DOUBLE_PRECISION,0,0,
!     $	MPI_COMM_WORLD,request,ierr)

	call MPI_ISEND(vort_t(ia),1,MPI_DOUBLE_PRECISION,0,0,
     $	MPI_COMM_WORLD,request,ierr)

!	EVEN NUMBER slice           
	   if(i.eq.ntelem) then
	   iline = (i/nelemz) + 1
	   if(mod(i,nelemz).eq.0) iline = iline-1
	   jj = 0
	   do inen = 4,6
		inode = ien(inen,ielem_pe(i))
		jj = jj+1
		do isd= 1,nsd
	   		xl1(isd,jj) = xn(isd,inode)
!			write(6,*) xl1(isd,inen)
	   	end do
	   	do idof =1,ndf
	      	dl1(idof,jj) = dg1(idof,inode)
	      	dlo(idof,jj) = dgo(idof,inode)
!		write(6,*)dl1(idof,jj), dlo(idof,jj) 
	   	end do
	   end do 

	   itemp = (ielem_t(i)/ne2d)+ 1
	   isl = itemp
	   ia  = iline+isl-1
	   umag_t(ia)=0.0d0
!.......Element Length
	  h1 = dsqrt((xl1(1,3)-xl1(1,1))*(xl1(1,3)-xl1(1,1))
     &             + (xl1(2,3)-xl1(2,1))*(xl1(2,3)-xl1(2,1)))
	  h2 = dsqrt((xl1(1,3)-xl1(1,2))*(xl1(1,3)-xl1(1,2))
     &             + (xl1(2,3)-xl1(2,2))*(xl1(2,3)-xl1(2,2)))
	  h3 = dsqrt((xl1(1,1)-xl1(1,2))*(xl1(1,1)-xl1(1,2))
     &             + (xl1(2,1)-xl1(2,2))*(xl1(2,1)-xl1(2,2)))

	  h_min = dmin1(h1,h2,h3)
	  h_max = dmax1(h1,h2,h3)
        
!.......he_trans()...Delta Xt = Grid Spacing Along the wall.
	  he_t(ia) =he_t(ia)+sqrt(h_max**2-h_min**2)

!.......magnitude of vorticity at Trip point
!.......Vorticity=(V,x-U,y).......V,x=0.0 and U,y.ne.0
           if(dl1(2,1)*dl1(2,2).eq.0.and.dl1(2,2)*dl1(2,3).eq.0
     $                             .and.dl1(2,3)*dl1(2,1).eq.0)then
              unp1(1)=(dl1(2,1)/dl1(1,1)+dl1(2,2)/dl1(1,2)
     $              +dl1(2,3)/dl1(1,3))
              un(1)  =(dlo(2,1)/dlo(1,1)+dlo(2,2)/dlo(1,2)
     $              +dlo(2,3)/dlo(1,3))
            else
              unp1(1)=(dl1(2,1)/dl1(1,1)+dl1(2,2)/dl1(1,2)
     $              +dl1(2,3)/dl1(1,3))/2.0
              un(1)  =(dlo(2,1)/dlo(1,1)+dlo(2,2)/dlo(1,2)
     $              +dlo(2,3)/dlo(1,3))/2.0
            endif

            u(1) = alpha*unp1(1)+oma*un(1)
            uy   = u(1)/h_min
            vx   = 0.0d0
            vort_t(ia)=vort_t(ia) + dabs(vx-uy)

	call MPI_ISEND(he_t(ia),1,MPI_DOUBLE_PRECISION,0,0,
     $	MPI_COMM_WORLD,request,ierr)

!	call MPI_ISEND(umag_t(ia),1,MPI_DOUBLE_PRECISION,0,0,
!     $	MPI_COMM_WORLD,request,ierr)

	call MPI_ISEND(vort_t(ia),1,MPI_DOUBLE_PRECISION,0,0,
     $	MPI_COMM_WORLD,request,ierr)
!	    write(6,*) umag_t(ia), he_t(ia),vort_t(ia),ia
	end if

	end if
	end do
!	WORKING PROCESSOR LOOP ENDS HERE 

	end if 
!	end do

#if defined(UNI)
#else
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

        return
        end

