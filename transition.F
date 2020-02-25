        subroutine transition(xn,ien,dg1,dgo)

        implicit none
        include "global.h"

        integer i,j,idf,isd,inode,idof,inen,iet,ii
	integer isl,jj
        real*8  xn(nsd,nnl),xl1(nsd,nen)
        integer ibase,i1,i2
	integer ien(nen,nec)
	real*8 dg1(ndf,nnl),dgo(ndf,nnl)
	real*8 dlo(ndf,nen),dl1(ndf,nen)
	real*8 cv,R,gma,velc,uy,vx
	real*8 h1,h2,h3,h_min,h_max
!	real*8 umag_t(100),he_t(100),vort_t(100)
!	real*8 un,unp1,u_temp
	real*8 unp1(ndf),un(ndf),u(ndf)
!	real*8 tmp1,tmp2,tmp3

        cv = 1.0d0/((1.4d0)*(0.4d0))
	gma= 1.40d0
        R    = 1.0d0/gma

	isl = 0

        do i=1, ntelem
           if (my_pe .eq. itrne_pe(i)) then
!	   write(6,*) itrne_pe(i),ielem_t(i)
!	ODD NUMBER SLICE
	   do inen = 1,3
		inode = ien(inen,ielem_t(i))
		do isd= 1,nsd
	   		xl1(isd,inen) = xn(isd,inode)
!			write(6,*) xl1(isd,inen)
	   	end do
	   	do idof =1,ndf
	      	dl1(idof,inen) = dg1(idof,inode)
	      	dlo(idof,inen) = dgo(idof,inode)
!		write(6,*)dl1(idof,inen), dlo(idof,inen) 
	   	end do
	   end do 
	   isl = isl +1

	   umag_t(isl,i)=0.0d0

!.......Element Length
	  h1 = dsqrt((xl1(1,3)-xl1(1,1))*(xl1(1,3)-xl1(1,1))
     &             + (xl1(2,3)-xl1(2,1))*(xl1(2,3)-xl1(2,1)))
	  h2 = dsqrt((xl1(1,3)-xl1(1,2))*(xl1(1,3)-xl1(1,2))
     &             + (xl1(2,3)-xl1(2,2))*(xl1(2,3)-xl1(2,2)))
	  h3 = dsqrt((xl1(1,1)-xl1(1,2))*(xl1(1,1)-xl1(1,2))
     &             + (xl1(2,1)-xl1(2,2))*(xl1(2,1)-xl1(2,2)))

	  h_min = 0.5*dmin1(h1,h2,h3)
	  h_max = 0.5*dmax1(h1,h2,h3)
        
!.......he_trans()...Delta Xt = Grid Spacing Along the wall.
	  he_t(isl,i)=sqrt(h_max**2-h_min**2)
!	  write(6,*)he_t(isl,i)

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
            vort_t(isl,i)=dabs(vx-uy)
!	  write(6,*) vort_t(isl,i)
!	    write(6,*) umag_t(isl,i),he_t(isl,i),vort_t(isl,i)

!	EVEN NUMBER slice           
	   jj     = 0
	   do inen = 7,9
		inode = ien(inen,ielem_t(i))
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

	   isl = isl + 1 
	   umag_t(isl,i)=0.0d0

!.......Element Length
	  h1 = dsqrt((xl1(1,3)-xl1(1,1))*(xl1(1,3)-xl1(1,1))
     &             + (xl1(2,3)-xl1(2,1))*(xl1(2,3)-xl1(2,1)))
	  h2 = dsqrt((xl1(1,3)-xl1(1,2))*(xl1(1,3)-xl1(1,2))
     &             + (xl1(2,3)-xl1(2,2))*(xl1(2,3)-xl1(2,2)))
	  h3 = dsqrt((xl1(1,1)-xl1(1,2))*(xl1(1,1)-xl1(1,2))
     &             + (xl1(2,1)-xl1(2,2))*(xl1(2,1)-xl1(2,2)))

	  h_min = 0.5*dmin1(h1,h2,h3)
	  h_max = 0.5*dmax1(h1,h2,h3)
        
!.......he_trans()...Delta Xt = Grid Spacing Along the wall.
	  he_t(isl,i)=sqrt(h_max**2-h_min**2)
!	  write(6,*)he_t(isl,i)

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
            vort_t(isl,i)=dabs(vx-uy)
!	    write(6,*) vort_t(isl,i)
          
        end if
	end do

!	Last slice
        do i=ntelem, ntelem
           if (my_pe .eq. itrne_pe(i)) then
!	   write(6,*) itrne_pe(i),ielem_t(i)
	   jj = 0
	   do inen = 13,15
		inode = ien(inen,ielem_t(i))
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
	   isl = isl +1

	   umag_t(isl,i)=0.0d0

!.......Element Length
	  h1 = dsqrt((xl1(1,3)-xl1(1,1))*(xl1(1,3)-xl1(1,1))
     &             + (xl1(2,3)-xl1(2,1))*(xl1(2,3)-xl1(2,1)))
	  h2 = dsqrt((xl1(1,3)-xl1(1,2))*(xl1(1,3)-xl1(1,2))
     &             + (xl1(2,3)-xl1(2,2))*(xl1(2,3)-xl1(2,2)))
	  h3 = dsqrt((xl1(1,1)-xl1(1,2))*(xl1(1,1)-xl1(1,2))
     &             + (xl1(2,1)-xl1(2,2))*(xl1(2,1)-xl1(2,2)))

	  h_min = 0.5*dmin1(h1,h2,h3)
	  h_max = 0.5*dmax1(h1,h2,h3)
        
!.......he_trans()...Delta Xt = Grid Spacing Along the wall.
	  he_t(isl,i)=sqrt(h_max**2-h_min**2)
!	  write(6,*)he_t(isl,i)

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
            vort_t(isl,i)=dabs(vx-uy)
!	  write(6,*) vort_t(isl,i)
        end if
	end do

!	#if defined(UNI)
!#else
!            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!C...FIND TOTAL NUMBER OF EQUATIONS ON FIRST PROCESSOR AND BROADCAST TO REST
!	    call MPI_BCAST(vort_t,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!	    call MPI_BCAST(umag_t,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!	    call MPI_BCAST(he_t,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!#endif

        return
        end
