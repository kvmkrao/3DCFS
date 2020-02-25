	subroutine fgmres_turb_sd(bg_t,x1,d1,do,dg1,ien,id_t,rng,
     &                     d1_t,do_t,dg1_t,avloc_t,wg_t,ynoden,
     $			   umag_trans,he_trans,
     &                     vort_trans, dist_trans,
!     &                     dist_trans,
     $			   ntrans)
	implicit none
	include "global.h"
	include "fmalloc.h"

        real* 8 bg_t(1),x1(1),d1(1)
        real* 8 do(1),id_t(1),dg1(1)
        real*8  d1_t(1),do_t(1),dg1_t(1),avloc_t(ndf_t,nnl),ynoden(1)
        integer ien(1),rng(1)

        real*8 umag_trans(1),he_trans(1)
        real*8 vort_trans(1),dist_trans(1)
!        real*8 inode_trans(1)
	integer ntrans


	real*8 zg1(nqc_t)

	real* 8 vloc1(ndf_t,nnl), v1(nqc_t,inner+1), wg_t(nqc_t)
	real* 8 h(inner+1, inner)
	real* 8 y(inner+1)
	real* 8 cc(inner), ss(inner)

	pointer (vloc1ptr,vloc1),(v1ptr,v1)
	pointer (hptr,h),(yptr,y),(ccptr,cc),(ssptr,ss)
        pointer (zg1ptr,zg1)
	
	integer iqc,i,ii,i1,j,jj,j1,k,k1,l,igmres
	real* 8 rnorm,gam,tmpo,ysave,hsave,rnorm0,tmp1,tmp2



C...ALLOCATE MEMORY
	vloc1ptr  = malloc(fsize*ndf_t*nnl)
	v1ptr     = malloc(fsize*nqc_t*(inner+1))
	hptr     = malloc(fsize*inner*(inner+1))
	yptr     = malloc(fsize*(inner+1))
	ssptr    = malloc(fsize*inner)
	ccptr    = malloc(fsize*inner)
	zg1ptr    = malloc(fsize*nqc_t)

C...USING TWO SIDED SCALING
	do iqc = 1,nqc_t
!                wg_t(iqc) = id_t(iqc)/sqrt(abs(wg_t(iqc)))
C		wg1(iqc) = id(iqc)/wg1(iqc)
		wg_t(iqc) = id_t(iqc)
	enddo


C...COMPUTE RESIDUAL r = W**(-1/2) * (b - A * d)

	do iqc = 1,nqc_t
		bg_t(iqc) = bg_t(iqc)*wg_t(iqc)
	enddo
	call equal(bg_t,v1(1,1),nqc_t)
	call fclear(bg_t,nqc_t)

C...OUTER GMRES LOOP (igmres)
	igmres = 0 
 10	igmres = igmres + 1

C...CONVERGENCE CHECK
c		call tstamp(time,8)
		call getnorm(v1(1,1),v1(1,1),tmp1,nqc_t)
c		call tstamp(time,7)

		rnorm = sqrt(tmp1)


		if (igmres.eq.1) then
			rnorm0 = rnorm

		endif

		if (my_pe.eq.0) then
			write(6,1000) its, iit, igmres, rnorm/rnorm0 
			write(7,1000) its, iit, igmres, rnorm/rnorm0
		endif


		if (rnorm.le.eps(2).or.igmres.gt.outer) goto 700


C...FIRST KRYLOV VECTOR
		do iqc = 1,nqc_t
			v1(iqc,1) = v1(iqc,1)/rnorm
		enddo

C...ARNOLDI PROCESS
		do j=1,inner
			j1 = j + 1

C...RIGHT SCALING
			do iqc = 1,nqc_t
 				zg1(iqc) = v1(iqc,j)*wg_t(iqc)
C				zg1(iqc) = v1(iqc,j)
			enddo
C...MATVEC PRODUCT 
			if (debug) call error("fgmres:gather",-999,.false.)
c			call tstamp(time,8)
                	call gather (vloc1, zg1, ndf_t)
c			call tstamp(time,4)

			
			call blkmvp_turb_sd(ien,x1,d1,do,vloc1,avloc_t,
     &                     d1_t,do_t,ynoden,
!     &                     d1_t,do_t,ynoden,
     $                     umag_trans,he_trans,
     &                     vort_trans,dist_trans,
     $                     ntrans)

			if (debug) call error("fgmres:scatter",-999,.false.)
c			call tstamp(time,11)
			call scatter(avloc_t, v1(1,j1), ndf_t,1)
C			call tstamp(time,5)
C...LEFT SCALING
			do iqc = 1,nqc_t
 				v1(iqc,j1) =v1(iqc,j1)*wg_t(iqc)
			enddo



C...GRAM_SCHMIDT PROCESS
			do i=1,j
c				call tstamp(time,8)
				call getnorm(v1(1,j1),v1(1,i),tmp1,nqc_t)
				tmpo = tmp1 
c				call tstamp(time,7)
				h(i,j) = tmpo
				tmpo = -tmpo
				call saxpy(nqc_t, tmpo, v1(1,i), 1, v1(1,j1), 1)
			enddo

c			call tstamp(time,8)
			call getnorm(v1(1,j1),v1(1,j1),tmp1,nqc_t)
			tmpo = tmp1 
!			write(6,*) "fgmres_turb", tmpo,(1.0/tmpo)
c			call tstamp(time,7)
			tmpo = sqrt(tmpo)
			h(j1,j) = tmpo
!	murali
			tmpo = 1.0/tmpo
!	murali
			do iqc = 1,nqc_t
				v1(iqc,j1) = v1(iqc,j1) *tmpo
			enddo
		enddo			

C...SOLVE HESSENBERG SYSTEM ON ALL PROCESSORS 	H_m * y = beta * e_1

C...INITIALIZE REDUCED RESIDUAL	
		call fclear(y,inner+1)
		y(1) = rnorm

C...GIVENS ROTATIONS
		do j=1,inner
			j1 = j + 1
c			previously computed rotations on column j
			do i=2,j
				i1 = i - 1
				hsave = h(i1,j)
				h(i1,j) = + cc(i1) * hsave + ss(i1) * h(i,j)
				h(i ,j) = - ss(i1) * hsave + cc(i1) * h(i,j)
			enddo

c			new rotation on column j
			gam = sqrt(h(j,j)**2 + h(j1,j)**2)
			cc(j) = h(j,j) / gam
			ss(j) = h(j1,j) / gam
			h(j,j) = cc(j) * h(j,j) + ss(j) * h(j1,j)
c			note: under-diagonal term h(j+1,j) becomes 0

			y(j1) = - ss(j) * y(j)
			y(j ) = + cc(j) * y(j)
c			note: till now y(j+1) = 0
			rnorm = abs(y(j1))
		enddo


c		back substitution
		j = inner !should reach here straight from rotation loop
		y(j) = y(j)/h(j,j)
		do jj=2,j
			k = j - jj + 1
			k1 = k + 1
			ysave = y(k)
			do l=k1,j
				ysave = ysave - h(k,l) * y(l)
			enddo
			y(k) = ysave / h(k,k)
		enddo


C...COMPUTE GLOBAL ITERATE dg = dg + Z_m * y
		j = inner !(PVM only fix)
		do jj=1,j
			tmpo = y(jj)
			call saxpy(nqc_t, tmpo, v1(1,jj), 1, bg_t, 1)
		enddo


C...RECOVER RESIDUAL FOR RESTART
		if (igmres.le.outer) then
			do jj=1,j
				ii = j - jj + 2
				y(ii-1) = -ss(ii-1)*y(ii)
				y(ii) = cc(ii-1)*y(ii)
			end do

			do jj=1,j+1
				tmpo = y(jj)
				if (jj.eq.1) tmpo = tmpo - 1.0
				call saxpy(nqc_t, tmpo, v1(1,jj), 1, v1, 1)
			enddo
		endif

	goto 10

700   continue


C...GO BACK TO UNSCALED SYSTEM
	do iqc = 1,nqc_t
 		bg_t(iqc) =bg_t(iqc)*wg_t(iqc)
	enddo

 1000	format("fgmres:its ",i4,":iit ",i4,":iouter ",i4,": |r|/|ro| = ",e12.6)

	call free(vloc1ptr)
	call free(v1ptr)  
	call free(hptr)
	call free(yptr)
	call free(ssptr)
	call free(ccptr) 
	call free(zg1ptr)  


	return
	end
