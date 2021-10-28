        subroutine fgmres_smag(bg1,x1,d1,avloc1,do,dg1,ien,id,rng,wg1,
     &                            d1_t,do_t,zlen,ang)
        implicit none
        include "global.h"
        include "fmalloc.h"
C............EXTERNAL VARIABLES

        real* 8 bg1(1),x1(1),d1(1),avloc1(ndf,nnl)
        real* 8 do(1),id(1),dg1(1),dg2(1)
        integer ien(1),rng(1)
        real*8  d1_t(1),do_t(1)

C............INTERNAL VARIABLES

        real*8 zg1(nqc)
        real*8 ang(nnc)

        real* 8 vloc1(ndf,nnl), v1(nqc,inner+1), wg1(nqc)
        real* 8 h(inner+1, inner)
        real* 8 y(inner+1)
        real* 8 cc(inner), ss(inner)

        pointer (vloc1ptr,vloc1),(v1ptr,v1)
        pointer (hptr,h),(yptr,y),(ccptr,cc),(ssptr,ss)
        pointer (zg1ptr,zg1)
        
        integer iqc,i,ii,i1,j,jj,j1,k,k1,l,igmres
        integer ip,iq
        real* 8 rnorm,gam,tmpo,ysave,hsave,rnorm0,tmp1,tmp2
        real *8 zlen


C...ALLOCATE MEMORY

        vloc1ptr  = malloc(fsize*ndf*nnl)
        v1ptr     = malloc(fsize*nqc*(inner+1))
        hptr     = malloc(fsize*inner*(inner+1))
        yptr     = malloc(fsize*(inner+1))
        ssptr    = malloc(fsize*inner)
        ccptr    = malloc(fsize*inner)
        zg1ptr    = malloc(fsize*nqc)

C...USING TWO SIDED SCALING
        do iqc = 1,nqc
!                wg1(iqc) = id(iqc)/sqrt(abs(wg1(iqc)))
!               wg1(iqc) = id(iqc)/wg1(iqc)
               wg1(iqc) = id(iqc)
        enddo


C...COMPUTE RESIDUAL r = W**(-1/2) * (b - A * d)
        do iqc = 1,nqc
                bg1(iqc) = bg1(iqc)*wg1(iqc)
        enddo
        call equal(bg1,v1(1,1),nqc)
        call fclear(bg1,nqc)

C...OUTER GMRES LOOP (igmres)
        igmres = 0 
 10     igmres = igmres + 1

C...CONVERGENCE CHECK
                call tstamp(time,8)
                call getnorm(v1(1,1),v1(1,1),tmp1,nqc)
                call tstamp(time,7)
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
                do iqc = 1,nqc
                        v1(iqc,1) = v1(iqc,1)/rnorm
                enddo

C...ARNOLDI PROCESS
                do j=1,inner
                        j1 = j + 1

C...RIGHT SCALING
                        do iqc = 1,nqc
                                zg1(iqc) = v1(iqc,j)*wg1(iqc)
C                               zg1(iqc) = v1(iqc,j)
                        enddo
C...MATVEC PRODUCT 
                        if (debug) call error("fgmres:gather",-999,.false.)
                        call tstamp(time,8)
                        call gather (vloc1, zg1, ndf)
                        call tstamp(time,4)
                        call blkmvp_iter_sd(ien,x1,d1,do,vloc1,
     $                      avloc1,ang,rng,d1_t,do_t)
!!                        call blkmvp_smag(ien,x1,d1,do,vloc1,avloc1,
!!     &                                d1_t,do_t,zlen)
                        if (debug) call error("fgmres:scatter",-999,.false.)
                        call tstamp(time,11)
                        call scatter(avloc1, v1(1,j1), ndf,1)
                        call tstamp(time,5)
C...LEFT SCALING
                        do iqc = 1,nqc
!!                                write(88,*) v1(iqc,j1),j1
!!                                write(90,*) wg1(iqc),j1
                                v1(iqc,j1) =v1(iqc,j1)*wg1(iqc)
!!                                write(89,*) v1(iqc,j1),j1
                        enddo



C...GRAM_SCHMIDT PROCESS
                        do i=1,j
                                call tstamp(time,8)
                                call getnorm(v1(1,j1),v1(1,i),tmp1,nqc)
                                tmpo = tmp1 
                                call tstamp(time,7)
                                h(i,j) = tmpo
                                tmpo = -tmpo
                                call saxpy(nqc, tmpo, v1(1,i), 1, v1(1,j1), 1)
                        enddo

                        call tstamp(time,8)
                        call getnorm(v1(1,j1),v1(1,j1),tmp1,nqc)
!!                        write(91,*) 'tmp1',tmp1,j1
!!                        write(91,*) (v1(ip,j1),ip=1,nqc)
                        tmpo = tmp1 
                        call tstamp(time,7)
                        tmpo = sqrt(tmpo)
                        h(j1,j) = tmpo
                        tmpo = 1.0/tmpo
                        do iqc = 1,nqc
                                v1(iqc,j1) = v1(iqc,j1) *tmpo
                        enddo
                enddo                   

C...SOLVE HESSENBERG SYSTEM ON ALL PROCESSORS   H_m * y = beta * e_1

C...INITIALIZE REDUCED RESIDUAL 
                call fclear(y,inner+1)
                y(1) = rnorm

C...GIVENS ROTATIONS
                do j=1,inner
                        j1 = j + 1
c                       previously computed rotations on column j
                        do i=2,j
                                i1 = i - 1
                                hsave = h(i1,j)
                                h(i1,j) = + cc(i1) * hsave + ss(i1) * h(i,j)
                                h(i ,j) = - ss(i1) * hsave + cc(i1) * h(i,j)
                        enddo

c                       new rotation on column j
                        gam = sqrt(h(j,j)**2 + h(j1,j)**2)
                        cc(j) = h(j,j) / gam
                        ss(j) = h(j1,j) / gam
                        h(j,j) = cc(j) * h(j,j) + ss(j) * h(j1,j)
c                       note: under-diagonal term h(j+1,j) becomes 0

                        y(j1) = - ss(j) * y(j)
                        y(j ) = + cc(j) * y(j)
c                       note: till now y(j+1) = 0
                        rnorm = abs(y(j1))
                enddo


c               back substitution
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
                        call saxpy(nqc, tmpo, v1(1,jj), 1, bg1, 1)
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
                                call saxpy(nqc, tmpo, v1(1,jj), 1, v1, 1)
                        enddo
                endif

        goto 10

700   continue

C...GO BACK TO UNSCALED SYSTEM
        do iqc = 1,nqc
                bg1(iqc) =bg1(iqc)*wg1(iqc)
        enddo

 1000   format("fgmres:its ",i4,":iit ",i4,":iouter ",i4,": |r|/|ro| = ",e12.6)

        call free(vloc1ptr)
        call free(v1ptr)  
        call free(hptr)
        call free(yptr)
        call free(ssptr)
        call free(ccptr) 
        call free(zg1ptr)  


        return
        end
