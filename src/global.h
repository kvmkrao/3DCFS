	include "mpf.h"

C...MACHINE SET UP
	integer request,ierr,my_pe,numprocs,numthreads,istat(MPI_STATUS_SIZE)
	common/cmpi/request,ierr,my_pe,numprocs,numthreads,istat
C...TIMER
	real*8 time(30)
	common/ctimer/time

	integer ncolor,my_color
	character*255 colorfile
	common/ccolor/ncolor,my_color,colorfile

C...CONSTANTS
      	real* 8    pi,epsilon,res_max, gamma
	integer mxeps
      	parameter (pi=3.141592654,epsilon=1.0e-8,res_max=1.0e+8)
	parameter(mxeps=9)
	parameter(gamma=1.4)

C...CONVERGENCE CRITERIA
	real*8 eps(mxeps)
	common/ceps/eps

C...PADDING SIZES
      	integer  nrngpad,ndfpad,ndf_tpad,nsdpad,nenpad,nquadpad,ndelpad
	parameter(nrngpad=24,ndfpad=5,nsdpad=4,nenpad=20,ndf_tpad=2)
	parameter(nquadpad=16,ndelpad=3)

C...QUADRATURE  DATA
	integer iquad,nquad3,nquad4
	real*8 sq3(0:nsdpad,nenpad,nquadpad),xq3(nsdpad,nquadpad),wq3(nquadpad)
	real*8 sq4(0:nsdpad,nenpad,nquadpad),xq4(nsdpad,nquadpad),wq4(nquadpad)
	common /cquad1/ iquad,nquad3,nquad4
	common/cquad2/ xq3,wq3,sq3,wq4,xq4,sq4

C...STABILITY PARAMETERS
	real* 8 delta(0:ndelpad) 
	common /cdelta/ delta

C...LOGICAL FLAGS
	logical debug, timing, restart, incdisp, ic_zero, steady, stokes, part
	logical dtintau, residual, freezstr, freezmesh, freezfluid,
     $	fsi_press_current
        logical sptime
	common /cexec/ debug, timing, restart, incdisp, ic_zero, steady, stokes, part, sptime
	common /cls/ dtintau,residual,freezstr,freezmesh,freezfluid,
     $	fsi_press_current
        real* 8 iunsteady, inostokes
	common /cflagslog/ iunsteady, inostokes

C...GMRES PARAMETERS
	integer inner,outer,iscaling,iprecond,isolver,numcl,ngmr,ikg
	common /cgmres1/ inner,outer,iscaling,iprecond,isolver,numcl,ngmr,ikg
	logical gmres
	common /lgmres/ gmres
	real*8 epsgm,epsrt
	common /reps/ epsgm,epsrt

C...MESH INFO
      	integer   nn,ne,nq,nq_t,nen,ndf,ndf_t,nsd,nc,nec,nnc,nqc,
     $	nqc_t,nrng,ncsr,nee,nnl,nnli,nnci,nn_s,ne_s,nen_s,nslices
      	common /cinfo/ nn,ne,nq,nq_t,nen,ndf,ndf_t,nsd,nc,nec,nnc,nqc,
     $	nqc_t,nrng,ncsr,nee,nnl,nnli,nnci,nn_s,ne_s,nen_s,nslices

C...PADDINGS
      	integer maxnec,maxnnc,maxnsd
      	common /cmaxi/ maxnec,maxnnc,maxnsd 


	integer force_surface
	common /force/  force_surface

C...TIME STEP AND ITERATION COUNTERS
	integer nts,nit,iit,its,nsrest,ntsbout_r,nsolution,idisk,ndisk,
     $	nitbref,ntss,itss,nsbdsk
	common /cproc/ its,nts,iit,nit,nsrest,ntsbout_r,nsolution, 
     $	idisk,ndisk,nitbref,ntss,itss,nsbdsk

	integer ntstart
	common ntstart

C...ELEMENT TYPES AND INFORMATION
	integer elem_tria_sd, elem_tria_st, elem_quad_sd, elem_quad_st
	integer elem_tetr_sd, elem_tetr_st, elem_hexa_sd, elem_hexa_st
	integer elem_wedg_sd,elem_wedg_quad_sd
	parameter (elem_tria_sd = 1, elem_tria_st = 5)
	parameter (elem_quad_sd = 2, elem_quad_st = 6)
	parameter (elem_tetr_sd = 3, elem_tetr_st = 7)
	parameter (elem_hexa_sd = 4, elem_hexa_st = 8)
	parameter (elem_wedg_sd = 9)
	parameter (elem_wedg_quad_sd = 10)
	integer           elem_type
	common /celement/ elem_type


C...REYNOLDS NUMBER
 	real* 8          reynolds, rho
	common /cdenvis/ reynolds, rho

C...BOUNDARY CONDITION AND INITIAL CONDITION INFO
	integer      bc(ndfpad,nrngpad),bc_t(ndf_tpad,nrngpad)
	real*8       ic(ndfpad,nrngpad),ic_t(ndf_tpad,nrngpad)

	common /c13/ bc,bc_t
	common /c14/ ic,ic_t
        integer  neum(nrngpad)
        common /c15/neum


C...GRAVITY VECTOR
	real* 8 gravity(3)
	common /cforprops/ gravity

C...FACEMAP DATA (LOC FACE NODE -> LOCAL ELEMENT NODE)
	integer maxneface, maxnnface, nelemtype
	parameter (maxneface = 6, maxnnface = 10, nelemtype = 15)
	integer facemap(maxneface,maxnnface,nelemtype), neface, nnface
	common /cfacemap/ facemap, neface, nnface

C...TIME INTEGRATION PARAMETER
	real* 8        alpha,oma
	common /ctdep/ alpha,oma

C...PRESSURE BOUNDARY CONDITION
	integer pnode
	common/cpresbc/pnode
C...PROBLEM DESCRIPTION
	character*255 title
	common /ctitle/ title

C...TIMESTEP AND ASSOCIATED PARAMETERS
	real* 8      tt, dt,dtinv,t_start
	common /ctt/ tt, dt,dtinv,t_start

C...CONVERGENCE VARIABLES
	real *8 del,residue
	common/cconv/del,residue

C...DIMENSIONAL ALIASES
	integer    xsd, ysd, zsd
	parameter (xsd = 1, ysd = 2, zsd = 3)
C...DOF ALIASES
	integer    udf, vdf, wdf, pdf
	parameter (udf = 1, vdf = 2, wdf = 3, pdf = 4)

C...TURBULENCE MODEL COEFFICIENT
	integer turb_model
	common /cturb/ turb_model
	real*8 he(1),turb(1),hwall(1),blyr
	pointer (heptr,he),(turbptr,turb),(hwallptr,hwall)
	common/cgeom/heptr,turbptr,hwallptr, blyr

C...VOLUME DATA
	real*8   volume(1),vmin,vmax
	pointer (volumeptr,volume)
	common/cvol/volumeptr,vmin,vmax

C...FILE VARIABLES
	character*255    mienfile, mxyzfile, mrngfile, dinfile1,dinfile2,
     $			 xinfile1,xinfile2,doutfile1,doutfile2,xoutfile1,
     $			 xoutfile2,dfile,xfile,sienfile,sxyzfile,doutfilevis,
     $			 minfile,moutfile,dinfile_t,doutfile_t,doutfilesum,doutfilepro,
     $			 mangfile,msangfile,meshbfile,mylfile,mwallnodefile,
     $                   mtransnfile,mtransefile,mtransdistfile
	common /c1_file/ mienfile, mxyzfile, mrngfile, dinfile1,dinfile2,
     $			 xinfile1,xinfile2,doutfile1,doutfile2,xoutfile1,
     $			 xoutfile2,dfile,xfile,sienfile,sxyzfile,doutfilevis,
     $			 minfile,moutfile,dinfile_t,doutfile_t,doutfilesum,doutfilepro,
     $			 mangfile,msangfile,meshbfile,mylfile,mwallnodefile,
     $                   mtransnfile,mtransefile,mtransdistfile

	integer  lookunit
	data     lookunit /7/

C...MESHMOVING STUFF
        logical meshmv, ffmeshmv, mesh_inc, nodalres
        integer mnq,mnqm,mouter,minner,mov_surface,sprng
        real*8 ym,pr,chi
        common/cmeshmv1/meshmv, ffmeshmv, mesh_inc, nodalres
	common/cmeshmv2/mnq,mnqm,mouter,minner
	common/cmeshmv3/ym,pr,chi,sprng

	integer      bcm(ndfpad,nrngpad)
	real*8       icm(ndfpad,nrngpad)
	common /cmeshmv4/ bcm
	common /cmeshmv5/ icm


C...INTERNAL PRESSURE
	real*8 inf_pres
	common/cinfp/inf_pres

C...Porosity Model
	logical porosity
	real*8 por_coeff
	common/por_model/ porosity,por_coeff

C....PROBES FOR VELOCITY AND PRESSURE AT CERTAIN NODES
	integer probepts,pmpts,ntranspts
	parameter(probepts=1000,pmpts=70000,ntranspts=51)
	integer nprobe, node_probe(probepts), probe_pe(probepts),
     $  node_pe(probepts),nely,nelz,nmeshb,cel(pmpts),
     $  cel_pe(pmpts),el_pe(pmpts)
        common/cprobe/nprobe, node_probe, probe_pe, node_pe,
     $  nely,nelz,nmeshb,cel_pe,el_pe
	integer ntelem,ielem_t(ntranspts),itrne_pe(ntranspts),
     $	ielem_pe(ntranspts)
	real*8 umag_t(51),he_t(51),vort_t(51)
	common/ctrip1/umag_t,he_t,vort_t
	common/ctelem/ntelem,ielem_t,itrne_pe,ielem_pe
	real*8 dy(pmpts),dz(pmpts),msh(18,pmpts),flowrate_pe(probepts),
     $	flowrate(probepts)
	common/cmesh/ dy,dz,msh,flowrate_pe,flowrate

C....ROTATION BOUNDARY CONDITIONS
	integer irotate
	real *8 vel,perc
	common irotate
	common/vel/perc
C....DOMAIN PARTITION VARIABLE
	integer ipartition
	common ipartition
	integer ipart(64)
	common ipart
	integer jpart(300)
	common jpart

C....THERMODYNAMIC CONDITIONS
	real*8 ShockCapFactor,WallTemp,BackPresInitial,BackPresFinal
	real*8 viscosity, rmach, density
	common /rthermo/ ShockCapFactor, WallTemp, BackPresInitial,
     $		BackPresFinal,viscosity,rmach
	logical shock_cap, Sutherland_Law, AdiabaticWall
	logical ApplyBackPressure, ApplyBleed
	common /lthermo/ shock_cap, Sutherland_Law, AdiabaticWall,
     $		ApplyBackPressure, ApplyBleed
	integer flag_he
	common /flaghe/ flag_he

C....MISC
	logical turbulent, formatted, explicit,trip
	common /misc/ turbulent, formatted, explicit,trip

C....DIAGONAL OR BLOCK DIAGONAL
	logical diag, block_diag
	common/block/ diag, block_diag


