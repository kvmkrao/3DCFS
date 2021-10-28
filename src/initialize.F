!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!       INITIALIZE DATA
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        subroutine initialize
        implicit none

        include "global.h"

        integer i, j, idel
!.......WHAT?
!.......epsilons, but do individual indexes have special meaning?
        do i = 1,mxeps
            eps(i) = 1.0e-8
        enddo

!.......WHAT?
        blyr = 0.0

!FIXME  probably same as t
        tt = 0.0
        dt = 1.0
        t_start = 0.0

!.......WHAT?
        inf_pres = 0.0

!.......WHAT?
        minner   = 5
        mouter   = 1
        inner    = 5
        outer    = 1
        iprecond = 1
        isolver  = 1
        iscaling = 1
        
!.......WHAT?
        ncolor    = 1
        my_color  = 1
        colorfile = "mcolor"

        ndf    = 0
        nsd    = 0
        nn     = 0
        ne     = 0
        neface = 0      ! number of faces of an element (for wedge = 5)

!.......WHAT?
        nq    = 0       ! no of equations
        nn_s  = 0
        ne_s  = 0
        nen_s = 4
        nec   = 0       ! no of elements per (cluster)

        iquad      = 2
        ntss       = 1
        nts        = 1
        ntstart    = 1
	nslices    = 1
        nsbdsk     = 1
!        nt_start   = 1
        nit        = 1
        nsrest     = 1
        nsolution  = 1  
        nitbref    = 1
        idisk      = 0
        ndisk      = 0

        sprng = 0
        nrng  = 0

        force_surface = 0

        debug      = .false.
        meshmv     = .false.
        ffmeshmv   = .false.
        nodalres   = .false.
        mesh_inc   = .false.
        timing     = .false.
        restart    = .false.
        incdisp    = .false.
        freezstr   = .false.
        freezmesh  = .false.
        freezfluid = .false.
        ic_zero    = .false.
        steady     = .false.
        stokes     = .false.
        part       = .false.
        sptime     = .false.

        diag       = .false.
        block_diag = .false.

        gmres      = .false.
        epsgm      = 0.0
        epsrt      = 0.0
        numcl      = 1
        ngmr       = 0
        ikg        = 0

        WallTemp        = 1.50   ! Wall Temp, Implemented When AdiabaticWall=.false.
        ShockCapFactor  = 0.5
        BackPresInitial = 1.0
        BackPresFinal   = 1.0
        viscosity       = 0.0

        shock_cap         = .false.
        ApplyBackPressure = .false.
        AdiabaticWall     = .false.
        ApplyBleed        = .false.
        Sutherland_Law    = .false.

        rmach             = 0.0

        flag_he           = 1

        formatted         = .false.
        explicit          = .false.
        turbulent         = .false.
	trip              = .false.

        fsi_press_current= .true.

        pnode      = 0
        nprobe     = 0

        porosity=.false.
        por_coeff=0.0;

        reynolds = 1.0
        rho      = 1.0

        alpha = 0.5
        chi   = 1.0

        turb_model = 0

        mienfile    = "mien"
        mxyzfile    = "mxyz"
        mrngfile    = "mrng"
        dinfile1    = "data1.in"
        dinfile2    = "data2.in"
        dinfile_t   = "dataturb.in"
        xinfile1    = "mxyz1.in"
        xinfile2    = "mxyz2.in"
        doutfile1   = "data1.out"
        doutfile2   = "data2.out"
        doutfile_t  = "dataturb.out"
        doutfilesum = "datasum.out"     !***********
        doutfilepro = "datapro.out"     !***********
        doutfilevis = "eddyvis.out"     !***********
        xoutfile1   = "mxyz1.out"
        xoutfile2   = "mxyz2.out"
        dfile       = "mdata"
        xfile       = "mcoord"
        sienfile    = "surf.ien"
        sxyzfile    = "surf.xyz"
        minfile     = "disp.in"
        moutfile    = "disp.out"

        title = "Unknown problem"

        do j=1,nrngpad   
            do i=1,ndfpad   
                bc(i,j) = 0
                bcm(i,j) = 0
                ic(i,j) = 0.0
                icm(i,j) = 0.0
            enddo
            neum(j) = 0
        enddo
        do j=1,nrngpad   
            do i=1,ndf_tpad  
                bc_t(i,j) = 0
                ic_t(i,j) = 0.0
            enddo
            neum(j) = 0
        enddo

        

        gravity(1) = 0.0
        gravity(2) = 0.0
        gravity(3) = 0.0

        do idel=0,ndelpad
            delta(idel) = 1.0  
        end do

        return
        endsubroutine initialize

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SETUP FACEMAP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        block data


c       facemap maps local face node numbers to local element node numbers
        integer maxneface, maxnnface, nelemtype
        parameter (maxneface = 6, maxnnface = 10, nelemtype = 15)
        integer facemap(maxneface,maxnnface,nelemtype), neface, nnface
        data facemap /
c       triangle         |                 |                 |                 |
     & 2, 3, 1, 0, 0, 0, 3, 1, 2, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
c       quadrilateral    |                 |                 |                 |
     & 1, 2, 3, 4, 0, 0, 2, 3, 4, 1, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
c       tetrahedron      |                 |                 |                 |
     & 3, 1, 2, 3, 0, 0, 2, 2, 3, 1, 0, 0,
     & 1, 4, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
c       hexahedron       |                 |                 |                 |
     & 1, 1, 2, 3, 4, 5, 4, 2, 3, 4, 1, 6,
     & 3, 6, 7, 8, 5, 7, 2, 5, 6, 7, 8, 8,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
c       space-time triangle                |                 |                 |
     & 1, 2, 3, 0, 0, 0, 2, 3, 1, 0, 0, 0,
     & 5, 6, 4, 0, 0, 0, 4, 5, 6, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
c       space-time quadrilateral           |                 |                 |
     & 1, 2, 3, 4, 0, 0, 2, 3, 4, 1, 0, 0,
     & 6, 7, 8, 5, 0, 0, 5, 6, 7, 8, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
c       space-time tetrahedron             |                 |                 |
     & 3, 1, 2, 3, 0, 0, 2, 2, 3, 1, 0, 0,
     & 1, 4, 4, 4, 0, 0, 7, 5, 6, 7, 0, 0,
     & 6, 6, 7, 5, 0, 0, 5, 8, 8, 8, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
c       space-time hexahedron              |                 |                 |
     & 1, 1, 2, 3, 4, 5, 4, 2, 3, 4, 1, 6,
     & 3, 6, 7, 8, 5, 7, 2, 5, 6, 7, 8, 8,
     & 9, 9,10,11,12,13,12,10,11,12, 9,14,
     &11,14,15,16,13,15,10,13,14,15,16,16,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
c       wedge
     & 1, 4, 1, 1, 2, 0, 2, 6, 3, 4, 5, 0,
     & 3, 5, 6, 5, 6, 0, 0, 0, 4, 2, 3, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
c       wedge-18node                       |                 |                 |
     & 1,13, 1, 1, 2, 0, 4,18, 6, 7, 8, 0,
     & 2,15, 3,13,14, 0, 5,17, 9,16,17, 0,
     & 3,14,15,14,15, 0, 6,16,18, 8, 9, 0,
     & 0, 0,13, 2, 3, 0, 0, 0, 7, 4, 5, 0,
     & 0, 0,12,10,11, 0, 0, 0, 0, 0, 0, 0/

!     & 5, 8, 1, 1, 2,13, 3,11, 6,13,14,14,
!     & 6, 9, 3, 7, 8,15, 1,12,15,10,11, 0,
!     & 4, 7, 9, 8, 9, 0, 2,10,12,14,15, 0,
!     & 0, 0, 7, 2, 3, 0, 0, 0,13, 4, 5, 0/

       common /cfacemap/ facemap, neface, nnface

        end
