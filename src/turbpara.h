	real*8 cb1,cb2,ct1,ct2,ct3,ct4,cv1
	real*8 sigma,r_kappa,cw1,cw2,cw3
	real*8 cr1,cr2,cr3,crs,cs
	integer mtrans


        parameter (cb1=0.1355)
        parameter (cb2=0.622)
        parameter (ct1=1.0)
        parameter (ct2=2.0)
        parameter (ct3=1.2)
        parameter (ct4=0.5)
        parameter (cv1=7.1)
        parameter (sigma=0.66666666666667)
        parameter (r_kappa=0.41)
        parameter (cw1=cb1/r_kappa/r_kappa+(1.0+cb2)/sigma)
        parameter (cw2=0.3)
        parameter (cw3=2.0)
        parameter (cr1=1.00)
        parameter (cr2=12.00)
        parameter (cr3=1.00)
        parameter (cs=0.10d0)
	parameter (mtrans=30)
