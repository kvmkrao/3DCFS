C...MESHMOVING STUFF
        logical meshmv, ffmeshmv
        integer mnq,mnqm,mouter,minner,mov_surface
        real*8 ym,pr
        common/cmeshmv1/meshmv, ffmeshmv
	common/cmeshmv2/mnq,mnqm,mouter,minner
	common/cmeshmv3/ym,pr

	integer      bcm(ndfpad,nrngpad)
	real*8       icm(ndfpad,nrngpad)
	common /cmeshmv4/ bcm
	common /cmeshmv5/ icm
