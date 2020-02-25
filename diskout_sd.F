c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine diskout_sd(d2)
	include "global.h"
	real* 8 d2(ndf,nnc)
        character*255 filexout, filedout
        integer i1,i2,i3,i4
	integer iolock,  iwho, itag, ilen, ipstat

        ibase = ichar('0')    !! integer value for char '0'

        filedout = "data.0000"
        filexout = "mxyz.0000"

        i4 = idisk/1000
        i3 = (idisk-i4*1000)/100
        i2 = (idisk-i4*1000-i3*100)/10
        i1 = (idisk-i4*1000-i3*100-i2*10)/1

        i4 = i4 + ibase
        i3 = i3 + ibase
        i2 = i2 + ibase
        i1 = i1 + ibase

        filedout(6:6) = char(i4)
        filedout(7:7) = char(i3)
        filedout(8:8) = char(i2)
        filedout(9:9) = char(i1)
        call fwdf(filedout,d2,ndf,nn)

	return
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine ndiskout_sd(d2)
        include "global.h"
        real* 8 d2(ndf,nnc)
        character*255 filexout, filedout
        integer i1,i2,i3,i4
        integer iolock,  iwho, itag, ilen, ipstat

        ibase = ichar('0')    !! integer value for char '0'

        filedout = "data.0000"
!!        filexout = "mxyz.0000"

        i4 = ndisk/1000
        i3 = (ndisk-i4*1000)/100
        i2 = (ndisk-i4*1000-i3*100)/10
        i1 = (ndisk-i4*1000-i3*100-i2*10)/1

        i4 = i4 + ibase
        i3 = i3 + ibase
        i2 = i2 + ibase
        i1 = i1 + ibase

	filedout(6:6) = char(i4)
        filedout(7:7) = char(i3)
        filedout(8:8) = char(i2)
        filedout(9:9) = char(i1)
        call fwdf(filedout,d2,ndf,nn)
c	write (6,*)'I M IN NDISKOUT',ndisk
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
