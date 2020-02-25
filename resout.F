c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine resout(d2)
	include "global.h"
	real* 8 d2(ndf,nnc)
        character*255 filerout
        integer i1,i2,i3,i4
	integer iolock,  iwho, itag, ilen, ipstat

        ibase = ichar('0')    !! integer value for char '0'

        filerout = "res.0000.00"

c
c dump number
c
        i4 = idisk/1000
        i3 = (idisk-i4*1000)/100
        i2 = (idisk-i4*1000-i3*100)/10
        i1 = (idisk-i4*1000-i3*100-i2*10)/1
c      
        i4 = i4 + ibase
        i3 = i3 + ibase
        i2 = i2 + ibase
        i1 = i1 + ibase
c
        filerout(5:5) = char(i4)
        filerout(6:6) = char(i3)
        filerout(7:7) = char(i2)
        filerout(8:8) = char(i1)
c
c iteration number
c
        i2 = iit/10
        i1 = (iit-i2*10)/1
c
        i2 = i2 + ibase
        i1 = i1 + ibase
c
        filerout(10:10) = char(i2)
        filerout(11:11) = char(i1)
c
        call fwdf(filerout,d2,ndf,nn)

	return
	end
