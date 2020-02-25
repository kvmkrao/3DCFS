c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	malloc.h                                                             c
c	FORTRAN declarations for malloc/realloc/calloc/free wrappers         c
c	ISIZE - INTEGER size in bytes                                        c
c	FSIZE - REAL size in bytes                                           c
c	920601 - standardized                                                c
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer mymalloc, realloc, calloc
	integer null, isize, fsize
	parameter (null = 0, isize = 8, fsize = 8)
