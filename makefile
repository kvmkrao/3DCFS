SHELL=/bin/sh
#COM=obliquex
COM=AI_9bld5


OBJ =  main.o initialize.o parseinput.o echoinput.o diskout_sd.o \
       hypo2k.o gen.o norm.o \
       Ncomm.o archinfo.o gquad.o mkshlocal.o\
       blk1_iter_sd_mf.o \
       fgmres_smag.o tree.o \
       pn.o ckvol.o\
       baldwin.o f2snew.o \
       blas.o  error.o  filter.o  \
       updated.o  mylen.o \
       yon.o tstamp.o  mymalloc.o  csr.o keywr.o utility.o resout.o probe.o \
       blkrhs_turb_sd.o blkmvp_turb_sd.o  gen_turb.o fgmres_turb_sd.o  \
       test.o diskout_turb_sd.o rotation.o transform.o blk1_iter_sd.o \
       trip.o shape.o frifwd.o utility1.o

include ~/ARCH/makefile.i386
LDR=$(FC)
.SUFFIXES: .F
 
.F.o:
	$(FC) $(FFLAGS)  -c  $(?:.o=.F) 
.c.o:
	$(CC) $(CFLAGS)  -c $(?:.o=.c)

LFLAGS=-L ~/PARMEM -L ~/PARLIB -lfemC -lsvml

$(COM): $(OBJ)
	$(LDR) -o $(COM) $(OBJ)

install: $(COM)
	mv $(COM) $(HOME)/bin/$(ARCH)/$(COM)
uni:
	$(MAKE) COMMONLIB='$(COMMONLIBU)' FFLAGS='$(FFLAGSU)' CFLAGS='$(CFLAGSU)' CC='$(CCU)' FC='$(FCU)'

clean:
	rm -f $(OBJ) AI_*
