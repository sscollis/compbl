#==============================================================================
#
#  Makefile for compbl (SGI)
#
#  Author:  Scott Collis
#
#  Revised: 10-13-94
#
#  -n32 compiles in the new 32 bit mode (faster than old ucode mode)
#  -d16 compiles in quad precision since code is in double precision
#
#==============================================================================
NPROC  = 5
NAME   = compbl
DEBUG  =
FFLAGS = -n32 -O2 -col120 -c $(DEBUG)
OFLAGS = -n32 -O2 $(DEBUG) -o $(NAME)
LIB    =
COMP   =  f77

OBJECTS = compbl.o  blkdta.o  derivs.o  input.o  macprc.o  \
rk4.o  getmat.o  grid.o calc.o spline.o

$(NAME): $(OBJECTS)
	$(COMP) $(OFLAGS) $(OBJECTS) $(LIB)

$(OBJECTS): common.h

.f.o:
	$(COMP) $(FFLAGS) $*.f

clean:
	/bin/rm *.o
