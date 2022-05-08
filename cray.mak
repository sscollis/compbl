#==============================================================================
#
#  Makefile for compbl (Cray)
#
#  Author:  Scott Collis
#
#  Revised: 10-13-94
#
#==============================================================================
NPROC  = 5
NAME   = compbl
DEBUG  =
FFLAGS = -c -N80 $(DEBUG)
OFLAGS = $(DEBUG) -o $(NAME)
LIB    =
COMP   = cf77

OBJECTS = compbl.o  blkdta.o  derivs.o  input.o  macprc.o  \
rk4.o  getmat.o  grid.o calc.o spline.o

$(NAME): $(OBJECTS)
	$(COMP) $(OFLAGS) $(OBJECTS) $(LIB)

$(OBJECTS): common.h

.f.o:
	$(COMP) $(FFLAGS) $*.f

clean:
	/bin/rm *.o
	
