#==============================================================================
#
#  Makefile for compbl (GCC)
#
#  Author:  Scott Collis
#
#  Revised: 2-12-2020 
#
#==============================================================================
NAME   = compbl
OPT    = -O2 
DEBUG  = -g
FFLAGS = -cpp -fdefault-real-8 -fdefault-double-8 -ffixed-line-length-120 \
-std=legacy -c $(OPT) $(DEBUG)
OFLAGS = $(OPT) $(DEBUG) -o $(NAME)
LIB    =
FC     = gfortran 
F77    = gfortran 

OBJECTS = compbl.o blkdta.o derivs.o input.o macprc.o getmat.o grid.o \
					calc.o spline.o

ifdef USE_NR
	OBJECTS += nr_rk4.o
endif

$(NAME): $(OBJECTS)
	$(F77) $(OFLAGS) $(OBJECTS) $(LIB)

$(OBJECTS): common.h

.f.o:
	$(F77) $(FFLAGS) $*.f

clean:
	$(RM) *.o $(NAME)
