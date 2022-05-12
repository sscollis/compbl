#==============================================================================
#
#  Makefile for compbl (GCC)
#
#  Author:  Scott Collis
#
#  Revised: 2-12-2020 
#
#==============================================================================
NAME    = compbl
OPT     = -O2 
DEBUG   = -g
FFLAGS  = -cpp -fdefault-real-8 -fdefault-double-8 -ffixed-line-length-120 \
          -std=legacy -ffpe-trap=invalid,zero,overflow -ffpe-summary=none \
          -c $(OPT) $(DEBUG)
DEFINES = -DUSE_RUNGE
OFLAGS  = $(OPT) $(DEBUG) -o $(NAME)
LIB     =
FC      = gfortran 
F77     = gfortran 
#
# Default is currently to use Numerical-Recipes (if you have a valid license)
#
USE_NR  = 1
ifdef USE_NR
  ifeq ($(LIBNR_DIR),)
    LIBNR_DIR = $(HOME)/git/NR-utilities
  endif
  LIB += -L$(LIBNR_DIR) -lnr 
else
  $(error Currently must link with Numerical-Recipes)
endif

OBJECTS = compbl.o blkdta.o derivs.o input.o macprc.o \
getmat.o grid.o calc.o spline.o

$(NAME): $(OBJECTS)
	$(F77) $(OFLAGS) $(OBJECTS) $(LIB)

$(OBJECTS): common.h

.f.o:
	$(F77) $(FFLAGS) $(DEFINES) $*.f

clean:
	$(RM) *.o $(NAME)
