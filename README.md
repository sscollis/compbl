# Compressible Boundary Layer Solver

Solves the 2D compressible boundary-layer similarity equations.

Note:  Numerical Recipes `nr_rk4.f` is used which is commercially licensed 
and cannot be made public.  This can easily be replaced and should be.

To build use:

    ln -s gcc.mak Makefile
    make USE_NR=1 

Note that the build with fail if you don't supply `nr_rk4.f`.

S. Scott Collis\
Sun May  8 15:58:55 MDT 2022
