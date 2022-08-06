# Compressible Boundary Layer Similarity Solver

Solves the 2D compressible boundary-layer similarity equations.

Note:  Numerical Recipes `nr_rk4.f` is used which is commercially licensed 
and cannot be made public.  This can easily be replaced and should be.

To build use:
```bash
ln -s gcc.mak Makefile
make USE_NR=1 
```

Note that the build with fail if you don't supply `nr_rk4.f`.

# Running

A simple test case is run using
```bash
./compbl < compbl.inp
```

# Output files

File          |   Variables
--------------|--------------------------------
`output.dat`  | `sqrt(2)*eta, rho, u, v, 0, T`
`profile.dat` | `y, rho, u, v, 0, T`
`first.dat`   | `y, dudy, dvdy, 0, dTdy`
`second.dat`  | `y, d2udy2, d2vdy2, 0, d2Tdy2`

# Output on a grid

For this to work, you need to supply a `grid.dat` file (2d).

On output you will get:

File        |  Description
------------|------------------------------------------
`bl.dat`    | which is a Plot3d file
`mean.dat`  | which is an LNS mean flow file
`delta.dat` | displacement thickness as a function of `x`

Notes:
  1. I have noticed the use of the FORTRAN SNGL() function
     in writing the `bl.dat` and `mean.dat` files and this
     may need to be removed.
  2. I'm so confident that these SNGL() need to be removed
     I have done so but left the original lines as comments

S. Scott Collis\
flow.physics.simulation@gmail.com
