# COMPBL:  Compressible Boundary Layer Similarity Solver

Solves the 2D compressible boundary-layer similarity equations.  Note that
this is a much simplier implementation than the Falkner-Skan-Cooke 
[FSC](https://github.com/sscollis/fsc) that is part of the 
[Flow Physics and Simulation (FPS)](https://github.com/flow-physics-simulation/flow-physics-simulation/edit/gh-pages/index.md) 
code suite.  As such, `compbl` is a good place to start exploring 
but it has limited capabilities in practise.

Note:  Numerical Recipes `nr_rk4.f` is used which is commercially licensed 
and cannot be made public.  This can easily be replaced and should be.

To build use:
```bash
ln -s gcc.mak Makefile
make USE_NR=1 
```
Note that the build with fail if you don't supply `nr_rk4.f` and it
would be great (and easy) to strip this dependency out in future 
versions.

## Running

A simple test case is run using
```bash
./compbl < compbl.inp
```

## Input File

The input file itself is a flat text file of the form:
```
0       ! 0 = constant viscosity, 1 = Sutherland
0.001   ! Viscosity, mu (makes delta=1)
0.3     ! Mach number
1.0     ! Prandtl number
0.0     ! Beta (Falkner-Skan parameter), 0=flat plate
1000    ! Re_delta (virtual origin)
n       ! write on a grid (y,n)
t       ! (t)emporal or (s)patial
```
Note that comments following the values, these comments
are ignored on input, and the input file itself it piped
to `compbl` as stdin as shown above.

## Output files

File          |   Variables
--------------|--------------------------------
`header.dat`  ! Basic info on the profiles
`output.dat`  | `sqrt(2)*eta, rho, u, v, 0, T`
`profile.dat` | `y, rho, u, v, 0, T`
`first.dat`   | `y, dudy, dvdy, 0, dTdy`
`second.dat`  | `y, d2udy2, d2vdy2, 0, d2Tdy2`

## Output on an `LNS3D` grid

For this to work, you need to supply a `grid.dat` file (2d)
that you use for LNS3D (which is actually in Plot3d format).
If that file is present, on output you will get:

File        |  Description
------------|------------------------------------------
`bl.dat`    | which is a Plot3d solution file
`mean.dat`  | which is an LNS3D mean flow file
`delta.dat` | displacement thickness as a function of `x`

With this capability, you can have a similarity solution 
BL flow as the base (mean) flow for subsequent `LNS3D` analysis.

Notes:
  1. I have noticed the use of the FORTRAN SNGL() function
     in writing the `bl.dat` and `mean.dat` files and this
     may need to be removed.
  2. I'm so confident that these SNGL() need to be removed
     I have done so but left the original lines as comments

S. Scott Collis\
flow.physics.simulation@gmail.com
