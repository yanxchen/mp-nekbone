Mixed-precision Nekbone
=======

> Nekbone solves a standard Poisson equation using a conjugate gradient iteration with a simple or spectral element multigrid preconditioner on a block or linear geometry. It exposes the principal computational kernel to reveal the essential elements of the algorithmic-architectural coupling that is pertinent to Nek5000.
> Original repository: https://github.com/Nek5000/Nekbone

It's an implementation of mixed-precision Nekbone.

## Usage
### Supported compilers
To compile mixed-precision Nekbone, both Fortran compiler and C compiler is required. Tested compilers include:
- Classic Flang: flang & clang (https://github.com/flang-compiler/flang)
- Intel compiler: ifort & icc
- Intel oneAPI compiler: ifx & icx
- Cray compiler: ftn & cc 
- GNU compiler: gfortran & gcc (old version like gfortran-5)

### Run test cases
- Test case w/o preconditioner: `test/nek`
- Test case w preconditioner: `test/mgrid`

Steps:

1. Specify compilers in `makenek` file: `F77` and `CC`, flang/clang is default.
2. (Optional) For parallel building: comment `IFMPI="false"` and specify MPI compiler.
3. (Optional) For changing optimization flags: modify `OPT_FLAGS_STD` and `OPT_FLAGS_MAG`.
4. Run the script: `./makenek`
5. Define the parameters for simulation in `data.rea`.
6. Run the program `./nekbone`.
