#!/bin/bash
module load intel/parallel-studio-2020.4
#to compile fortran codes
# 1) use command ifort and flag -mkl
# 2) all extensions need to be .f90
# 3) example
#	ifort fund_const.f90 GaAs_parameters.f90 matrix_diag_MKL.f90 -mkl
for f in *.f; do
	ifort -c $f
done
ifort -mkl=parallel *.o -o d2mc #whatever name you want
