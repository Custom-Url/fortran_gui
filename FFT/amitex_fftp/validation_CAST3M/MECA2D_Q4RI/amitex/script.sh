#!/bin/bash

module load amitex_fftp/6.2.0-gcc482 

MATVTK="disque_1_R0.3_20.vtk"
MAT="mat_lin.xml"
CHAR="char_def_imp.xml"
ALGO="algo_default.xml"
OUT="r20"

mpirun amitex_fftp -nm $MATVTK -m $MAT -c $CHAR -a $ALGO -s tmp/$OUT
