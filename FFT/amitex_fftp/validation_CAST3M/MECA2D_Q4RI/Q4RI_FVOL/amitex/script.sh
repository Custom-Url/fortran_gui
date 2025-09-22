#!/bin/bash

#module load amitex_fftp/dev

AMITEX="/home/lgelebart/amitex_fftp/libAmitex/bin/amitex_fftp"

MATVTK="disque_1_R0.3_20.vtk"
MATVTK="disque_1_R0.3_129.vtk"
#MATVTK="fissure_CT_129.vtk"
MAT="mat_lin_pore.xml"
#MAT="mat_lin_homo.xml"
#MAT="mat_lin_fiss.xml"
#CHAR="char_def_imp.xml"
CHAR="char_sig_imp0.xml"
CHAR="char_sig_imp.xml"
ALGO="algo_default.xml"
OUT="r20"

rm tmp/*
mpirun -np 12 $AMITEX -nm $MATVTK -m $MAT -c $CHAR -a $ALGO -s tmp/$OUT
