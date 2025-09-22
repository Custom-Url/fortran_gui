#!/bin/bash

# CAHN-HILLIARD PHASE-FIELD (J. Boisse + L.G.)
#=======================================================================

#--VARIABLES D'ENVIRONNEMENT POUR AMITEX
source ../env_amitex.sh

#-- OPTIMISATION openMPI - reseau OMNIPATH
#MPIRUNalias='mpirun -mca pml cm -mca mtl psm2'
#MPIRUNalias1='mpirun -np 1 -mca pml cm -mca mtl psm2'

MPIRUNalias='mpirun'
MPIRUNalias1='mpirun -np 1'

#DIR= TO_BE_FILLED
DIR=..

cd $DIR/cas_tests_PF

# version user_CahnHilliard3 : use of the STANDARD ALGORITHM with NLOC (the best way to do in that case)
#--------------------------------------------------------------------------------------------------------
AMITEX_PF=$DIR/libAmitex/src/user_CahnHilliard3/amitex_fftp

# 
# initial field : inclusion
#-----------------------------------------------------------------------
TEST="CahnHilliard"
TEST3="CahnHilliard_3"
ZVTK="microstructures/zone_CahnHilliard.vtk"
MATE="material/mat_CahnHilliard_3.xml"
LOAD="loading/load_CahnHilliard.xml"
ALGO="algorithm/algo_CahnHilliard_3.xml"

mkdir $DIR/resultats_PF/$TEST
$MPIRUNalias $AMITEX_PF  -nz $ZVTK -m $MATE -c $LOAD -a $ALGO -s ../resultats_PF/$TEST/$TEST3

#
# Decomposition spinodale (initial field : random around 0.5)
#-----------------------------------------------------------------------
TEST="CahnHilliard_spinodal"
TEST3="CahnHilliard_spinodal_3"
ZVTK="microstructures/zone_CahnHilliard_spinodal.vtk"
MATE="material/mat_CahnHilliard_spinodal_3.xml"
#LOAD same as before
#ALGO same as before

mkdir $DIR/resultats_PF/$TEST
$MPIRUNalias $AMITEX_PF  -nz $ZVTK -m $MATE -c $LOAD -a $ALGO -s ../resultats_PF/$TEST/$TEST3

# 
# initial field : inclusion   +   elastic coupling
#-----------------------------------------------------------------------
TEST="CahnHilliard_elast"
TEST3="CahnHilliard_elast_3"
ZVTK="microstructures/zone_CahnHilliard_elast.vtk"
MATE="material/mat_CahnHilliard_elast_3.xml"
LOAD="loading/load_CahnHilliard_elast.xml"
#ALGO same as before

mkdir $DIR/resultats_PF/$TEST
$MPIRUNalias $AMITEX_PF  -nz $ZVTK -m $MATE -c $LOAD -a $ALGO -s ../resultats_PF/$TEST/$TEST3

exit 0

