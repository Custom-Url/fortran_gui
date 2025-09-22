#!/bin/bash

# Phase Field implemented by Y. Chen

#--VARIABLES D'ENVIRONNEMENT POUR AMITEX
source ../env_amitex.sh

MPIRUNalias='mpirun'
MPIRUNalias1='mpirun -np 1'

#DIR= TO_BE_FILLED
DIR=..

cd $DIR/cas_tests_PF

# Choice Miehe implementation context 
#-------------------------------------
#    step 1 (less integrated) = "USER ALGORITHM"
#    step 2 (medium integration) = "STANDARD ALGORITHM" with user-defined modeling (before_unpas and after_unpas)
#    step 3 (high integration) = "STANDARD ALGORITHM with 'NON-LOCAL' FRAMEWORK"
# here step 2
#====================================================================================

AMITEX_PF=$DIR/libAmitex/src/user_Miehe2/amitex_fftp
ALGO="algorithm/algo_Miehe2.xml"





# ====================================================================
# Cracked body under shear loads (xy=Mode II, xz/yz=ModeIII)
TEST="shearCrack2"
mkdir ../resultats_PF/$TEST

rm ../resultats_PF/$TEST/SC_xy*
$MPIRUNalias $AMITEX_PF -nm microstructures/shearCrack.vtk -m material/mat_shearCrack.xml -c loading/load_shearCrack_xy.xml -a $ALGO -s ../resultats_PF/$TEST/SC_xy

rm ../resultats_PF/$TEST/SC_xz*
$MPIRUNalias $AMITEX_PF -nm microstructures/shearCrack.vtk -m material/mat_shearCrack.xml -c loading/load_shearCrack_xz.xml -a $ALGO -s ../resultats_PF/$TEST/SC_xz

rm ../resultats_PF/$TEST/SC_yz*
$MPIRUNalias $AMITEX_PF -nm microstructures/shearCrack.vtk -m material/mat_shearCrack.xml -c loading/load_shearCrack_yz.xml -a $ALGO -s ../resultats_PF/$TEST/SC_yz


