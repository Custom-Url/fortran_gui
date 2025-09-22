#!/bin/bash

# Phase Field implemented by Y. Chen

#--VARIABLES D'ENVIRONNEMENT POUR AMITEX
source ../env_amitex.sh

MPIRUNalias='mpirun -np 16'
MPIRUNalias1='mpirun -np 1'

#DIR= TO_BE_FILLED
DIR=..

cd $DIR/cas_tests_PF

# Choice Miehe implementation context 
#-------------------------------------
#    step 1 (less integrated) = "USER ALGORITHM"
#    step 2 (medium integration) = "STANDARD ALGORITHM" with user-defined modeling (before_unpas and after_unpas)
#    step 3 (high integration) = "STANDARD ALGORITHM with 'NON-LOCAL' FRAMEWORK"
# here step 3
#====================================================================================

AMITEX_PF=$DIR/libAmitex/src/user_Miehe3/amitex_fftp
ALGO="algorithm/algo_Miehe3.xml"


# ====================================
# Single Edge Notched Tensile specimen

TEST="SENT3"
mkdir ../resultats_PF/$TEST

rm ../resultats_PF/$TEST/SENT*
$MPIRUNalias1 $AMITEX_PF -nm microstructures/SENT.vtk -m material/mat_SENT3.xml -c loading/load_SENT.xml -a $ALGO -s ../resultats_PF/$TEST/SENT


