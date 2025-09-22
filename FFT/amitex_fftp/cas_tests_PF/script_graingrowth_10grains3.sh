#!/bin/bash

# Phase Field J. Boisse

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

# version user_graingrowth3 : use of the STANDARD ALGORITHM with NLOC (the best way to do in that case)
#--------------------------------------------------------------------------------------------------------
AMITEX_PF=$DIR/libAmitex/src/user_graingrowth3/amitex_fftp

# 
# initial fields : 10 inclusions
#-----------------------------------------------------------------------
TEST="graingrowth_10grains"
TEST3="graingrowth_10grains_3"
ZVTK="microstructures/zone_graingrowth_10grains.vtk"
MATE="material/mat_graingrowth_10grains_3.xml"
LOAD="loading/load_graingrowth_10grains.xml"
ALGO="algorithm/algo_graingrowth_3.xml"

mkdir $DIR/resultats_PF/$TEST
$MPIRUNalias $AMITEX_PF  -nz $ZVTK -m $MATE -c $LOAD -a $ALGO -s ../resultats_PF/$TEST/$TEST3
