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

# version user_Kochman2 : use of the STANDARD ALGORITHM 
#----------------------------------------------------------------------------------------

AMITEX_PF=$DIR/libAmitex/src/user_Kochman2/amitex_fftp

# transformation martensitique 2D [J. Kochmann et al. / comput. Methods Appl. Mech. Engrg. 305 (2016) 89-110]
TEST="martensite_Kochmann2016"
TEST2="martensite_Kochmann2016_2"
mkdir $DIR/resultats_PF/$TEST
$MPIRUNalias $AMITEX_PF  -nz microstructures/zone_${TEST}.vtk -m material/mat_${TEST}.xml -c loading/load_${TEST}.xml -a algorithm/algo_${TEST2}.xml -s ../resultats_PF/$TEST/$TEST2


TEST="martensite_Kochmann2016_polygrain"
TEST2="martensite_Kochmann2016_polygrain_2"
mkdir $DIR/resultats_PF/$TEST
$MPIRUNalias $AMITEX_PF  -nz microstructures/zone_${TEST}.vtk -m material/mat_${TEST}.xml -c loading/load_${TEST}.xml -a algorithm/algo_${TEST2}.xml -s ../resultats_PF/$TEST/$TEST2


