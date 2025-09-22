#!/bin/bash

# Phase Field J. Boisse / adaptation LG

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

# version user_Kochman3 : use of the STANDARD ALGORITHM with NLOC (the best way to do in that case)
#----------------------------------------------------------------------------------------

AMITEX_PF=$DIR/libAmitex/src/user_Kochman3/amitex_fftp

# transformation martensitique 2D [J. Kochmann et al. / comput. Methods Appl. Mech. Engrg. 305 (2016) 89-110]
TEST="martensite_Kochmann2016"
TEST3="martensite_Kochmann2016_3"
mkdir $DIR/resultats_PF/$TEST
$MPIRUNalias $AMITEX_PF  -nz microstructures/zone_${TEST}.vtk -m material/mat_${TEST3}.xml -c loading/load_${TEST}.xml -a algorithm/algo_${TEST3}.xml -s ../resultats_PF/$TEST/$TEST3


TEST="martensite_Kochmann2016_polygrain"
TEST3a=$TEST3 #on censerve le meme algo
TEST3="martensite_Kochmann2016_polygrain_3"
mkdir $DIR/resultats_PF/$TEST
$MPIRUNalias $AMITEX_PF  -nz microstructures/zone_${TEST}.vtk -m material/mat_${TEST3}.xml -c loading/load_${TEST}.xml -a algorithm/algo_${TEST3a}.xml -s ../resultats_PF/$TEST/$TEST3


