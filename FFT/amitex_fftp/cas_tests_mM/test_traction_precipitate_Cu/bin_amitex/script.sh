#!/bin/bash

#--- VARIABLES D'ENVIRONNEMENT POUR AMITEX
source ../env_amitex.sh

#non cubic grid because the mM grid must be a multiple of 8 (because of the CFC symmetry) 
#    and cannot be cubic to avoid annihilation
NX="115"
NY="127"
NZ="154"
# DO NOT CHANGE these values without changing half_thickness in DCM/DCM_ContCu

Namitex="1"

# RESOLUTION PB MPI_SPAWN AVEC OPENMPI(PSM2) et OMNIPATH 
#------------------------------------------------------- 
# Test si PBS_NODEFILE n'est pas vide (-n) :
#	si ok : on utilise qsub -> option ofi
#	sinon : on lance sur un noeud -> pas d'option differente
if [ -n "$PBS_NODEFILE" ]; then
   export OMPI_MCA_mtl=ofi  
fi 

amitex="../../../libAmitex/src/user_ddd_mm2/amitex_fftp"

MATE="materials/mat_elasiso_eigs_1.xml"
CHAR="loadings/char_traction.xml"
ALGO="algorithms/algo_noacv_user.xml"

#Simulation with 1 zone per voxel
#---------------------------------
rm -rf results/*
touch results/empty
mpirun -np $Namitex $amitex -NX $NX -NY $NY -NZ $NZ -m $MATE -c $CHAR -a $ALGO -s results/res


