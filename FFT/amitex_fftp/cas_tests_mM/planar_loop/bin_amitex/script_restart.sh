#!/bin/bash

#--- VARIABLES D'ENVIRONNEMENT POUR AMITEX
source ../env_amitex.sh

if [ -z ${Namitex+x} ]; then Namitex="10"; fi
if [ -z ${NX+x} ]; then NX="32"; fi
echo "Namitex='$Namitex'"
echo "NX='$NX'"

# RESOLUTION PB MPI_SPAWN AVEC OPENMPI(PSM2) et OMNIPATH 
#------------------------------------------------------- 
# Test si PBS_NODEFILE n'est pas vide (-n) :
#	si ok : on utilise qsub -> option ofi
#	sinon : on lance sur un noeud -> pas d'option differente
if [ -n "$PBS_NODEFILE" ]; then
   export OMPI_MCA_mtl=ofi  
fi 

amitex="../../../libAmitex/src/user_ddd_mm2/amitex_fftp"

MATE="materials/mat_elasiso_eigs_1_restart.xml"
CHAR="loadings/char_traction.xml"
ALGO="algorithms/algo_noacv_user_restart.xml"

#Simulation with 1 zone per voxel
#---------------------------------
rm -rf results_restart/*
touch results_restart/empty
mpirun -np $Namitex $amitex -NX $NX -NY $NX -NZ $NX -m $MATE -c $CHAR -a $ALGO -s results_restart/res


