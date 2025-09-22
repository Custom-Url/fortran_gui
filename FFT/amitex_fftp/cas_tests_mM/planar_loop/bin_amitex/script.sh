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

MATE="materials/mat_elasiso_eigs_1.xml"
CHAR="loadings/char_traction.xml"
ALGO="algorithms/algo_noacv_user.xml"

#Simulation with 1 zone per voxel
#---------------------------------
rm -rf results/*
touch results/empty
mpirun -np $Namitex $amitex -NX $NX -NY $NX -NZ $NX -m $MATE -c $CHAR -a $ALGO -s results/res

# simulation with a unique zone
#------------------------------
#MATEVTK="microstructures/geom_box_32.vtk"
#MATEVTK="microstructures/geom_box_256.vtk"
#rm -rf results2/*
#mpirun -np $Namitex $amitex -nm $MATEVTK -m $MATE -c $CHAR -a $ALGO -s results2/res


