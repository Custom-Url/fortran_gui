#!/bin/bash

#--- VARIABLES D'ENVIRONNEMENT POUR AMITEX
source ../env_amitex.sh

# RESOLUTION PB MPI_SPAWN AVEC OPENMPI(PSM2) et OMNIPATH 
#------------------------------------------------------- 
# Test si PBS_NODEFILE n'est pas vide (-n) :
#	si ok : on utilise qsub -> option ofi
#	sinon : on lance sur un noeud -> pas d'option differente

if [ -n "$PBS_NODEFILE" ]; then
   export OMPI_MCA_mtl=ofi  
fi 

# DO NOT CHANGE these values without changing half_thickness in DCM/DCM_ContCu
NX="920"
NY="1016"
NZ="1232"

# ok on 3 nodes with
#NX="805"
#NY="889"
#NZ="1078"

Namitex="80"  # on 3 nodes (3*28-4) => PB
Namitex="108" # on 4 nodes (4*28-4) => OK



amitex="../../../libAmitex/src/user_ddd_mm2/amitex_fftp"

MATE="materiaux/mat_linisodefLVI_1.xml"
CHAR="chargements/char_traction.xml"
ALGO="algorithmes/algo_default_user.xml"
#ALGO="algorithmes/algo_default_user_noDD.xml"

# lancement calcul (attention a conserver "N" procs pour "mM"
rm -rf results/*
mpirun -np $Namitex $amitex -NX $NX -NY $NY -NZ $NZ -m $MATE -c $CHAR -a $ALGO -s results/res


