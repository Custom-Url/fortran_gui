#!/bin/sh

module purge
module load gcc/8.1.0
module load tfel/dev/gcc-8.1.0

# SURTOUT PAS!!
# l'ajout du module load openmpi pourrit l'execution de la loi!!
#source ~/env_amitex_gcc810-openmpi310.sh


LOI="FCC_charpagne_v2.mfront"

mkdir cast3m21
mkdir amitex

cp $LOI cast3m21
cd cast3m21
mfront --obuild --interface=castem21 $LOI
cd ..

cp $LOI amitex
cd amitex
mfront --obuild --interface=umat $LOI
#mfront --obuild=level2 --interface=umat $LOI


