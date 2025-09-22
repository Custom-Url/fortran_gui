#!/bin/bash
#module purge
source /home/lgelebart/env_amitex_gcc810-openmpi310.sh
module unload tfel
module load tfel/dev/gcc-8.1.0

#pour utilisation amitex perso ou non
#module load amitex/master/gcc-8.1.0
source /home/lgelebart/amitex_fftp/env_amitex.sh

# 1 voxel
ZVTK="../../../cas_tests/microstructures/homogene/homogene_1.vtk"

MATE="mate/mat_test_charpagne_1agregat.xml"
CHAR="load/traction_100_3D.xml"
ALGO="algo/algo_default_GD.xml"

rm -rf res1
mkdir res1

mpirun -n 1 amitex_fftp -nz $ZVTK -m $MATE -c $CHAR -a $ALGO -s res1/res

# agregat 3x3x3
ZVTK="../../../cas_tests/microstructures/arlequin/arlequin_N3_r1.vtk"
#ZVTK=" /usr/local/AMITEX/AMITEX_FFT/MICROSTRUCTURES/arlequin_N9_r5.vtk"
rm -rf res3
mkdir res3
mpirun  -n 4 amitex_fftp -nz $ZVTK -m $MATE -c $CHAR -a $ALGO -s res3/res

#exit 0

# agregat 10x10x10
ZVTK="../../../cas_tests/microstructures/arlequin/arlequin_N10_r1.vtk"

rm -rf res10
mkdir res10
mpirun amitex_fftp -nz $ZVTK -m $MATE -c $CHAR -a $ALGO -s res10/res

# agregat 9x9x9
ZVTK="../../../cas_tests/microstructures/arlequin/arlequin_N9_r1.vtk"

rm -rf res9
mkdir res9
mpirun amitex_fftp -nz $ZVTK -m $MATE -c $CHAR -a $ALGO -s res9/res



