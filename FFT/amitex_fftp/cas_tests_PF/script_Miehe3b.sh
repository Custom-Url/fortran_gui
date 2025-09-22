#!/bin/bash

# Phase Field implemented by Y. Chen

#--VARIABLES D'ENVIRONNEMENT POUR AMITEX
source ../env_amitex.sh

MPIRUNalias='mpirun'
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
ALGO="algorithm/algo_Miehe3b.xml"


# ====================================
# Single Edge Notched Tensile specimen

TEST="SENT3b"
mkdir ../resultats_PF/$TEST

rm ../resultats_PF/$TEST/SENT*
$MPIRUNalias $AMITEX_PF -nm microstructures/SENT.vtk -m material/mat_SENT3.xml -c loading/load_SENT.xml -a $ALGO -s ../resultats_PF/$TEST/SENT


# ==========================================================================
# Single Edge Notched Tensile specimen with inhomogeneous gc or lc 
TEST="SENT_VARgclc3b"
mkdir ../resultats_PF/$TEST

rm results/$TEST/SNmulti_FIXgc_VARlc_2*
$MPIRUNalias $AMITEX_PF -nm microstructures/SENT_VARgclc.vtk -m material/mat_SENT_FIXgc_VARlc3.xml -c loading/load_SENT_VARgclc.xml -a $ALGO -s ../resultats_PF/$TEST/SENT_FIXgc_VARlc_2

rm results/$TEST/SNmulti_FIXlc_VARgc_2*
$MPIRUNalias $AMITEX_PF -nm microstructures/SENT_VARgclc.vtk -m material/mat_SENT_FIXlc_VARgc3.xml -c loading/load_SENT_VARgclc.xml -a $ALGO -s ../resultats_PF/$TEST/SENT_FIXlc_VARgc_2


# ====================================================================
# Cracked body under shear loads (xy=Mode II, xz/yz=ModeIII)
TEST="shearCrack3b"
mkdir ../resultats_PF/$TEST

rm ../resultats_PF/$TEST/SC_xy*
$MPIRUNalias $AMITEX_PF -nm microstructures/shearCrack.vtk -m material/mat_shearCrack3.xml -c loading/load_shearCrack_xy.xml -a $ALGO -s ../resultats_PF/$TEST/SC_xy

rm ../resultats_PF/$TEST/SC_xz*
$MPIRUNalias $AMITEX_PF -nm microstructures/shearCrack.vtk -m material/mat_shearCrack3.xml -c loading/load_shearCrack_xz.xml -a $ALGO -s ../resultats_PF/$TEST/SC_xz

rm ../resultats_PF/$TEST/SC_yz*
$MPIRUNalias $AMITEX_PF -nm microstructures/shearCrack.vtk -m material/mat_shearCrack3.xml -c loading/load_shearCrack_yz.xml -a $ALGO -s ../resultats_PF/$TEST/SC_yz


# ===========================================================================================
# Double Edge Notched Tensile specimen with symmetric or asymmetric notch positions
# NOTE: These simulations are larger than the others
TEST="DENT3b"
mkdir ../resultats_PF/$TEST

rm ../resultats_PF/$TEST/DENTasym*
$MPIRUNalias $AMITEX_PF -nm microstructures/DENTasym.vtk -m material/mat_DENT3.xml -c loading/load_DENT.xml -a $ALGO -s ../resultats_PF/$TEST/DENTasym

rm ../resultats_PF/$TEST/DENTsym*
$MPIRUNalias $AMITEX_PF -nm microstructures/DENTsym.vtk -m material/mat_DENT3.xml -c loading/load_DENT.xml -a $ALGO -s ../resultats_PF/$TEST/DENTsym


# ============================================================
# Bi-material tensile specimen (crack branching & coalescence)
# NOTE: This simulation is larger than the others
TEST="bimat3b"
mkdir ../resultats_PF/$TEST

rm ../resultats_PF/$TEST/bimat*
$MPIRUNalias $AMITEX_PF -nm microstructures/bimat.vtk -m material/mat_bimat3.xml -c loading/load_bimat.xml -a $ALGO -s ../resultats_PF/$TEST/bimat



