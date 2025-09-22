#!/bin/bash

# Amitex Multi-Echelles - 1er test de lancement de 2 calculs

#--VARIABLES D'ENVIRONNEMENT POUR AMITEX
source ../env_amitex.sh

MPIRUNalias='mpirun'
MPIRUNalias1='mpirun -np 1'
AMITEX='../libAmitex/src/user_biscale2/amitex_fftp'

mkdir tmp
mkdir cmd_files
#                                                CALCUL "1 ECHELLE" GLOB
#-----------------------------------------------------------------------
ALGO="algorithmes/algo_default.xml"
CHAR="../cas_tests/chargements/char_def_imp.xml"
MATE="../cas_tests/materiaux/mat_lin_2.xml"
ZVTK="../cas_tests/microstructures/al2p5_64/mat_al2p5_64.vtk"

#$MPIRUNalias $AMITEX -nz $ZVTK -a $ALGO -m $MATE -c $CHAR -s tmp/glob0

#                                         CALCUL "2 ECHELLES" - LINEAIRE
#                                                       UNE SEULE GRILLE
#-----------------------------------------------------------------------
#  creation fichier de commande echelle loc 
#        (pour faire simple ici on reprend les donnees du calcul glob) 
CMDFILE="cmd_files/cmd_loc.in"
echo "&CMD" > $CMDFILE
echo "fic_numZ=\""$ZVTK"\"" >> $CMDFILE
echo "fic_mat=\""$MATE"\"" >> $CMDFILE
echo "fic_char=\""$CHAR"\"" >> $CMDFILE
echo "fic_algo=\""$ALGO"\"" >> $CMDFILE
echo "fic_vtk=\"tmp/loc\"" >> $CMDFILE
echo "/">> $CMDFILE

#ALGO2="algorithmes/algo_default_biscale2.xml"
#$MPIRUNalias $AMITEX -nz $ZVTK -a $ALGO2 -m $MATE -c $CHAR -s tmp/glob1

#                                         CALCUL "2 ECHELLES" - LINEAIRE
#                                                           DEUX GRILLES
#-----------------------------------------------------------------------
#  creation fichier de commande echelle loc 
#        (pour faire simple ici on reprend les donnees du calcul glob) 
ZVTK2="../cas_tests/microstructures/al2p5_65/mat_al2p5_65.vtk"
CMDFILE="cmd_files/cmd_loc.in"
echo "&CMD" > $CMDFILE
echo "fic_numZ=\""$ZVTK2"\"" >> $CMDFILE
echo "fic_mat=\""$MATE"\"" >> $CMDFILE
echo "fic_char=\""$CHAR"\"" >> $CMDFILE
echo "fic_algo=\""$ALGO"\"" >> $CMDFILE
echo "fic_vtk=\"tmp/loc\"" >> $CMDFILE
echo "/">> $CMDFILE

#ALGO2="algorithmes/algo_default_biscale2.xml"
#$MPIRUNalias $AMITEX -nz $ZVTK -a $ALGO2 -m $MATE -c $CHAR -s tmp/glob1

#exit 0

#                   CALCUL "2 ECHELLES" - LINEAIRE - CHARGEMENT COMPLEXE
#-----------------------------------------------------------------------
#  meme chose qu'avant mais avec chargement different
CHARCPLX="chargements/char_cplx.xml"

CMDFILE="cmd_files/cmd_loc.in"
echo "&CMD" > $CMDFILE
echo "fic_numZ=\""$ZVTK"\"" >> $CMDFILE
#echo "fic_numM=\"\"" >> $CMDFILE
echo "fic_mat=\""$MATE"\"" >> $CMDFILE
echo "fic_char=\""$CHARCPLX"\"" >> $CMDFILE
echo "fic_algo=\""$ALGO"\"" >> $CMDFILE
echo "fic_vtk=\"tmp/loc\"" >> $CMDFILE
echo "/">> $CMDFILE

#ALGO2="algorithmes/algo_default_biscale2.xml"
#$MPIRUNalias $AMITEX -nz $ZVTK -a $ALGO2 -m $MATE -c $CHARCPLX -s tmp/glob1



#                               CALCUL "2 ECHELLES" - TEST PB CHARGEMENT
#-----------------------------------------------------------------------
#  meme chose qu'avant mais avec chargement different loc/glob
CHAR2="../cas_tests/chargements/char_def_imp_polyx.xml"

CMDFILE="cmd_files/cmd_loc.in"
echo "&CMD" > $CMDFILE
echo "fic_numZ=\""$ZVTK"\"" >> $CMDFILE
#echo "fic_numM=\"\"" >> $CMDFILE
echo "fic_mat=\""$MATE"\"" >> $CMDFILE
echo "fic_char=\""$CHAR2"\"" >> $CMDFILE
echo "fic_algo=\""$ALGO"\"" >> $CMDFILE
echo "fic_vtk=\"tmp/loc\"" >> $CMDFILE
echo "/">> $CMDFILE

#ALGO2="algorithmes/algo_default_biscale2.xml"
#$MPIRUNalias $AMITEX -nz $ZVTK -a $ALGO2 -m $MATE -c $CHAR -s tmp/glob2


#                                     CALCUL "2 ECHELLES" - NON-LINEAIRE
#-----------------------------------------------------------------------
#  creation fichier de commande echelle loc 
#        (pour faire simple ici on reprend les donnees du calcul glob) 
ZVTK="../cas_tests/microstructures/al2p5_64/zone_al2p5_64.vtk"
MVTK="../cas_tests/microstructures/al2p5_64/mat_al2p5_64.vtk"
CHAR="../cas_tests/chargements/char_fluage.xml"
MATE="../cas_tests/materiaux/mat_beton.xml"

CMDFILE="cmd_files/cmd_loc.in"
echo "&CMD" > $CMDFILE
echo "fic_numZ=\""$ZVTK"\"" >> $CMDFILE
echo "fic_numM=\""$MVTK"\"" >> $CMDFILE
echo "fic_mat=\""$MATE"\"" >> $CMDFILE
echo "fic_char=\""$CHAR"\"" >> $CMDFILE
echo "fic_algo=\""$ALGO"\"" >> $CMDFILE
echo "fic_vtk=\"tmp/loc\"" >> $CMDFILE
echo "/">> $CMDFILE

#ALGO2="algorithmes/algo_default_biscale2.xml"
#$MPIRUNalias $AMITEX -nz $ZVTK -nm $MVTK -a $ALGO2 -m $MATE -c $CHAR -s tmp/glob3


#                                     CALCUL "2 ECHELLES" - NON-LINEAIRE
#                                                              2 GRILLES
#-----------------------------------------------------------------------
#  creation fichier de commande echelle loc 
ZVTK="../cas_tests/microstructures/al2p5_64/zone_al2p5_64.vtk"
MVTK="../cas_tests/microstructures/al2p5_64/mat_al2p5_64.vtk"
CHAR="../cas_tests/chargements/char_fluage.xml"
MATE="../cas_tests/materiaux/mat_beton.xml"

ZVTK2="../cas_tests/microstructures/al2p5_65/zone_al2p5_65.vtk"
MVTK2="../cas_tests/microstructures/al2p5_65/mat_al2p5_65.vtk"

CMDFILE="cmd_files/cmd_loc.in"
echo "&CMD" > $CMDFILE
echo "fic_numZ=\""$ZVTK2"\"" >> $CMDFILE
echo "fic_numM=\""$MVTK2"\"" >> $CMDFILE
echo "fic_mat=\""$MATE"\"" >> $CMDFILE
echo "fic_char=\""$CHAR"\"" >> $CMDFILE
echo "fic_algo=\""$ALGO"\"" >> $CMDFILE
echo "fic_vtk=\"tmp/loc\"" >> $CMDFILE
echo "/">> $CMDFILE

ALGO2="algorithmes/algo_default_biscale2.xml"
$MPIRUNalias $AMITEX -nz $ZVTK -nm $MVTK -a $ALGO2 -m $MATE -c $CHAR -s tmp/glob3

