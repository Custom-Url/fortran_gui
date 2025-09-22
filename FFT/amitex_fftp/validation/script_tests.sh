#!/bin/bash
if [ ! $(basename $(pwd)) = "validation" ]; then
   echo "ERROR : script_tests.sh must be launched in its own directory 'validation'"
   exit 0
fi
#-----------------------------------------------------------------------

DATE='/bin/date'
BEFORE=$($DATE +'%s')

#----------------------- DEFINITION DES VARIABLES D'ENVIRONNEMENT POUR AMITEX
source ../env_amitex.sh

./clean_results.sh
rm tests.log
rm time.log

cd ../cas_tests


MPIRUNalias='mpirun'
MPIRUNalias1='mpirun -np 1'
MPIRUNalias4='mpirun -np 4'

#-- MPI/OPENMP sur Maldives (exploratoire)
#export OMP_NUM_THREADS=6
#MPIRUNalias='mpirun -np 2 -ppn 6 -print-rank-map'
#MPIRUNalias1='mpirun -np 1 -ppn 1 -print-rank-map'

#-- OPTIMISATION openMPI - OMNIPATH
#MPIRUNalias='mpirun -mca pml cm -mca mtl psm2'
#MPIRUNalias1='mpirun -np 1 -mca pml cm -mca mtl psm2'

#-- AU CCRT (cobalt)
#MPIRUNalias='ccc_mprun'
#MPIRUNalias1='ccc_mprun -n 1'
#MPIRUNalias4='ccc_mprun -n 4'


AMITEX='amitex_fftp'

# Example to chose the decomposition p_row,p_col
#if user_2decomp.txt is not given or wrong => use default decomposition
#AMITEX='../libAmitex/bin/amitex_fftp -2decomp user'
#
# echo "4"  > user_2decomp.txt
# echo "7" >> user_2decomp.txt

# Example to chose a decomposition with low p_row
#AMITEX='../libAmitex/bin/amitex_fftp -2decomp low_p_row'

Nsim=0

# cas test : bi-couche croix (couches materiau et zone ont des orientations differentes)
#            4proc : chaque proc a un materiau absent ET une zone absente
TEST="bilayer_2mat_2zone_xmla"
$MPIRUNalias4 $AMITEX -nm microstructures/multilayers/2layersM.vtk -nz microstructures/multilayers/2layersZ.vtk -m materiaux/mat_2mat_2zone_xml.xml -c chargements/char_traction_1.xml -a algorithmes/algo_default.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log
#exit 0
# cas test : un materiau continu par voxel (sans rentrer de fichier zone.vtk)
TEST="linelas_continuous"
$MPIRUNalias $AMITEX  -m materiaux/mat_continuous_R10.xml -c chargements/char_traction_1.xml -a algorithmes/algo_default.xml -NX 10 -NY 10 -NZ 10 -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log

# cas test utilisation de voxels composite "laminate" : polycrystal + elasticite isotrope heterogene
TEST="laminate_polyx_elasiso"
$MPIRUNalias $AMITEX  -nz microstructures/voxcomp/polyX_27G_R21.vtk -m materiaux/mat_polyx_laminate_elasiso.xml -c chargements/char_traction_1.xml -a algorithmes/algo_default_voxcomplaminate.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log 

# cas test utilisation de voxels composite "reuss" : polycrystal + elasticite isotrope heterogene
TEST="reuss_polyx_elasiso"
$MPIRUNalias $AMITEX  -nz microstructures/voxcomp/polyX_27G_R21.vtk -m materiaux/mat_polyx_reuss_elasiso.xml -c chargements/char_traction_1.xml -a algorithmes/algo_default.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log

# meme cas test que le precedant MAIS sans voxels composites
TEST="polyx_elasiso"
$MPIRUNalias $AMITEX  -nz microstructures/voxcomp/polyX_27G_R21.vtk -m materiaux/mat_polyx_elasiso.xml -c chargements/char_traction_1.xml -a algorithmes/algo_default.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log  

# cas test diffusion lineaire flux impose sur maillage arlequin N3_5
TEST="diffusion_fluximp"
$MPIRUNalias $AMITEX  -nz microstructures/arlequin/arlequin_N3_r5.vtk -m materiaux/mat_lin_diffusion.xml -c chargements/char_diffusion_flux.xml -a algorithmes/algo_default_diffusion.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log  

# cas test diffusion lineaire gradQ impose sur maillage arlequin N3_5
TEST="diffusion_gradQimp"
$MPIRUNalias $AMITEX  -nz microstructures/arlequin/arlequin_N3_r5.vtk -m materiaux/mat_lin_diffusion.xml -c chargements/char_diffusion_gradQ.xml -a algorithmes/algo_default_diffusion.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log  

# cas test lineaire deformation imposee sur un maillage 65x65x65 avec la premiere methode (sans acceleration de convergence ni filtre)
TEST="lin_elas_def_imp_65_old"
$MPIRUNalias $AMITEX  -nz microstructures/al2p5_65/mat_al2p5_65.vtk -m materiaux/mat_lin_2.xml -c chargements/char_def_imp.xml -a algorithmes/algo_old.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log  

# cas test lineaire deformation imposee sur un maillage 64x64x64 avec la premiere methode (sans acceleration de convergence ni filtre)
TEST="lin_elas_def_imp_64_old"
$MPIRUNalias $AMITEX  -nz microstructures/al2p5_64/mat_al2p5_64.vtk -m materiaux/mat_lin_2.xml -c chargements/char_def_imp.xml -a algorithmes/algo_old.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log  

# cas test lineaire deformation imposee sur un maillage 65x65x65 avec acceleration de convergence et un filtre hexaedrique
TEST="lin_elas_def_imp_65"
$MPIRUNalias $AMITEX  -nz microstructures/al2p5_65/mat_al2p5_65.vtk -m materiaux/mat_lin_2.xml -c chargements/char_def_imp.xml -a algorithmes/algo_default.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log  

# cas test lineaire deformation imposee sur un maillage 64x64x64 avec acceleration de convergence et un filtre hexaedrique
TEST="lin_elas_def_imp_64"
$MPIRUNalias $AMITEX  -nz microstructures/al2p5_64/mat_al2p5_64.vtk -m materiaux/mat_lin_2.xml -c chargements/char_def_imp.xml -a algorithmes/algo_default.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log  

# cas test lineaire deformation imposee sur un maillage 65x65x65 avec un schema octaedrique
TEST="lin_elas_def_imp_65_octa"
$MPIRUNalias $AMITEX  -nz microstructures/al2p5_65/mat_al2p5_65.vtk -m materiaux/mat_lin_2.xml -c chargements/char_def_imp.xml -a algorithmes/algo_octa.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log  

# cas test lineaire deformation imposee sur un maillage 64x64x64 avec un schema octaedrique
TEST="lin_elas_def_imp_64_octa"
$MPIRUNalias $AMITEX  -nz microstructures/al2p5_64/mat_al2p5_64.vtk -m materiaux/mat_lin_2.xml -c chargements/char_def_imp.xml -a algorithmes/algo_octa.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log  

# cas test lineaire deformation imposee sur un maillage 65x65x65 sans acceleration de convergence
TEST="lin_elas_def_imp_65_noAC"
$MPIRUNalias $AMITEX  -nz microstructures/al2p5_65/mat_al2p5_65.vtk -m materiaux/mat_lin_2.xml -c chargements/char_def_imp.xml -a algorithmes/algo_no_AC.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log  

# cas test lineaire deformation imposee sur un maillage 64x64x64 sans acceleration de convergence
TEST="lin_elas_def_imp_64_noAC"
$MPIRUNalias $AMITEX  -nz microstructures/al2p5_64/mat_al2p5_64.vtk -m materiaux/mat_lin_2.xml -c chargements/char_def_imp.xml -a algorithmes/algo_no_AC.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log  

# Test d'une loi thermoelastique deformation imposee sur un maillage 64x64x64 avec acceleration de convergence et un filtre hexaedrique
TEST="thermoelas_def_imp_64"
$MPIRUNalias $AMITEX  -nz microstructures/al2p5_64/mat_al2p5_64.vtk -m materiaux/mat_thermo.xml -c chargements/char_def_imp_thermo.xml -a algorithmes/algo_default.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log  

# Test d'une loi thermoelastique deformation imposee nulle sur un materiau bicouche (2 voxels)
TEST="bicouche_thermoelas"
$MPIRUNalias1 $AMITEX  -nz microstructures/2couches1.vtk -m materiaux/mat_thermo_bicouche.xml -c chargements/char_thermo_bicouche.xml -a algorithmes/algo_default.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log  

# Test d'une loi thermoelastique, temperature en parametre exterieur, deformation imposee nulle sur un materiau bicouche (2 voxels)
TEST="bicouche_paramextelas"
TEST2="bicouche_thermoelas"
$MPIRUNalias1 $AMITEX  -nz microstructures/2couches1.vtk -m materiaux/mat_paramext_bicouche.xml -c chargements/char_paramext_bicouche.xml -a algorithmes/algo_default.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST2/reference/$TEST2 ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log  

# cas test beton deformation imposee sur un maillage 65x65x65 avec la premiere methode (sans acceleration de convergence ni filtre)
TEST="beton_relax_65_old"
$MPIRUNalias $AMITEX  -nm microstructures/al2p5_65/mat_al2p5_65.vtk -nz microstructures/al2p5_65/zone_al2p5_65.vtk -m materiaux/mat_beton.xml -c chargements/char_def_imp_relax.xml -a algorithmes/algo_old.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log  
#exit 0
# cas beton_relax_65 + variable interne initilisee par fichier VTK + acceleration de convergence et un filtre hexaedrique
# juste ecriture vtk de la variable interne initilaisee par fichier vtk
$MPIRUNalias $AMITEX  -nm microstructures/al2p5_65/mat_al2p5_65.vtk -nz microstructures/al2p5_65/zone_al2p5_65.vtk -m materiaux/mat_beton_intvarVTK.xml -c chargements/char_def_imp_relax.xml -a algorithmes/algo_default.xml -s ../resultats/$TEST/$TEST -varint2print 20

# cas beton_relax_65 avec acceleration de convergence et un filtre hexaedrique
TEST="beton_relax_65"
$MPIRUNalias $AMITEX  -nm microstructures/al2p5_65/mat_al2p5_65.vtk -nz microstructures/al2p5_65/zone_al2p5_65.vtk -m materiaux/mat_beton.xml -c chargements/char_def_imp_relax.xml -a algorithmes/algo_default.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log  

# cas beton_decharge_65 avec acceleration de convergence et un filtre hexaedrique
TEST="beton_decharge_65"
$MPIRUNalias $AMITEX  -nm microstructures/al2p5_65/mat_al2p5_65.vtk -nz microstructures/al2p5_65/zone_al2p5_65.vtk -m materiaux/mat_beton.xml -c chargements/char_decharge.xml -a algorithmes/algo_default.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log  

# cas fissure avec acceleration de convergence et un filtre hexaedrique
TEST="essai_fissure"
$MPIRUNalias $AMITEX -nz microstructures/fissure/fissure_64.vtk -m materiaux/mat_fissure.xml -c chargements/char_traction_zz.xml -a algorithmes/algo_default.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log  

# cas fissure avec acceleration de convergence et un filtre hexaedrique sur 1 processeur
TEST_1p="essai_fissure_1proc"
$MPIRUNalias1 $AMITEX -nz microstructures/fissure/fissure_64.vtk -m materiaux/mat_fissure.xml -c chargements/char_traction_zz.xml -a algorithmes/algo_default.xml -s ../resultats/$TEST_1p/$TEST_1p

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST_1p/$TEST_1p -n $TEST_1p -s ../validation/tests.log
let Nsim++ ;echo $TEST_1p >> ../validation/tests.log  

# cas fissure avec acceleration de convergence,filtre hexaedrique et modACV=2
TEST="essai_fissure_modACV"
$MPIRUNalias $AMITEX -nz microstructures/fissure/fissure_64.vtk -m materiaux/mat_fissure.xml -c chargements/char_traction_zz.xml -a algorithmes/algo_default_modACV.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log  

# cas beton_fluage_64 avec acceleration de convergence et un filtre hexaedrique
TEST="beton_fluage_64"
$MPIRUNalias $AMITEX  -nm microstructures/al2p5_64/mat_al2p5_64.vtk -nz microstructures/al2p5_64/zone_al2p5_64.vtk -m materiaux/mat_beton.xml -c chargements/char_fluage.xml -a algorithmes/algo_default.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log  

# meme simu que precedemment avec interruption sur la contrainte axiale de fluage issue du calcul precedent
TEST="beton_fluage_64_interrupt"
LD_LIBRARY_PATH0=$LD_LIBRARY_PATH
LD_LIBRARY_PATH=/u/q/dg765/amitex_fftp/cas_tests/lib_user_functions:$LD_LIBRARY_PATH

$MPIRUNalias $AMITEX  -nm microstructures/al2p5_64/mat_al2p5_64.vtk -nz microstructures/al2p5_64/zone_al2p5_64.vtk -m materiaux/mat_beton.xml -c chargements/char_fluage2.xml -a algorithmes/algo_default.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log  
LD_LIBRARY_PATH=$LD_LIBRARY_PATH0

# cas polyxCC avec acceleration de convergence et un filtre hexaedrique
TEST="polyxCC_def_imp"
$MPIRUNalias $AMITEX -nz microstructures/27G_65/zone_27G_65.vtk -m materiaux/mat_polyx_bin.xml -c chargements/char_def_imp_polyx.xml -a algorithmes/algo_default.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log  

# cas elastique en grandes transformations deformation imposee sur un maillage 65x65x65 avec acceleration de convergence et un filtre hexaedrique
TEST="elas_def_imp_65_GD"
$MPIRUNalias $AMITEX  -nz microstructures/al2p5_65/mat_al2p5_65.vtk -m materiaux/mat_lin_GD.xml -c chargements/char_def_imp_GD.xml -a algorithmes/algo_default_GD.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log  

# cas traction d'un polycristal en grandes transformations avec acceleration de convergence et un filtre hexaedrique
TEST="polyxGD_traction"
$MPIRUNalias $AMITEX  -nz microstructures/arlequin/arlequin_N2_r16.vtk -m materiaux/mat_polyx_GD.xml -c chargements/char_def_imp_polyx_GD.xml -a algorithmes/algo_default_GD.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log 

# cas traction d'un polycristal avec interruption (Attention, verifier que LD_LIBRARY_PATH a ete correctement ajuste ci-dessous)
TEST="polyxGD_traction_interrupt"
LD_LIBRARY_PATH0=$LD_LIBRARY_PATH
LD_LIBRARY_PATH=/u/q/dg765/amitex_fftp/cas_tests/lib_user_functions:$LD_LIBRARY_PATH

$MPIRUNalias $AMITEX  -nz microstructures/arlequin/arlequin_N2_r16.vtk -m materiaux/mat_polyx_GD.xml -c chargements/char_def_imp_polyx_GD_interrupt.xml -a algorithmes/algo_default_GD.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log  

LD_LIBRARY_PATH=$LD_LIBRARY_PATH0 

# cas traction d'un polycristal en grandes transformations (ACV/filtre hexa) + pilotage def33 a triaxialite (PK1) impoee
TEST="polyxGD_triax2p5_PK1"
$MPIRUNalias $AMITEX  -nz microstructures/arlequin/arlequin_N2_r16.vtk -m materiaux/mat_polyx_GD.xml -c chargements/char_triax2p5_PK1.xml -a algorithmes/algo_default_GD.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log  

# idem cas precedent MAIS : 1/ triaxialite impose sur CAUCHY 2/ critere 1e-3 sur Smacro
TEST="polyxGD_triax2p5_cauchy"
$MPIRUNalias $AMITEX  -nz microstructures/arlequin/arlequin_N2_r16.vtk -m materiaux/mat_polyx_GD.xml -c chargements/char_triax2p5_cauchy.xml -a algorithmes/algo_Smacro_GD.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log  

# cas gurtin (non local) hpp acceleration de convergence et prise en compte du gradient de u complet
# non-local utilisateur
AMITEX2='../libAmitex/src/user_gurtinSSnloc3/amitex_fftp'
TEST="gurtin_user_hpp_nsym"
$MPIRUNalias1 $AMITEX2 -nz microstructures/monocristal/MonoX_bizone_R21.vtk -m materiaux/mat_gurtin_user_hpp.xml -c chargements/char_traction_gurtin.xml -a algorithmes/algo_gurtin_hpp_nsym.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log  

# cas test : torsion poutre beam (inclusion cubique) 
TEST="torsion_beam"
$MPIRUNalias $AMITEX  -nm microstructures/beam/cube44_cylinder_mate.vtk -nz microstructures/beam/cube44_cylinder_zone.vtk -m materiaux/mat_renfort2_i_m0.xml -c chargements/char_defimp_torsion_beam.xml -a algorithmes/algo_default.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log

# cas test : flexion poutre beam (inclusion cubique) 
TEST="flexion_beam"
$MPIRUNalias $AMITEX  -nm microstructures/beam/cube44_cylinder_mate.vtk -nz microstructures/beam/cube44_cylinder_zone.vtk -m materiaux/mat_renfort2_i_m0.xml -c chargements/char_defimp_flexion_beam.xml -a algorithmes/algo_default.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log

# cas test : flexion plaque laminate (article IJSS-2008, Nguyen, Sab, Bonnet)
TEST="flex_laminate_NSB"
$MPIRUNalias $AMITEX  -nm microstructures/multilayers/laminate_NSB_n64_mate.vtk -nz microstructures/multilayers/laminate_NSB_n64_zone.vtk -m materiaux/mat_laminate_NSB.xml -c chargements/char_defimp_flexion_NSB.xml -a algorithmes/algo_default.xml -s ../resultats/$TEST/$TEST

../validation/compare ../resultats/$TEST/reference/$TEST ../resultats/$TEST/$TEST -n $TEST -s ../validation/tests.log
let Nsim++ ;echo $TEST >> ../validation/tests.log

#
AFTER=$($DATE +'%s')
ELAPSED=$(($AFTER - $BEFORE))

# TEMPS TOTAL
#------------
echo "Temps d'execution des tests (s) = " $ELAPSED > ../validation/time.log

# TEST FINAL QUE LE NOMBRE DE lignes OK de tests.log correspond au nombre de simulations
#---------------------------------------------------------------------------------------
NbOK=`grep "OK" ../validation/tests.log | wc -l`

echo "========================================================================" >> ../validation/tests.log
if [ $Nsim != $NbOK ]
then 
     echo "                                                             WARNING !!!" >> ../validation/tests.log
     echo "                                                Total "$Nsim" / Succeed "$NbOK >> ../validation/tests.log
else
     echo "                                                                GOOD !!!" >> ../validation/tests.log
     echo "                                                Total "$Nsim" / Succeed "$NbOK >> ../validation/tests.log
fi

echo "" >>../validation/tests.log
echo "   REMARKS : * Test 'gurtin_user_hpp_nsym' may fail if mfront is not installed" >>../validation/tests.log  
echo "" >>../validation/tests.log
	
echo "             * Tests 'flexion_beam' and 'torsion_beam' may fail " >>../validation/tests.log  
echo "                      BUT the error err_sig should remain very small  " >>../validation/tests.log  
echo "" >>../validation/tests.log
	

exit 0



