#!/bin/bash

#                                                     TEST WRONG VTK FILES
#=========================================================================
# cas test lineaire elastique anisotrope, def. imposée, algo. par défault

# Inversion spacing et origin : Message erreur sur la sortie standard - OK
TEST="vtklog1"
mkdir resultats/$TEST
mpirun -np 6 ../amitex_fftp  -nm ../microstructures/wrong_vtk/bicouche4_inverse_spacing_origin.vtk -m ../materiaux/mat_linaniso.xml -c ../chargements/char_def_imp.xml -a ../algorithmes/algo_default.xml -s resultats/$TEST/$TEST

# Dimensions et nombre de cellules incompatibles : Message erreur sur la sortie standard - OK
TEST="vtklog2"
mkdir resultats/$TEST
mpirun -np 6 ../amitex_fftp  -nm ../microstructures/wrong_vtk/bicouche4_incompatible_dimensionANDcell_data1.vtk -m ../materiaux/mat_linaniso.xml -c ../chargements/char_def_imp.xml -a ../algorithmes/algo_default.xml -s resultats/$TEST/$TEST

# Dimensions et nombre de cellules incompatibles : Message erreur sur la sortie standard - OK
TEST="vtklog3"
mkdir resultats/$TEST
mpirun -np 6 ../amitex_fftp  -nm ../microstructures/wrong_vtk/bicouche4_incompatible_dimensionANDcell_data2.vtk -m ../materiaux/mat_linaniso.xml -c ../chargements/char_def_imp.xml -a ../algorithmes/algo_default.xml -s resultats/$TEST/$TEST


#                                                     TEST WRONG MAT FILES
#=========================================================================
# Mauvais nom de librairie .so : message d'erreur explicite sur sortie standard - OK
TEST="vtklog4"
mkdir resultats/$TEST
cd ..
mpirun ./amitex_fftp  -nz microstructures/arlequin/arlequin_N3_r5.vtk -m materiaux/err1_mat_lin_diffusion.xml -c chargements/char_diffusion_flux.xml -a algorithmes/algo_default_diffusion.xml -s check_error_messages/resultats/$TEST/$TEST
cd check_error_messages


# Faute de frappe dans un champ xml : l'erreur renvoyée par FoX n'est pas évidente à détecter - BOF  (MAIS DIFFICILE DE FAIRE MIEUX...)
TEST="vtklog5"
mkdir resultats/$TEST
cd ..
mpirun ./amitex_fftp  -nz microstructures/arlequin/arlequin_N3_r5.vtk -m materiaux/err2_mat_lin_diffusion.xml -c chargements/char_diffusion_flux.xml -a algorithmes/algo_default_diffusion.xml -s check_error_messages/resultats/$TEST/$TEST
cd check_error_messages




