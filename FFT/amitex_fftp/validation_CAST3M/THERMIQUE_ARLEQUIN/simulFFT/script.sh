
module load amitex/8.1.0/gcc-8.1.0
AMITEX="../../../libAmitex/bin/amitex_fftp"
#AMITEX="amitex_fftp"

# ARLEQUIN N3_3
#mpirun -np 1 $AMITEX  -nz /usr/local/AMITEX/AMITEX_FFT/MICROSTRUCTURES/arlequin_N3_r3.vtk -m mat_lin_diffusion.xml -c char_diffusion.xml -a algo_default_diffusion_noACV.xml -s resultats/r3

# ARLEQUIN N3_9
#mpirun $AMITEX  -nz /usr/local/AMITEX/AMITEX_FFT/MICROSTRUCTURES/arlequin_N3_r9.vtk -m mat_lin_diffusion.xml -c char_diffusion.xml -a algo_default_diffusion_noACV.xml -s resultats/r9

# ARLEQUIN N3_4
#mpirun $AMITEX  -nz /usr/local/AMITEX/AMITEX_FFT/MICROSTRUCTURES/arlequin_N3_r4.vtk -m mat_lin_diffusion.xml -c char_diffusion.xml -a algo_default_diffusion_noACV.xml -s resultats/r4

# ARLEQUIN N3_10
#mpirun $AMITEX  -nz /usr/local/AMITEX/AMITEX_FFT/MICROSTRUCTURES/arlequin_N3_r10.vtk -m mat_lin_diffusion.xml -c char_diffusion.xml -a algo_default_diffusion_noACV.xml -s resultats/r10

#AVEC ACV==============
# ARLEQUIN N3_3
mpirun -np 1 $AMITEX  -nz /usr/local/AMITEX/AMITEX_FFT/MICROSTRUCTURES/arlequin_N3_r3.vtk -m mat_lin_diffusion.xml -c char_diffusion.xml -a algo_default_diffusion.xml -s resultats/acvr3

# ARLEQUIN N3_9
mpirun $AMITEX  -nz /usr/local/AMITEX/AMITEX_FFT/MICROSTRUCTURES/arlequin_N3_r9.vtk -m mat_lin_diffusion.xml -c char_diffusion.xml -a algo_default_diffusion.xml -s resultats/acvr9

# ARLEQUIN N3_4
mpirun $AMITEX  -nz /usr/local/AMITEX/AMITEX_FFT/MICROSTRUCTURES/arlequin_N3_r4.vtk -m mat_lin_diffusion.xml -c char_diffusion.xml -a algo_default_diffusion.xml -s resultats/acvr4

# ARLEQUIN N3_5
mpirun $AMITEX  -nz /usr/local/AMITEX/AMITEX_FFT/MICROSTRUCTURES/arlequin_N3_r5.vtk -m mat_lin_diffusion.xml -c char_diffusion.xml -a algo_default_diffusion.xml -s resultats/acvr5

# ARLEQUIN N3_10
mpirun $AMITEX  -nz /usr/local/AMITEX/AMITEX_FFT/MICROSTRUCTURES/arlequin_N3_r10.vtk -m mat_lin_diffusion.xml -c char_diffusion.xml -a algo_default_diffusion.xml -s resultats/acvr10

# ARLEQUIN N3_10
mpirun $AMITEX  -nz /usr/local/AMITEX/AMITEX_FFT/MICROSTRUCTURES/arlequin_N3_r10bis.vtk -m mat_lin_diffusion.xml -c char_diffusion.xml -a algo_default_diffusion.xml -s resultats/acvr10bis


