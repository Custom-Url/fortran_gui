#!/bin/bash

# Trap Ctrl+C and kill all child processes
trap 'echo "Interrupted. Killing child processes..."; kill 0; exit 1' INT

# Set environment
source /mnt/data/dg765/FFT
export LD_LIBRARY_PATH=/mnt/data/dg765/micro/FFT/amitex_fftp/libAmitex/lib:$LD_LIBRARY_PATH

# AMITEX executable and MPI
AMITEX="/mnt/data/dg765/micro/FFT/amitex_fftp/libAmitex/src/user_Miehe2/amitex_fftp"
MPIRUNalias='mpirun -np 18'

# Constant files
MATE="mate/mate_PF.xml"
LOAD="load/char.xml"
ALGO="algo/algo_Miehe2.xml"

# Loop through all generated VTK files in the nested directory structure
for vtk_file in micr/code_python/mesh/L*/h0.0002/iUC*.vtk; do
# for vtk_file in micr/code_python/mesh/L0.05_vf0.35/h0.0002/iUC*.vtk; do

    spacing_dir=$(basename "$(dirname "$vtk_file")")        # e.g. h0.0002
    vf_dir=$(basename "$(dirname "$(dirname "$vtk_file")")") # e.g. L0.05_vf0.35
    vtk_basename=$(basename "$vtk_file" .vtk)                 # e.g. iUC100_vpMIN0.09

    resu_dir="resu/dg765/${vf_dir}/${spacing_dir}/${vtk_basename}"
    mkdir -p "$resu_dir"


    echo "Running AMITEX for: $vtk_file"
    echo "Results in: $resu_dir"

    output_prefix="${resu_dir}/Load0.0"
    $MPIRUNalias $AMITEX -nm "$vtk_file" -m "$MATE" -c "$LOAD" -a "$ALGO" -s "$output_prefix"
done



