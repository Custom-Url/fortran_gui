AMITEX FFTP Program
===================

1 - Program Installation
------------------------
Set the desired compiler in the `FC` variable (e.g., `ifort` or `gfortran`).  
Specify the folder where the FFT library is located (nothing to do for Maldives).  
Ensure required modules are loaded beforehand (e.g., `gnu` or `intel`, `mpi`, and `tfel`).  
Start the installation with the command:  
`*./install*`

2 - Running the Program
-----------------------
In the *cas_tests* folder, the script *script.sh* allows you to run several test cases.  
To run it on the Maldives cluster, simply execute:

    qsub script_maldives

You can specify in this script the desired number of nodes or the chosen node.  
A script compatible with the Poincare cluster at the Maison de la Simulation is also available.  
In this case, run the following command from a login node:

    llsubmit script_poincare

3 - Plotting from Averaged Values
---------------------------------
The standard outputs (std, mstd, and zstd) contain average data  
(respectively over the whole cell, per material, and per zone).  
These are written in columns, and the quantity corresponding to each column is described  
in comments at the beginning of the file. Thus, results can be easily visualized using gnuplot.  
An example gnuplot script is provided in the *post/plot* directory.  
It can be executed after running the validation tests.

4 - Deformed Shape Visualization
--------------------------------
In the *post/deformed_shape* directory, the *deformedShape* program allows visualization of  
output fields in a VTK file where voxels are mapped to their deformed (current) configuration  
(for large deformations). This program is run using:

    ./deformedShape input_root [output_root]

The first parameter is the root of the input file(s),  
the second parameter is optional and defines the output file root (default: *sortie*).  
The program will search for *input_root_def.vtk* or, if not found,  
for *input_root_defi.vtk* (*i* from 1 to 9).  
This field (representing the displacement gradient) is used to build the displacement field  
and compute the coordinates of each voxel in the deformed configuration.  
The displacement gradient is then rewritten in this configuration  
to the file *output_root.vtk*.  
Additionally, if *input_root_sig* or *input_root_pi* files are present  
(representing Cauchy and Piola-Kirchhoff stress fields, respectively),  
those variables will also be included in *output_root.vtk*.

5 - For Developers: Validation and Result Updates
-------------------------------------------------
In the *validation* folder, the script *script_tests.sh* runs selected calculations  
and checks that the results remain consistent with previously obtained results.  
The script *update_results.sh* updates the reference results  
by taking those obtained with the new version of the program.

IMPORTANT
---------
Before each push, the developer must run the tests to ensure results are correct  
and that the new developments do not cause regressions in the code.  
These tests should be run with code compiled using `intel15`.  
You must then update the test results by running the *update_results* script.  
Before each commit, it is recommended to run the *clean* script  
to avoid pushing local file paths.

TRANSLATION PROGRESS
--------------------
Currently halfway through updating  NL_base_mod.f90

