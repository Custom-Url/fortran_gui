.. _FAQ:

.. highlight:: xml

FAQ
===========================

Below are some Frequently Asked Questions or troubles faced by users. 

`Execution`_

`Installation`_

`Input vtk files`_

---------------------------

.. _Execution:

Execution
^^^^^^^^^

* **The code runs interactively on the host machine but not in batch mode on remote machines.**

This problem, observed on our cluster, was solved by enabling ssh connections by the command below: ::
 	
	enablessh

---------------------------

.. _installation:

Installation
^^^^^^^^^^^^

* **Advice.**

Try to be 'homogeneous' when using compiled external libraries: *openMPI* or *intelMPI*, *fftw* and *mfront* (optionnal):
	
	- use the same type of compiler (*intel* or *gcc*), and if possible the same version, as the compiler used to compile AMITEX_FFTP.

  
* **Troubles when installing AMITEX_FFTP at the installation of 2decomp&fft.**

Check with your linux administrator that *fftw* has been installed twice : first with double precision (default) and second, with single precision.
If installing *fftw* by your own, see http://www.fftw.org/doc/Installation-on-Unix.html. An example, here with the *intel* compiler, is given below: ::

        ./configure CC=icc --prefix=/home/toto/fftw-3.3.4_build 
        make
        make install
		
        ./configure CC=icc --prefix=/home/toto/fftw-3.3.4_build --enable-float 
        make
        make install

---------------------------

.. _input-vtk-files:

Input vtk files
^^^^^^^^^^^^^^^

* **My vtk files can be read by paraview but not by AMITEX_FFTP.**

In AMITEX_FFTP, the vtk-reader is 'home-made' and the header described in :ref:`geometry` must be strictly respected.

* **With very large images (over 2 billions of voxels), paraview doesn't read my '.vtk' file while AMITEX_FFTP accept it.**
 
Actually, paraview doesn't support such large files with the vtk legacy file format.
To be read in paraview, you must convert your '.vtk' files in '.vti' files.
It is planned to add a 'vti' reader in AMITEX_FFTP.


