.. _install:

.. highlight:: bash

Installation 
--------------

The code should be installed on **Linux systems**, on which it is developped and tested.

**AMITEX_FFTP** has been installed successfully on **macOS** a few users: see :ref:`mac_install` if you want to try it.

Installation on **WINDOWS** has been tested once with **CYGWIN** (https://www.cygwin.com/)... try it if you have time ;-)

------------------------------------------------------------------------------------------------------------

Prerequisites
^^^^^^^^^^^^^

If possible, install GCC, opemMPI and FFTW with package manager (such as yum, apt... depending on you linux distribution).

If really not possible, you can install them by your own in your home directory (see :ref:`gcc_openmpi_fftw_mfront`). 

**1. Fortran compiler** : Intel 15.0.0 (minimum tested) or GCC > 5.0.0 

The program (v8.17.1) has been tested successfully with Intel 19 and GCC 9.2.0. 

**2. MPI library** : The program has been tested successfully with at least : Intel MPI 5.0.1 (2015)...2019 as well as Open MPI 2.0.4 and 3.1.0

.. Warning::
   Observations made on intel-based clusters
	**Performances obtained with [Intel compiler + Intel MPI] are significantly better than with [GCC compiler + Open MPI]**.

        **While waiting for a clarification, use [Intel compiler + Intel MPI] if possible**.



**3. FFT libraries** :

As mentioned, **AMITEX_FFTP** comes with the **2decomp&fft** library which provides interfaces with various FFT libraries.

AMITEX_FFTP has been tested successfully with libraries **fftw** (see [FFTW32005]_), versions 3.3.3 and 3.3.4. If **2decomp&fft** offers the possibility to choose different FFT engines (acml, ffte, fftpack5, mkl), only **fftw** is used for the development of **AMITEX_FFTP**.

Note that **fftw** (http://www.fftw.org/) must be installed with both *double* (default) and *simple precision*.


**4. The MFront library (optional)**

The **MFront** library (http://tfel.sourceforge.net/documentations.html) is useful to easily develop behavior laws. These laws are available in a dynamic library and are based on the standard *UMAT* interface which is compatible with **AMITEX_FFTP**. 

The **MFront** library should be installed and compiled with the same compiler than the one used to compile **AMITEX_FFTP**.

------------------------------------------------------------------------------------------------------------

.. _installation_guide:

Installation guide
^^^^^^^^^^^^^^^^^^
.. Admonition:: Install! 

	The installation is straightforward : **adust** the *install* file and **launch it**.

	Then the executable program **amitex_fftp** is located in *(amitex_fftp)/libAmitex/bin* folder

1. Adjust the following first lines of the *install* file located in the main folder *(amitex_fftp)*::


	# Compiler : ifort, gfortran
	FC=ifort
	# FC=gfortran

	# MPI compiler
	MPIFC=mpif90
	#MPIFC=mpiifort


	# FFTW Library
	#   fftw3      - FFTW version 3.x
	#   fftw3_f03  - FFTW 3.3-beta1 or later (with Fortran 2003 interface)
	#   generic    - A general FFT algorithm (no 3rd-party library needed)

	FFT=fftw3_f03

	# FFTW library path
	FFT_inc=/usr/include
	FFT_lib=/usr/lib64


The user must choose :

	* the compiler: ifort (Intel) or gfortran (GCC) (prefer ifort if possible),

	* the wrapper compiler for the compilation with MPI : mpif90 or mpiifort 

	* the FFT library and inform the corresponding library path. Prefer fftw_f03 used during the development.


2. Launch install::

	./install


In practice, ./install :

	* fills correctly *(amitex_fftp)/lib_extern/2decomp_fft-1.5.8472/src/Makefile.inc* and various Makefiles and scripts,

	* compiles the **2decomp&fft** library,

	* compiles the **FoX** library, 

	* compiles **AMITEX_FFTP**,

	* builds the executable program: the file *amitex_fftp* in the *(amitex_fftp)/libAmitex/bin* folder.

.. Warning::

	If the optional **MFront** library is not installed, the compilation of  *(amitex_fftp)/cas_tests/comportements/mazars* and *(amitex_fftp)/cas_tests/comportements/mgurtin_hpp_glissement* will fail, obviously. 

	If you do not wish to use **MFront** do not pay attention to this error. 


------------------------------------------------------------------------------------------------------------

Validation
^^^^^^^^^^^^^

The file *script_tests.sh* within the folder *(amitex)/validation* provides a list of simulations based on geometries, loadings, algorithms and material properties defined in xml files respectively in directories *microstructures*, *chargements*, *algorithmes* and *materiaux* under *(amitex_fftp)/cas_tests*. 

The script *script_tests.sh* can be run either interactively ``./script_tests.sh`` or in batch mode. In this cas, different batch scripts for different batch schedulers (Torque, LoadLeveler...) are provided. Of course, you will have to create your own batch script in agreement with the batch scheduler available on your platform.

Finally, *script_test.sh* compares simulations results to reference results stored within each *(amitex)/resultats/simulation_name/reference* folder and give the non-regression comparison in the file *tests.log* and total time in *time.log*.


------------------------------------------------------------------------------------------------------------

.. _gcc_openmpi_fftw_mfront:

GCC, openMPI, fftw, mfront 
^^^^^^^^^^^^^^^^^^^^^^^^^^^

If really impossible to install GCC, openMPI, fftw with a package manager, try the following installation procedures in you home directory.

**For GCC** ::

	wget ftp://ftp.irisa.fr/pub/mirrors/gcc.gnu.org/gcc/releases/gcc-8.3.0/gcc-8.3.0.tar.gz
	tar xvfz gcc-8.3.0.tar.gz
	cd gcc-8.3.0
	./contrib/download_prerequisites
	mkdir build
	cd build
	../configure --enable-languages=c,c++,fortran --disable-multilib --prefix=<PATH TO NEW VERSION>
	make -j 4     # distributed compilation on 4 process
	make install

Then you must adjust the PATH and LD_LIBRARY_PATH, for example in a file env_gcc8.3.0.sh ::
	
	#!/bin/bash
	export PATH=<PATH TO NEW GCC VERSION>/bin:$PATH
	export LD_LIBRARY_PATH=<PATH TO NEW GCC VERSION>/lib64:$LD_LIBRARY_PATH

To use te compiler ::
	
	source env_gcc8.3.0.sh
	gfortran -v  #to check the version

**For openMPI**

See https://www.open-mpi.org/faq/. 

Below an example of standard (no specific options) installation of openMPI ::

	gunzip -c openmpi-4.0.2.tar.gz | tar xf -
	cd openmpi-4.0.2
	./configure --disable-multilib --prefix=<PATH TO NEW VERSION> 
	make all install

Then you must adjust the PATH and LD_LIBRARY_PATH, for example in a file env_openmpi.sh ::

	#!/bin/bash
	export PATH=<PATH TO NEW OPENMPI VERSION>/bin:$PATH
	export LD_LIBRARY_PATH=<PATH TO NEW OPENMPI VERSION>/lib:$LD_LIBRARY_PATH

To use the compiler ::
	
	source env_openmpi.sh
	mpif90 -v  #to check the version


**For FFTW**

See http://www.fftw.org/fftw3_doc/Installation-on-Unix.html.

Below an example of standard (no specific options) installation of FFTW with intel compiler (single and double precisions required for **AMITEX_FFTP**) in the main folder ::

	#install double-precision
	./configure CC=icc --prefix=<PATH TO NEW FFTW VERSION> 
	#remove CC=icc for install with gcc instead of intel
	make
	make install

	#install single-precision
	./configure CC=icc --prefix=<PATH TO NEW FFTW VERSION> --enable-float
	#remove CC=icc for install with gcc instead of intel
	make
	make install


Then you must adjust the PATH and LD_LIBRARY_PATH, for example in a file env_fftw.sh ::

	#!/bin/bash
	export PATH=<PATH TO NEW FFTW VERSION>/bin:$PATH
	export LD_LIBRARY_PATH=<PATH TO NEW FFTW VERSION>/lib:$LD_LIBRARY_PATH


**For MFRONT**

Download the code at https://sourceforge.net/projects/tfel/files/

For a complete and up-to-date description refer to the manual page http://tfel.sourceforge.net/install.html

Below is a summary that should work with you native gcc compiler or with intel compiler (choose the correct lines) ::

	mkdir tfel-3.2.1-gcc	# for installation with gcc (native gcc compiler)
	mkdir tfel-3.2.1-intel	# for installation with intel

	tar -xf tfel-3.2.1.tar.bz2
	cd tfel-3.2.1

	mkdir build
	cd build

	#for native gcc
	cmake ../ -DCMAKE_BUILD_TYPE=Release -Dlocal-castem-header=ON -Denable-fortran=ON -Denable-aster=ON -DCMAKE_INSTALL_PREFIX=path_to_tfel3.2.1-gcc & >cmake.log

	#for intel
	CXX=icpc CC=icc FC=ifort F77=ifort cmake ../ -DCMAKE_BUILD_TYPE=Release -Dlocal-castem-header=ON -Denable-fortran=ON -Denable-aster=ON -DCMAKE_ INSTALL_PREFIX=path_to_tfel-3.2.1-intel

	make -j 8		# 8 is the number of core for a parallel compilation
	make check -j 8		# launch tests
	make install		# finalize installation in tfel-3.2.1-gcc or tfel-3.2.1-intel

	# For a non-native gcc compiler : add CXX=g++ CC=gcc FC=gfortran F77=gfortran 
        #                                 to the cmake command line above


Then you must adjust the PATH and LD_LIBRARY_PATH, for example in a file env_mfront.sh ::

	#!/bin/bash
	export PATH=<PATH TO NEW TFEL VERSION>/bin:$PATH
	export LD_LIBRARY_PATH=<PATH TO TFEL VERSION>/lib:$LD_LIBRARY_PATH

------------------------------------------------------------------------------------------------------------

.. _mac_install:

macOS specificity
^^^^^^^^^^^^^^^^^^^

macOS installation is not tested by developers but **AMITEX_FFTP** has been successfully installed by different users on different configurations.

An example of successfull configuration is, on a Catalina-OS (01/2020):

	* Xcode + command line tools, 
	* gcc9, 
	* open-mpi for gcc9,
	* FFTW 3.3.8.

A problem comes with macOS because the commande 'sed' used in the scripts *install* and *clean*, has not the same behavior on macOS and Linux.

To overcome this issue :
	
	* install 'gsed' on macOS (sudo port install gsed)

	* Adjust the variable SED in scripts *install* **AND** *clean* so that 'gsed' is used instead of 'sed' :: 

		#------- For MAC USERS : uncomment #SED=gsed
		SED=sed
		SED=gsed

Once modified *install* **AND** *clean*, the installation procedure is the same as described in :ref:`installation_guide`. 


