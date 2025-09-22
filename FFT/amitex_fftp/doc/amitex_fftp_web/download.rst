.. _download:

Download
========

Here you can freely download AMITEX_FFTP for Education & Reserach purposes.

With the downloaded files you can:

	* **USE** the code, on a standalone machine as well as on HPC clusters,
	* **DEVELOP** your own algorithm or modifications as most of the fortran sources are provided,
	* If you want to DEVELOP, feel free to contact us for technical questions. Note also that you are welcome if you want to **SHARE** your developments in the distributed current version. 

Just keep in mind that AMITEX_FFTP is no a "big" code, there is not so many people behind it ;-)

---------------------------

Version 4.0.0
^^^^^^^^^^^^^

`Download v4.0.0 <_static/download_amitex_fftp.html>`_

Important modifications have been introduced with respect to version 2.3.5.
The most important for the users are :

	* a simplified way to introduce user defined behaviours, still compatible with **MFront** (http://tfel.sourceforge.net/index.html)
	* a new possiblity to take into account diffusion problems (stationary problems with a unique variable for the moment, like thermal diffusion)
	* the dimensions of the voxels (SPACING in *vtk* files) are now taken into account
	* bug corrections (especially the problem with the recent versions of openMPI)

.. warning::

	**Version 2.3.5 input xml files are no more compatible**

	**These files have to be slightly adapted**


Validation tests have been run on various clusters, or standalone PC, but only with intel XEON processors, with the following results:

 
---------------------------

Version 2.3.5
^^^^^^^^^^^^^

Only minor modifications have been introduced with respect to version 2.3.4.
The most important for the users are :

	* bug fixed : writing of messages on the standard output when using gcc compiler
	* introduction of an anisotropic elastic behaviour (Law_Number : 3)

Validation tests have been run on the CEA cluster (maldives).
Compared to version 2.3.4, validation tests on the other platforms should no be affected.
 
---------------------------

Version 2.3.4
^^^^^^^^^^^^^

Validation tests have been run on various clusters, or standalone PC, but only with intel XEON processors, with the following results:

   +------------------------+------------+------------------+--------+
   | CCRT cluster (airain)  | intel 14   | bullxmpi 1.2.8.4 | OK     |
   +------------------------+------------+------------------+--------+
   | MDS cluster (poincare) | gcc 4.9.0  | openMPI 1.8.4    | OK     | 
   +------------------------+------------+------------------+--------+
   |                        | gcc 4.7.2  | openMPI 1.6.3    | OK     |
   +------------------------+------------+------------------+--------+
   |                        | intel 13   | intelMPI 4       | OK     |
   +------------------------+------------+------------------+--------+
   |                        | intel 15   | intelMPI 5       | OK     |
   +------------------------+------------+------------------+--------+
   | CEA cluster (maldives) | intel 15   | intelMPI 5       | OK     |
   +------------------------+------------+------------------+--------+
   |                        | gcc 4.8.2  | openMPI 1.6.5    | OK     |
   +------------------------+------------+------------------+--------+
   |                        | gcc 4.8.2  | openMPI 1.8.1    | **PB** |
   +------------------------+------------+------------------+--------+
   | Standalone PC          | gcc 4.9.2  | openMPI 1.6.5    | OK     |
   +------------------------+------------+------------------+--------+


The problem with openMPI 1.8.1 was observed on the validation test *beton_relax_65_old*. 


.. warning::
   Performances obtained with [Intel compiler + IntelMPI] are significantly better than with [GCC compiler + openMPI].

