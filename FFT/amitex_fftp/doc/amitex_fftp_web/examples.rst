.. _examples:

.. _example_polyx_mfront:

Example: polyx_mfront
======================

.. Examples
.. ==========

.. **********************
.. Example : polyx_mfront
.. **********************

The example described below demonstrates on a large polycrystalline simulation:

* The ability of AMITEX_FFTP to take advantage of massive parallelism (simulation performed on 1024 cores)
* The ability of AMITEX_FFTP to take into account a behavior law elaborated with the code generator **MFront** developed at **CEA** (http://tfel.sourceforge.net/documentations.html)

The algorithm converged with a total number of 102 iterations and the computation time, on 1024 cores (64 nodes, bi-processors Sandy Bridge E5-2670), was approximatively 4 hours.
  
The microstructure (A), the average behavior (B) and the local stress (C) and strain (D) fields are presented below. 	

+------------------------------+------------------------------+
| .. figure:: _static/Ex1.png  | .. figure:: _static/Ex2.svg  |
|    :width: 100%              |    :width: 100%              |
+------------------------------+------------------------------+
| *Microstructure (A)*         | *Average behavior (B)*       |
+------------------------------+------------------------------+

|

+------------------------------+------------------------------+
| .. figure:: _static/Ex3.png  | .. figure:: _static/Ex4.png  |
|    :width: 100%              |    :width: 100%              |
+------------------------------+------------------------------+
| *Local stress (C)*           | *Local Strain (D)*           |
+------------------------------+------------------------------+

----


Geometry
##########
The unit-cell consists of 42875 grains (Voronoï tessellation).

The unit-cell is discretized with a 1024x1024x1024 grid (almost 30x30x30 voxels per grain).

----

Behavior law
##############
The behavior law accounts for crystalline plasticity for Cubic Face Centered crystals in a small perturbation framework (description in the *Code_Aster* documentation R5.03.11). The implementation is provided with the code generator **MFront** (http://tfel.sourceforge.net/documentations.html) using the so-called "UMAT" interface, also compatible with the finite element code **CAST3M** (http://www-cast3m.cea.fr/). The material parameters are those furnished with the **MFront** code with 42875 crystalline orientations (one per grain) randomly distributed with an isotropic crystallographic texture.

----

Specificities
##############
The number of internal variables (required at each voxel) is 54.

The numerical integration of such behaviors is quite "heavy".

----

Loading 
#########
The applied loading is a tensile test loading (axial strain component and zero-stress on the five other components are imposed), at a strain rate of :math:`10^{-4}s^{-1}` until 1%.

The loading is uniformly discretized in 100 steps (1 second per step).

----

Algorithm 
###########

The algorithm is the classical fixed-point algorithm.
The discrete Green operator is based on a cubic filter (of size one voxel).

The criterion (stress divergence and applied load) is :math:`10^{-4}`.
