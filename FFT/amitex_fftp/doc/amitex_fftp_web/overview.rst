.. _overview:


Introduction
============


The first purpose of **AMITEX_FFTP** is to implement an efficient **distributed** solver based on **FFTs** for **non-linear** mechanical simulations on **heterogeneous unit-cells** (discretized by 3D structured meshes). Following the idea of multi-physics problems, one variable diffusion (e.g. thermal diffusion) has ben added as a first step.

AMITEX_FFTP strongly relies on the **2DECOMP&FFT** library (http://www.2decomp.org) that provides:

* a user-friendly programming interface to work with a 2D pencil decomposition for data distribution on distributed-memory platforms,
* an interface with most popular external FFT libraries,
* parallel I/O.

For the user interface with XML files, AMITEX_FFTP relies on **FoX** (http://www1.gly.bris.ac.uk/~walker/FoX/) a Fortran library for XML.

----

Features
^^^^^^^^

Since AMITEX_FFTP is built on 2DECOMP&FFT, it can largely benefit from similar advantages: scalable, flexible, user-friendly and portable.

The user interface of AMITEX_FFTP allows to:

* "read" the heterogeneous microstructure from VTK files (legacy file format),
* "assign" and "distribute" material behaviors and their corresponding properties over the 2D pencil decomposition proposed by 2decomp.

Fortran objects (within Fortran modules) associated to the material, the loading, the algorithm parameters, are created at the beginning of the main program and all the implementation of the solver relies on these objects.

The mechanical (and one-species diffusion) FFT-based solver within AMITEX_FFTP allows to:

* Solve non-linear mechanical problems on unit-cells with prescribed average stress or strain components,
* Solve one species diffusion problems on unit-cells with prescribed average flux or gradient,
* Choose between a small perturbation or a finite strain framework (mechanics)
* Choose between two algorithms: the classical fixed-point algorithm (see Moulinec 1994) and an accelerated algorithm (to be published),
* Choose between classical or filtered discrete Green operators.

The mechanical behavior law is evaluated within a standard *umat* procedure which ensures compatibility with both:

* The **CAST3M** (http://www-cast3m.cea.fr) Finite Element software, simplifying the cross comparison between the two codes,
* The **MFront** (http://tfel.sourceforge.net/index.html) code generator, simplifying the process of behavior law implementation.


----

Distribution
^^^^^^^^^^^^

AMITEX_FFTP should be available soon with a free license for research and education purposes similar to the CAST3M license (see `here <http://www-cast3m.cea.fr/CT_INTERNET_Cast3M_recherche_FR_ANG_2011_version_23032011_SP.pdf>`_).

For now, it is accessible through a **GitLab** project https://gitlab.maisondelasimulation.fr/jderouil/amitex_fftp (restricted access).

----

Contributors
^^^^^^^^^^^^

.. todo : voir les accents en html

* Lionel Gélébart, Engineer - Researcher, CEA Saclay/DEN/DMN/SRMA
* Julien Derouillat, Engineer - Researcher, Maison de la Simulation
* Nicolas Doucet, Internship (12 months), CEA Saclay/DEN/DMN/SRMA
* Franck Ouaki, Post-doctoral (19 months), CEA Saclay/DEN/DMN/SRMA
* Fabien Bernachy-Barbé, Engineer - Researcher, CEA Saclay/DEN/DPC/SECR

----

Contact
^^^^^^^

* lionel.gelebart at cea.fr	
* julien.derouillat at cea.fr


