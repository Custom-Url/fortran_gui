.. _user:

.. highlight:: xml


User interface 
===============

.. _prior-definitions:

Prior definitions
-----------------

A heterogeneous *Unit-cell* consists of one (or more) *Material*, and each *Material* consists of one (or more) *Zone*.
They are defined based on the following properties:

* All the voxels of the 3D uniform grid define the *Unit-cell*.
* Within a *Material*, all the voxels have the same behaviour law.
* Within a *Zone* of a *Material*, all the voxels have the same material properties (coefficients of the behavior law).

.. Note::

	The description of a unit-cell **is not unique**.

	For example a concrete can be described as :
	
	1. **2** *materials* (matrix and "inclusions"), all the "inclusions" belonging to a **unique** *zone*,
	
	2. **2** *materials* (matrix and "inclusions"), dividing "inclusions" into **N** *zones* (one *zone* per inclusion)

	3. **N+1** *materials* (matrix and N inclusions, with one *material* per inclusion and a **unique** *zone* per *material*)

	If multiple descriptions can be used, the choice among them is mainly related to the outputs proposed by *amitex_fftp* (per material and per zone statistical quantities, see below)

	Finally, it is recommended to reduce the number of materials (avoid proposition 3. if N is large)

--------------------------------------------------------------------------------------------------------------------

.. _how-to-launch:

How to launch AMITEX_FFTP
-------------------------

To use **amitex_fftp**, you must add (amitex)/libAmitex/bin to the PATH, for example in a file env_amitex.sh ::

	#!/bin/bash
	export PATH=(amitex)/libAmitex/bin:$PATH

And then ::
	
	source env_amitex.sh

If you have compiled GCC and/or openMPI and/or fftw3, and used them for compiling **amitex_ftp** (see :ref:`gcc_openmpi_fftw_mfront`), be sure that their respectives *bin* and *lib* (or *lib64*) folders are added to the PATH and LD_LIBRARY_PATH variables in the env_amitex.sh 


Interactively or within a batch script (for the job manager available on your machine), a simulation is launched with one of the following command lines ::
	
	MATEVTK="path_to/materialID.vtk"
	ZONEVTK="path_to/zoneID.vtk"
	MATEXML="path_to/material.xml"
	LOADXML="path_to/loading.xml"
	ALGOXML="path_to/algorithm.xml"
	RESDIR="path_to/results/

	#WARNING : you must create the folder "results" before launching amitex_fftp

	0. help 
	mpirun -np n amitex_fftp -help  ! OR SIMPLY :  amitex_fftp -help

	1. General case
	mpirun amitex_fftp  -nm $MATEVTK -nz $ZONEVTK -m $MATEXML -c $LOADXML -a $ALGOXML -s $RESDIR/output_file
 
	2. One material (assumes materialID.vt full of 1)
	mpirun amitex_fftp  -nz $ZONEVTK -m $MATEXML -c $LOADXML -a $ALGOXML -s $RESDIR/output_file

	3. One zone per material (assumes zoneID.vtk full of 1)
	mpirun amitex_fftp  -nm $MATEVTK -m $MATEXML -c $LOADXML -a $ALGOXML -s $RESDIR/output_file

	4. One material with one zone per voxel (assumes zoneID.vtk varying from 1 to the number of voxels) with dx=dy=dz=1.
	mpirun amitex_fftp    -m $MATEXML -c $LOADXML -a $ALGOXML -NX nx -NY ny -NZ nz -s $RESDIR/output_file

	5. One material with one zone per voxel (assumes zoneID.vtk varying from 1 to the number of voxels)
	mpirun amitex_fftp  -m $MATEXML -c $LOADXML -a $ALGOXML -NX nx -NY ny -NZ nz -DX dx -DY dy -DZ dz -s $RESDIR/output_file


Depending on the environment variables, you may have to specify the number of processes and replace ``mpirun`` by ``mpirun -np number_of_processes``. 

The five possible command lines listed above correspond to the three following cases :
 	
1. The *materials* and their corresponding *zones* are defined by *materialID.vtk* and *zoneID.vtk*.

2. The unit-cell consists of a **unique** *material* with **different** *zones* defined in *zoneID.vtk*.

3. The **different** *materials* defined by *materialID.vtk* have a **unique** *zone* (all the voxel within a material have the same properties).

4. The unit-cell consists of a **unique** *material* with one *zone* per voxel (the voxel dimensions are set to 1x1x1).

5. Same as 4. but the user provides the voxel dimensions



The material properties, the loading and the algorithm parameters are defined in the xml files : *material.xml*, *loading.xml* and  *algorithm.xml*.

The rootname for the output file is given by *output_file*.

--------------------------------------------------------------------------------------------------------------------

.. _geometry:

Geometry
------------

The input file(s) *materialID.vtk* and/or *zoneID.vtk* are integer fields defined on a uniform grid according to the *vtk simple legacy format* (http://www.vtk.org/VTK/img/file-formats.pdf).

All the voxels with a common materialID value define a *material* domain.
Within a given *material*, all the voxels with a common zoneID value define a *zone*.


**CONSTRAINTS FOR THE GENERATION OF VTK FILES**

**Data definition**

The materialID field must contain integers starting from 1 to N (N=number of *materials*), increasing without any skipping value.

**Within each material** domain,the zoneID field must contain integers starting from 1 to N (N=number of *zones*), increasing without any skipping value.

In both cases, definitions starting from 0 to N-1 are accepted but not encouraged.

**VTK file headers**

The VTK file header must **strictly contain the following 10 lines**::

	  	# vtk DataFile Version 4.5
		Materiau
	  	BINARY
	  	DATASET STRUCTURED_POINTS
	  	DIMENSIONS    66   66   66
	  	ORIGIN    0.000   0.000   0.000
	  	SPACING    1.000000    1.000000   1.000000
	  	CELL_DATA   274625
	  	SCALARS MaterialId unsigned_short
	 	LOOKUP_TABLE default

DIMENSIONS defines the number of voxels per side plus 1. Here, the grid consists of 65x65x65 voxels.

SPACING defines the voxel size and especially the voxel shape. Here voxels are cubic.

CELL_DATA is the number of voxels, it must be in agreement with DIMENSIONS, here 274625=65x65x65.

SCALARS defines the data type (*char*, *short*, *int*, *long*, *unsigned_char*, *unsigned_short*, *unsigned_int*, *unsigned_long*, respectively coded on 1, 2, 4 and 8 bytes), here the type is *unsigned_short* (2 bytes).

.. _Binary data:

**Binary data**

Binary data, given just after the header, **must satisfy the following constraints**:

* byte ordering used to write the binary data is *big endian*
* the data type when writing binaries must be consistent with the data type given in the header
* the maximum field value (especially for zoneID) must be less than the limit value corresponding to the data type

.. caution::
	
	In case of **unsigned** data types, the limit value used by AMITEX_FFTP is the limit value for corresponding **signed** data types (**~ 0.5 x limit value for unsigned**).
	
	For example, the maximum value for both *char* and *unsigned_char* data types (coded on 1 byte) is 127.


.. admonition:: TIPS

	1. You can use specific :ref:`tools` to generate vtk files for AMITEX_FFTP

	2. **if PARAVIEW can't read your vtk file don't go further with AMITEX_FFTP.**



--------------------------------------------------------------------------------------------------------------------

.. _XML-files:

XML files
--------------------

Of course, XML files must follow the XML syntax (use a XML editor with syntax coloring to detect syntax errors).

A comment is written like this::

	<!-- A comment -->

There is no need to follow any specific order for the various inputs given in XML files.

**Case sensitivity rules**

XML is a case sensitive language. The name of XML **nodes** and **attributes** are case-sensitive.

However, the name of the string variables that can be input between quotation marks are case insensitive. 

Hence, in the following XML line ::

      <Coeff Index="1" Type="Constant_Zone">	 

Here, *Coeff* is a **node**, *Index* and *Type* are **attributes** : they require the uppercase "C", "I" and "T". 

However, the string *Constant_Zone* can be read *constant_zone* in AMITEX_FFTP.

**Exception** : path names are case-sensitive.  

-----------------------------------------------------------------------------------------------

.. _Material-properties:

Material properties
--------------------

.. WARNING::

	**The number of materials and the number of zones (per material), provided in the "material.xml" file, must be consistent with the description given by the vtk file(s).**


**EXAMPLE**

This *material.xml* file below provides in a single example different ways to associate the geometry to material properties ::

	<?xml version="1.0" encoding="UTF-8"?>
	<Materials>

	    <!-- REFERENCE MATERIAL -->
	    <Reference_Material Lambda0="2.0952e+10" Mu0="1.5014e+10"/>


	    <!-- MATERIAL "2", elastic isotropic behavior                                      -->
	    <!--           prefer the syntax Lib="" for native behavior Law (see Material "1") -->
	    <Material numM="2" Lib="(amitex)/libAmitex/src/materiaux/libUmatAmitex.so" Law="elasiso">
	        <!-- Two constant coefficients (Lamé coefficients Lambda and Mu) -->
	        <Coeff Index="1" Type="Constant" Value="4.0385e+10"/>
	        <Coeff Index="2" Type="Constant" Value="2.6923e+10"/>
	    </Material>


	    <!-- MATERIAL "1", linear viscoelastic behavior                                     -->
            <!--               Lib="" possible for 'native' behavior Law                        -->
            <!--               possible if AMITEX_PATH is defined (export AMITEX_PATH=(amitex)) -->
	    <Material numM="1" Lib="" Law="viscoelas_maxwell">  

	    	<!-- coefficient 1: constant per zone -->
	        <Coeff Index="1" Type="Constant_Zone">
	            <!-- value within zone 1 -->
	            <Zone numZ="1" Value="5.235e+10" />
	            <!-- value within zone 2 -->
	            <Zone numZ="2" Value="3.912e+10" />
	        </Coeff>

	    	<!-- coefficient 2: constant per zone, values defined in the ASCII file "Coeff2" -->
	        <Coeff Index="2" Type="Constant_Zone " File="Coeff2" Format="ASCII" />

	    	<!-- coefficient 3: constant per zone, values defined in the binary file "Coeff3.bin" -->
	        <Coeff Index="3" Type="Constant_Zone " File="Coeff3.bin" Format="binary" />


	        <!-- Internal variables are defined similarly-->
	        <IntVar Index="1" Type="Constant" Value="0."/>

	        <IntVar Index="2" Type="Constant_Zone">

	            <Zone numZ="1" Value="0." />
	            <Zone numZ="2" Value="10." />

	        </IntVar>

	        <IntVar Index="3" Type="Constant_Zone" File="Var3"/>

	        <IntVar Index="4" Type="Variable" File="(rep_to_vtk)/coeff.vtk" Format="vtk"/> 

	    </Material>

	</Materials>


**BEGIN AND END**

The file must begin and end with ::
	
	<?xml version="1.0" encoding="UTF-8"?>
	<Materials>
		.
		.
	</Materials>


**REFERENCE MATERIAL** ::

	<Reference_Material Lambda0="2.0952e+10" Mu0="1.5014e+10"/>

The reference material is defined by its Lamé coefficients. For an optimal choice with the basic scheme, see [Moulinec1998]_\ . If convergence acceleration is used (see :ref:`algorithm-parameters`), the same choice can be used but the influence of this choice on the convergence is much less significant.

**MATERIALS**

Each *material* is associated to the geometry by the material number *numM* which directly refers to a material domain in *materialID.vtk*. 

Materials numbers *numM* starts to one (if *ID* in *materialID.vtk* starts to zero : *ID* = 0 corresponds to *numM* = 1, *ID* = 1 to *numM* = 2 and so on) ::

	<Material numM="2" Lib="(amitex)/libAmitex/src/materiaux/libUmatAmitex.so" Law="elasiso">
	   .
	   .
	</Material>
	<Material numM="1" Lib="" Law="elasiso">
	   .
	   .
	</Material>


One ``<Material.. > .. </Material>`` section must be defined per *material* present in the **materialID.vtk** (*i.e.* a single section if the file is omitted). 

The material behaviour laws are implemented through *UMAT* compatible procedures gathered in a dynamic library.

The fields ``Lib`` and ``Law`` provide the complete path to the library and the name of the procedure. 

In the standard version of **amitex_fftp**, only a few native behavior laws are provided within the library *(amitex)/libAmitex/src/materiaux/libUmatAmitex.so*.

For native behaviors, if the environment variable AMITEX_PATH is defined, the complete path can be avoided, using ``Lib=""``.

The available native behavior laws are :: 

	elasiso 		: isotropic elasticity 
	elasiso_GD 		: Lagrangian isotropic elasticity 
	elasaniso 		: orthotropic elasticity 
	thermoelasiso 		: isotropic thermoelasticity
	paramextelasiso         : isotropic 'thermo'elasticity but with an external parameter governing the dilatation (instead of the temperature)	
	contrainte_imposee 	: imposed stress
	viscoelas_maxwell	: Maxwell linear visco-elasticity


However, introducing a new behaviour law is quite simple and discribed below (see :ref:`user-behavior-label`).   

.. Caution ::
	The *numM* values, **must be defined between 1 and N** (N=number of *materials*), while the integer values distributed in the file *materialID.vtk* can be defined between 1 and N or 0 and N-1 (see `Geometry`_). 
	**If generating your own 'vtk' files, prefer a distribution between 1 and N.** 


**MATERIAL COEFFICIENT**

The material coefficient's index is the position of each coefficient within the vector COEFF used as an input of the *UMAT* compatible procedure for the behaviour law (see :ref:`user-behavior-label`). 

For an isotropic behavior (``Law="elasiso"``), the coefficients COEFF(1) and COEFF(2) are the Lamé coefficients :math:`\lambda` and :math:`\mu`, respectively.


In the example above, four different ways are used to define these values ::

	1. <Coeff Index="1" Type="Constant" Value="4.0385e+10"/>

	2. <Coeff Index="1" Type="Constant_Zone">
	          <Zone numZ="1" Value="5.235e+10" />
	          <Zone numZ="2" Value="3.912e+10" />
	   </Coeff>

	3. <Coeff Index="2" Type="Constant_Zone " File="Coeff2" Format="ASCII"/>

	4. <Coeff Index="3" Type="Constant_Zone " File="Coeff3.bin" Format="binary"/>


1. is used if the coefficient is constant within the *Material*
2. 3. and 4. are used if the coefficient is constant within each *Zone* defining the *Material*
3. is used to import coefficients from a text file, which can be useful when considering a relatively large number of zones (see :ref:`examples`). To define a coefficient for two zones, the file must contain at least two lines (one value per line) see an example of the file Coeff2: ::

		1.429430417
		3.905850691

   If the *Format* tag is not specified, the program will assume that the file given is written in the ASCII format.

.. _binary-files:

4. is used to import coefficients from a binary file. This is very useful when there is a very large number of zones. The header must contain the following **2 lines**::

               2
               double

   The first line is the number of binary values within the file, it must greater or equal to the number of zones defined in the VTK file.

   The second line defines the data type (*char*, *short*, *int*, *long*, *unsigned_short*, *unsigned_int*, *unsigned_long*, respectively coded on 1, 2, 4 and 8 bytes, or *float* and *double* on 4 and 8 bytes), here the type is *double*. The rest of the file should be binary data satisfying the constraints defined in `Binary data`_ (*big endian* ordering and maximum values).

**INITIAL INTERNAL VARIABLES** 

Initial internal variables are assigned according to the same procedure. Here, the index is the position of each internal variable within the vector STATEV (used as an input/output of the *UMAT* compatible procedure for the behaviour law, see :ref:`user-behavior-label`).

A fifth possibility is proposed to assign internal variables with binary *vtk* files (format *float* or *double*) ::

	5. <IntVar Index="4" Type="Variable" File="(rep_to_vtk)/variable.vtk" Format="vtk"/> 

The *vtk* files definition is the same as the definition used for the geometry (see :ref:`geometry`), but here with format *float* or *double*.

This is especially interesting to assign a continuously variable field (a random field for example).

If the *Material* do not cover the whole unit-cell (which occurs with more than one *Material* in the unit-cell), only the corresponding part of the complete field given in 'variable.vtk' is used. The rest of the field is not considered.


-----------------------------------------------------------------------------------------------

.. _Composite_voxels:

Composite voxels (advanced users)
------------------------------------------------------------
.. Important::

	- The important idea is that **the same vtk files can be used with or without composite voxels**. 
	
	- The introduction of composite voxels in the simulation is done by **additionnal informations within the material.xml file**.  

	- The definition of composite voxels requires a **non-straightforward pre-treatment** to precisely identify where they are located and what they consist of (which materials, which volume fractions...)

Introduced progressively since version 5.0.0, the final goal is to be able to account for:

	1. different averaging rules (at least Voigt, Reuss and Laminate)
	2. any number of phases
	3. any linear or non-linear behavior
	4. Small Strains and Finite Srains frameworks 

In v8.17.1, Finite Strain extension is not yet available.
The other points are fullfilled.

**GENERAL DESCRIPTION**

The introduction of composite voxels in the file *material.xml* must follow the three different points : 

	**1. Sections** ``<Material.. > .. </Material>``
	 
	All the materials present in the unit-cell are defined as previously (see :ref:`Material-properties`), 
	whether they are present or not as homogeneous voxels. 

	A node ``<Coeff_composite .. />`` must be added to provide an elastic isotropic behavior
	for the numerical integration of the average behavior (Laminate or Reuss).
	The two parameters are respectively the Lamé coefficients :math:`(\lambda,\mu)`
	
	**2. Section** ``<Material_composite> .. </Material_composite>``

	The path of the directory were the complete description of composite voxels is given 
	according to a precise format given below.
	
	**3. Section** ``<Interphase> .. </Interphase>`` (if required)

	List of the materials and/or zones which do not appear as homogeneous voxels.    

	Actually, sometimes, a specific ID of *materialID.vtk* or a specific zone of *zoneID.vtk* can be completely 
	overlaid by the composite voxels so that it doesn't appear anymore as homogeneous voxels.
	This is especially the case when dealing with thin interphases.
	
.. `Format-for-composite-voxels-definition`_:
	
**FORMAT FOR COMPOSITE VOXELS DEFINITION**

The whole set of composite voxels is divided into families, each being defined by:
	- a common set of material IDs (``numM`` in :ref:`Material-properties`)  
	- an averaging rule (`voigt`, `reuss` or `laminate`)

The complete definition of composite voxels is given in a single directory containing:

	**1. One file** ``list_composite_materials.txt`` defining the different families of composite voxels, one line per family.
	Each line has the form ``I J K ... P averaging_rule`` where I, J, K ... P are the material IDs present in this family.
	In the example below, two families are defined::
	
		1 1 voigt
		1 2 4 6 laminate
	
	The first one gathers all the composite voxels with :
		* two phases associated to the same material ID 1 (``numM`` in :ref:`Material-properties`). Here the mechanical contrast comes from the coefficients which can be different from one zone ID to another
		* the voigt avering rule,
	and the second one, with : 
		* 4 phases, associated to the material IDs 1, 2, 4 and 6 (``numM`` in :ref:`Material-properties`)
		* the laminate averaging rule.
		
	**2. One directory** ``rep_I_J_K..._P`` **per family** (``rep_1_1`` and ``rep_1_2_4_6`` in the example).
	Within each directory, a local phase renumbering is considered : I :math:`\leftrightarrow` 1, J :math:`\leftrightarrow` 2, K :math:`\leftrightarrow` 3...
	and the following files are contained 
	
		* *pos.bin* (for *i* =1 .. number of phases) : the linearized positions of the composite voxels  
		* *zonei.bin* (for *i* =1 .. number of phases) : their zone IDs for phase *i* (local phase number)  
		* *fvi.bin*   (for *i* =1 .. number of phases) : their volume fractions of phase *i* (local phase number)
		* *Nijx.bin*, *Nijy.bin*, *Nijz.bin* (for *i*, *j* =1 .. number of phases and *i* :math:`\ne` *j*) : their components of the normal vector to the interface between phases *i* and *j*.
		* *Sij.bin* (for *i*, *j* =  .. number of phases and *i* :math:`\ne` *j*) : their surfaces of intersection between the voxel and the interface between phases *i* and *j*.
		
Each file contains a list of binary values, one for each voxel composite of the family, respecting the format binary-files_.

.. Note ::

	In practice, the mandatory files are:
		
		- For all averaging rules : *pos.bin*, *zonei.bin* (for *i* =1 ..  number of phase), *fvi.bin* (for *i* =1 ..  number of phase -1)
	
		- For laminate : *N12x.bin*, *N12y.bin*, *N12z.bin* (an arbitrary choice is made to keep the normal *N12* for the laminate model)
		
	However, all the normals *Nij* and surfaces *Sij* could be used in the futur, and can be usefull for evaluating interfaces average quantities in post-treatment.

		

**INTERPHASE DEFINITION**

The section ``<Interphase> .. <\Interphase>`` **must be added** to the file *material.xml* **IF AND ONLY IF, after adding the composite voxel's definition**:

	1. A material ID, *matID*, (``numM`` in :ref:`Material-properties`) do not appear anymore as a homogeneous voxel.
	Such an "Interphase" material **must be** identified through::
	
		<Interphase_material  numM="matID" Nzones="number_of_zones"/>
	
	As all the materials in **amitex_fftp**, its coefficients can be associated to different zones 
	so that the number of zones ``Nzones`` **must be specified** (unlike the other materials, 
	the number of zones can't be deduced from the files *materialID.vtk* and *zoneID.vtk* associated to homogeneous voxels).
	
	2. At least one zone ID of a given material ID, *matID*, do not appear anymore as a homogeneous voxel.
	Such "Interphase" zones **must be** identified through::
	
		<Interphase_zone_list numM="matID" Nzones="number_of_zones">
			zoneID1 zoneID2 ... zoneIDnumber_of_zones
		</Interphase_zone_list>  
	
	

Below is an example where material 2, which consists of 4 zones (associated to different material coefficients) 
does not appear anywhere as a homogeneous voxel. In addition, material 1 appears as a homogeneous voxel but not its three zones 2, 4 and 6.::

	<Interphase>
		<Interphase_material  numM="2" Nzones="4">  
		<Interphase_zone_list numM="1" Nzones="3">
			2 4 6
		</Interphase_zone_list>  
	</Interphase>


**EXAMPLE 1 : a voronoï polycrystal with isotropic grains**

.. Note::

	In this exemple, and in any case where the definition of **"interphase" is not required**
	(all the material and material's zones are present as homogeneous voxels), the simulation **can be run** 
	without composite voxels by simply commenting the section 	``<Material_composite   > .. </Material_composite>``.


In this example, there is only one material and the file *zoneID.vtk* consists of a voronoï tessellation (one grain per zone):
all the grains have an elastic isotropic behavior but different Lamé coefficients.  
In the *material.xml* file below, the unique material (``numM=1``) is defined with constant per zone Lamé coefficients 
and the Lamé coefficients for the numerical integration of the averaging rule are the same. The introduction
of composite voxels is simply done by adding the section ``<Material_composite> .. <\Material_composite>`` to specify the directory
where the composite voxels are fully defined (see _`Format-for-composite-voxels-definition`)  :: 

	<?xml version="1.0" encoding="UTF-8"?>
	<Materials>
	
	<!-- REFERENCE MATERIAL -->
	<Reference_Material Lambda0="5.76923076923077e8" Mu0="3.84615384615385e8"/>
	
	<!-- MATERIAL 1  -->
	<Material numM="1" Lib="/home/gelebart/amitex_fftp/libAmitex/src/materiaux/libUmatAmitex.so" Law="elasiso" > 

			<Coeff Index="1" Type="Constant_Zone" File="materiaux/coefficients/Lambda1_polyx27G_R21.bin" Format="binary"/>
			<Coeff Index="2" Type="Constant_Zone" File="materiaux/coefficients/Mu1_polyx27G_R21.bin" Format="binary"/>

			<Coeff_composite Index="1" Type="Constant_Zone" File="materiaux/coefficients/Lambda1_polyx27G_R21.bin" Format="binary" />  
			<Coeff_composite Index="2" Type="Constant_Zone" File="materiaux/coefficients/Mu1_polyx27G_R21.bin" Format="binary"/>  
	</Material>

	<!-- DIRECTORY FOR THE DEFINITION OF COMPOSITE VOXELS -->
	<Material_composite> 
			<Coeff_composite directory="microstructures/voxcomp/polyX_27G_R21_reuss"/> 
	</Material_composite>

	</Materials>


**EXAMPLE 2 : a single composite voxel**

.. Note::

	In this exemple, and in any case where the definition of **"interphase" is required**, the simulation **can not be** run 
	without composite voxels by simply comenting the section ``<Material_composite   > .. </Material_composite>``.


In this example, the file *materialID.vtk* consists of a single voxel. 
The composite voxel definition (not described here) identifies this voxel as
a composite voxel between two different material (1 and 2). Obviously, 
in this specific case, there is no homogeneous voxel at all and these two
materials must be introduced as "interphase" materials.
The corresponding *material.xml* file is given below::  

	<?xml version="1.0" encoding="UTF-8"?>
	<Materials>
	
		<!-- REFERENCE MATERIAL -->
	    <Reference_Material Lambda0="22.21175e+10" Mu0="14.80765e+10"/>

	    <!-- MATERIAL 1 -->
	    <Material numM="1" Lib="(amitex)/libAmitex/src/materiaux/libUmatAmitex.so" Law="elasiso" >
			<Coeff Index="1" Type="Constant" Value="4.0385e+10"/>
			<Coeff Index="2" Type="Constant" Value="2.6923e+10"/>
			<Coeff_composite Index="1" Type="Constant" Value="4.0385e+10"/>
			<Coeff_composite Index="2" Type="Constant" Value="2.6923e+10"/>
	    </Material>
	    
	    <!-- MATERIAL 2 -->
	    <Material numM="2" Lib="(amitex)/libAmitex/src/materiaux/libUmatAmitex.so" Law="elasiso" >
			<Coeff Index="1" Type="Constant" Value="4.0385e+11"/>
			<Coeff Index="2" Type="Constant" Value="2.6923e+11"/>
			<Coeff_composite Index="1" Type="Constant" Value="4.0385e+11"/>
			<Coeff_composite Index="2" Type="Constant" Value="2.6923e+11"/>
	    </Material>

	    <!-- DIRECTORY FOR THE DEFINITION OF COMPOSITE VOXELS -->
	    <Material_composite>
			<Coeff_composite directory="microstructures/VoxComp"/> 
	    </Material_composite>

	    <!-- IDENTIFICATION OF MATERIALS WHICH ARE NOT PRESENT AS "HOMOGENEOUS" VOXELS -->
	    <Interphase>
			<Interphase_material  numM="1" Nzones="1"/>    
			<Interphase_material  numM="2" Nzones="1"/>
	    </Interphase>

	</Materials>

The two materials are desribed in the two sections ``<Material.. > .. </Material>``.
Note that Lamé coefficients are given, if necessary, for the numerical integration of the averaging rule
(``<Coeff_composite ../>``). 

As already mentionned, the two materials do not appear as homogeneous voxels and they are
identified as "interphase" ``<Interphase> .. </Interphase>``.



 
 


--------------------------------------------------------------------------------------------------------------------

.. _user-loading-outputs:

Loading and outputs: pure mechanical loading
------------------------------------------------

**EXAMPLE 1**

Here is an example of a *loading.xml* file where the loading and the stress-strain fields output are defined. This loading reproduces an experimental creep procedure : at first the load is applied proportionnaly (tensile test) and then the macroscopic stress is kept constant (creep test)::

	<?xml version="1.0" encoding="UTF-8"?>
	<Loading_Output>

	    <!-- ADDITIONNAL OUPUT QUANTITIES         -->
	    <!-- specify the quantities of interest   -->
	    <Output> 

	        <!-- Stress output (stress= "1") -->
		<vtk_StressStrain Strain = "0" Stress = "1"/>
		<vtk_IntVarList numM="1">
		  1
		</vtk_IntVarList>
		<Zone numM="1">
		    <VarIntList>  1  </VarIntList>
		</Zone>

	    </Output>

	    <!-- SUCCESSIVE LOADINGS AND OUTPUT TIMES -->

            <!-- Partial loading 1 -->
	    <Loading Tag="1">

		<!-- User defined time discretization, 10 increments -->
		<Time_Discretization Discretization="User" Nincr="10" />
		<!-- Increment times (user definition) -->
		<Time_List>
		32832 83808 143424 211680 289440 380160 483840 604800 738720 898560 
		</Time_List>
		<Output_zone Number ="10"/>
		<!--  No field output required --> 
		<!-- Strain driven on xx component --> 
		<xx Driving="Strain" Evolution="Linear" Value="-5e-4"/>
		<yy Driving="Stress" Evolution="Constant"  />
		<zz Driving="Stress" Evolution="Constant"  />
		<xy Driving="Stress" Evolution="Constant"  />
		<xz Driving="Stress" Evolution="Constant"  />
		<yz Driving="Stress" Evolution="Constant"  />
                <!-- For finite strains simulations       -->
		<yx Driving="Strain" Evolution="Constant"  />
		<zx Driving="Strain" Evolution="Constant"  />
		<zy Driving="Strain" Evolution="Constant"  />
	    </Loading>

            <!-- Partial loading 2 -->
	    <Loading Tag="2">

		<!-- Linear time discretization -->
		<Time_Discretization Discretization="Linear" Nincr="22" Tfinal="25920000" />
		<Output_zone Number ="5"/>
		<Output_cell Number="10"/> <!-- usefull for very high number of increments --> 
		<!-- Field output for the last increment (22)  -->
		<Output_vtkList>
		  22
		</Output_vtkList>
	        <!-- Full stress driving -->
		<xx Driving="Stress" Evolution="Constant" />
		<yy Driving="Stress" Evolution="Constant" />
		<zz Driving="Stress" Evolution="Constant" />
		<xy Driving="Stress" Evolution="Constant" />
		<xz Driving="Stress" Evolution="Constant" />
		<yz Driving="Stress" Evolution="Constant" />
                <!-- For finite strains simulations      -->
		<yx Driving="Strain" Evolution="Constant" />
		<zx Driving="Strain" Evolution="Constant" />
		<zy Driving="Strain" Evolution="Constant" />

	    </Loading>

	</Loading_Output>


.. Note::

	This *loading.xml* file example is designed for a finite strain simulation and must be consistent with `algorithm-parameters`_.

	In case of a Small Perturbation assumption, discard the ``yx``, ``zx`` and ``zy`` components.  
 

**EXAMPLE 2 : loading with a stress direction**

This second example imposes a ``xx`` strain loading while maintaining a constant stress direction given by DirStress. 
Here it corresponds to a constant triaxiality ratio of 2.5. ::
   
	<!-- OUPUT QUANTITIES -->
	<Output>
	    <vtk_StressStrain Strain = "0" Stress = "0"/>
	</Output>

	<!-- SUCCESSIVE LOADINGS AND OUTPUT TIMES -->
	<Loading Tag="1">
	    <Time_Discretization Discretization="Linear" Nincr="40" Tfinal="40"/>
	    <xx Driving="Stress" DirStress="1." />
	    <yy Driving="Stress" DirStress="1.3" />
	    <zz Driving="Strain" Evolution="Linear" Value="0.004" DirStress="1.6" />
	    <xy Driving="Stress" DirStress="0" />
	    <xz Driving="Stress" DirStress="0" />
	    <yz Driving="Stress" DirStress="0" />
	    <!-- below : Finite Strain specific -->     
	    <yx Driving="Stress" DirStress="0" />
	    <zx Driving="Stress" DirStress="0" />
	    <zy Driving="Stress" DirStress="0" />
	    <DirStress Type="cauchy" />  <!-- Cauchy or PK1 -->
	</Loading>



**BEGIN AND END**

The file must begin and end with ::

	<?xml version="1.0" encoding="UTF-8"?>
	<Loading_Output>
	.
	.
	</Loading_Output>


.. _user-output:

**DEFAULT OUTPUT**  

The section *<output> ... </Output>* should be optionnal. Up to now, the minimum requested is ::

	    <Output> 
		<vtk_StressStrain Strain = "0" Stress = "0"/>
	    </Output>

In that context, the default output reduces to three ouput files : 

	*.std* file : stress and strains average and standard deviations within the whole unit-cell. 
  
	*.mstd* file : per material stress and strains average and standard deviations. 

	*.log* file : informations on the simulation (number of iterations, criterion values, computation times etc...).
 
These quantities are output for each loading increment (= time step).


**SPECIFIC OUTPUT** 

The section *<output> ... </Output>* allows to specify two types of additionnal output : fields (*vtk* files) and per zone statistical quantities (*zstd*). Because these files can be very heavy (and their writing time consuming), the times to write these outputs files must be limited and chosen carefully. This is done in the section *<Loading> ... </Loading>*.

	* *.vtk files* - Strain and/or stress fields (VTK files) can be chosen through the ``<vtk_StressStrain>`` node. If no field are required, simply complete "0" and "0" for *Strain* and *Stress*. Similarly, the ``<vtk_IntVarList>`` node, gives the possibility to write in a VTK file the value of a list of internal variables related to the material *numM*. One ``<vtk_IntVarList>`` node per *material* node can be added, if necessary. The example below asks for *.vtk* output fields for the internal variables "1" and "5" of material "1" ::

		<vtk_IntVarList numM="1">
		  1 5
		</vtk_IntVarList>

	* *.zstd files* - If the *Zone* node is present, per zone average and standard deviations of the stress, strain and requested internal variables, are printed for the material associated to the *numM* value. If nothing else is precised, only the means and standard deviations of the stress and strain tensors will be printed. However, a ``<VarIntList>`` node can be added to choose the list of internal variables to be output. One ``<Zone>`` node per *material* can be added, if necessary. The example below asks for per zone output *.zstd* files for the stress, strain and internal variable "3" of material "2", and only for the stress and strain of material "1" ::  

		<Zone numM="1">
		</Zone>
		<Zone numM="2">
		    <VarIntList>  3  </VarIntList>
		</Zone>

The *.std*, *.mstd* and *.zstd* files are written so that the value of each variable is given in column. Comments at the beginning of these files describe the variable associated to each column.  More details about how to plot these results are given in the :ref:`plot` section.


**LOADING** 

The loading is defined by successive partial loadings, tagged 1, 2, 3, etc...

**Time_Discretization** within each partial loading can be chosen:

	* linear (loading 1 in the example): the user gives the final time and the number of increments, 

	* user defined (loading 2): the user gives a time list, excluding the initial time.

In both cases, the initial time is the final time of the previous loading or equal to 0 for the first partial loading (*i.e.* Tag="1").

**Output_vtkList** is the list of increment numbers, within the partial loading, for which stress, strain and/or internal variable fields (*.vtk* files) will be output (see the choice in the node ``<output>``).

**Output_zone** gives the number of loading steps in the partial loading for which per zone quantities (*.zstd* files) will be output. For each loading, the output increments are equally distributed among the increasing increments, including the final increment.
		
**Output_cell** gives the number of loading steps in the partial loading for which per unit-cell and per material quantities (*.std* and *.mstd* files) will be output. This is especially useful if the number of increment is very large.

**xx, yy, zz, xy, xz, yz** (and yx, zx, yz in finite strain) allows to define the loading, for each partial loading, and each component  :

	* Driving="Stress" or "Strain" if average stress or strain component is applied

	* Evolution="Linear" Value="val": the component evolves linearly as a function of the time from the initial value until "val". Voigt notation is assumed for small strains (the shear components of the strain tensor are multiplied by 2).

	* Evolution="Constant": the component is constant and equal to its initial value 

	* DirStress="val"(optional): if present must be present on all the components to impose a constant stress direction. At finite strains, an additional node must be given to specify wether the Cauchy or PK1 stress tensor is considered ::
	
		<DirStress Type="cauchy" />  <!-- Cauchy or PK1 -->

For "Constant" and "Linear" evolutions, the initial value is the final value of the previous partial loading, or equal to 0 for the first partial loading *i.e. Tag="1"*.


.. _warning-HPP:
.. warning::
	**With the small perturbation hypothesis** (see :ref:`algorithm-parameters`), the "Strain" and "Stress" components are related to the component of :

	* the average **linearized strain tensor**, 

	* the average **Cauchy stress tensor**. 

	Since the **Voigt notation** is used (in the inputs, the outputs and within the code itself), the value given for xy, xz or yz must be the double of the xy, xz or yz strain components.


	**Without the small perturbation hypothesis** (see :ref:`algorithm-parameters`), corresponding to the finite strain framework, the "Strain" and "Stress" components are related to the components of :

	* the average **first Piola-Kirchoff stress tensor**,

	* the average **displacement gradient** (which is not strictly speaking a measure of the strain).

--------------------------------------------------------------------------------------------------------------------

Mechanical and external loading
-------------------------------

If the user wishes to use a material behaviour law depending on the temperature and/or external parameters it has to be implemented in the **UMAT** compatible procedure (see `here <http://www-cast3m.cea.fr/index.php?page=sources&source=umat>`_). It is possible to impose the temperature and external parameters, constants on the whole unit-cell, during the loading. The temperature corresponds to the scalar **UMAT** variable *TEMP*. The other external parameters are gathered within the vector **UMAT** variable *PREDEF*. **These external parameters must be indexed continuously starting from 1.**

Here is an example of a *loading.xml* file where the temperature and two other external parameters are imposed::

	<?xml version="1.0" encoding="UTF-8"?>
	<Loading_Output>

	    <!-- ADDITIONNAL OUPUT QUANTITIES (optionnal) -->
	    <!-- <Output> </Output> -->

            <!-- Initialization of the temperature and the external parameters -->
            <InitLoadExt>
              <T Value="20" />
	      <Param Index ="1" Value="0" />
	      <Param Index ="2" Value="10" />
            </InitLoadExt>
	    <!-- SUCCESSIVE LOADINGS AND OUTPUT TIMES -->

            <!-- Partial loading 1 -->
	    <Loading Tag="1">

	        <Time_Discretization Discretization="Linear" Nincr="40" Tfinal="40"/>
		<xx Driving="Strain" Evolution="Linear" Value="-5e-4"/>
		<yy Driving="Stress" Evolution="Linear" Value="0." />
		<zz Driving="Stress" Evolution="Linear" Value="0." />
		<xy Driving="Stress" Evolution="Linear" Value="0" />
		<xz Driving="Stress" Evolution="Linear" Value="0" />
		<yz Driving="Stress" Evolution="Linear" Value="0" />

		<!-- Temperature is constant, equal to the initial value -->
		<T Evolution="Constant" />
		<!-- External parameters 1 and 2 have a linear evolution -->	
		<Param Index="1" Evolution="Linear" Value ="50"/>
		<Param Index="2" Evolution="Linear" Value ="0"/>

	    </Loading>

	</Loading_Output>


In order to be imposed during the loading the temperature and each external parameter **must** :

1. **be initialized** inside the section *<InitLoadExt> ... </InitLoadExt>* ::

	<InitLoadExt>
              <T Value="20" />
	      <Param Index ="1" Value="0" />
	      <Param Index ="2" Value="10" />
	</InitLoadExt>


2. be described **within each section** *<Loading Tag="i">* as described in the example. The *Evolution* can be linear or constant as a function of time. In both cases, the initial value is the final value of the previous partial loading or equal to the value defined in the *InitLoadExt* tag for the first partial loading (*i.e.* Tag="1"). ::

	<!-- Temperature is constant, equal to the initial value -->
	<T Evolution="Constant" />
	<!-- External parameters 1 and 2 have a linear evolution -->	
	<Param Index="1" Evolution="Linear" Value ="50"/>
	<Param Index="2" Evolution="Linear" Value ="0"/>
		

--------------------------------------------------------------------------------------------------------------------

.. _algorithm-parameters:

Algorithm parameters
---------------------

This section describes the algorithm parameters given the file *algorithm.xml* (see :ref:`how-to-launch`).

**DEFAULT ALGORITHM FILE**

The default algoithm file given below should be used first and modified if needed  ::

	<?xml version="1.0" encoding="UTF-8"?>
	<Algorithm_Parameters>

   	<Algorithm Type="Basic_Scheme">                     <!-- "Default" (Basic_Scheme) -->
            <Convergence_Criterion Value="Default"/>        <!-- "Default" (1e-4) or positive value <1e-3 -->
            <Convergence_Acceleration Value="true"/>        <!-- "True" or "False" -->
    	</Algorithm>

    	<Mechanics>
            <Filter Type="Default"/>                        <!-- "Default" (hexa) ou "no_filter" ou "hexa" ou "octa" -->
            <Small_Perturbations Value="true"/>             <!-- "True" ou "False" -->
    	</Mechanics>

	</Algorithm_Parameters>



**Filter** - Different discrete Green operators are implemented in the code :
	
	* the classical one [Moulinec1998]_
	* the hexaedral and octaedral filtered Green operators, which are Finite Difference based operators

The hexaedral filtered Green Operator is the "Default". Equivalent to the so-called "rotated" Green operator proposed by Willot [Willot2015]_, it is :
	
	* very efficient to reduce the sensitivity to the mechanical contrast,
	* efficient to reduce spurious oscillations


**Small_Perturbation** - The user can choose to perform the simulation with or without the small perturbation assumption. Of course, the loading file must be consistent with this choice (see :ref:`Warning <warning-HPP>`).


**Algorithm** - Here, the basic scheme [Moulinec1998]_ is implemented with an additionnal convergence acceleration procedure which can be used or not. This convergence acceleration procedure (see procedure *ACT3* in the finite element code CAST3M), described in [Chen2019]_, is :

	* very efficient to reduce the sensitivity to the mechanical contrast,
	* very efficient to reduce the sensitivity to the reference material choice,
	* more memory consuming.

The convergence criterion is given by an equilibrium condition (*i.e.* :math:`div(\sigma) = 0`) together with an average condition on the applied load. The default value is :math:`10^{-4}` (a maximum value is set to :math:`10^{-3}`).

**If memory is not a problem, do not hesitate and use the convergence acceleration!** 


**ADDITIONAL PARAMETERS**

By default, choice is made in *amitex_fftp* to perform at least one iteration (even if the initial guess satisfies convergence criteria), and at least one iteration after a convergence acceleration. This is not mandatory, and it can be overcome this by setting **Nitermin** and **Nitermin_acv** to 0 as follows ::  

   	<Algorithm Type="Basic_Scheme">                     
             <Convergence_Criterion Value="Default"/>
             <Convergence_Acceleration Value="true"/>
             <Nitermin Value="0"/>            <!-- "0" or "1" : 1 = at least one iteration per loading increment-->
             <Nitermin_acv Value="0"/>        <!-- "0" or "1" : 1 = at least one iteration after each accelerated solution-->
    	</Algorithm>

Nitermin=Nitermin_acv=0 will probably become the default behavior in the future.

**In case of convergence issue**, and only in that case, the following parameters can be added ::

   	<Algorithm Type="Basic_Scheme">                     
             <Convergence_Criterion Value="Default"/>
             <Convergence_Acceleration Value="true"/>
             <Convergence_Criterion_Smacro Value="1e-3"/>      <!--  positive value <1e-1 -->
                                                               <!--  default = value used for Convergence_Criterion -->
             <Convergence_Criterion_Compatibility Value="1e-9"/>    <!-- "Default" (1e-10) or positive value <1e-3 -->
             <Nitermax Value="2000"/>                               <!-- Default is 1000 -->
             <Initialize Value="previous"/>     <!-- "default" or "previous" -->
                                                <!-- if previous : uses the solution field for the previous step 
                                                                   as the initial guess for the current step -->
    	</Algorithm>

**Convergence_Criterion_Smacro**, the convergence criterion on the macroscropic applied stresses can be relaxed (it can be usefull when applying a loading with a loading with an applied stress direction).

**Convergence_Criterion_Compatibility**, the compatibility criterion fixed at :math:`10^{-10}`, can be relaxed. Note that this value must remain very small because the strain field being an output of the Green operator it shoud be naturally compatible. As an example, it has been necessary, in case of a crack simulation to reduce this value to :math:`10^{-9}`.

**Nitermax**, the maximum number of iterations can be increased.

**Initialize** can be set to "previous" so that the strain field initial guess at the begining of a new step is the strainfield obtained at the end of the previous step. Instead, the initial guess is by default a linear extrapolation of the two previous strain fields. 

**Finite Strains with non-symmetric reference material behavior**

In the finite strain frameworks, the reference material behavior (relating the first Piola Kirchoff stress tensor  to the displacement gradient) can be symmetrized or not. By default in *amitex_fftp*, it is symmetrized. This choice results from simulations with small strain values, leading to identical convergence properties for the finite strains and the small strain frameworks (this is no more the case when using a non-symmetrized behavior). However, non-symmetrized reference behavior can be used as follows ::

	<Mechanics>
	    <Filter Type="Default"/>                                         
	    <C0sym Value="false"/>      <!-- default is true -->
	</Mechanics>

**Small strain simulation with non symmetrized displacement gradients**

It is possible within *amitex_fftp*, to perform simulations with the small strains assumption but taking into account the complete displacement gradient (and not only the symmetric part). This allows to use behaviors that depend on the non symmetric part of the displacement gradient through the variables DFGRD0 and DGRD1 of the *umat* procedure. This option can be used as follows ::

	<Mechanics>
	    <Filter Type="Default"/>                                         
	    <Small_Perturbations Value="true" Displacement_Gradient="nsym"/> 
	</Mechanics>

Note that up to now, the average applied displacement gradient is assumed symmetric with the strain components given in the *loading.xml* file.


--------------------------------------------------------------------------------------------------------------------

.. _user-behavior-label:

User behavior law 
------------------

General case
^^^^^^^^^^^^

**amitex_fftp** offers the possibility to introduce user material behaviors through a *UMAT* compatible procedure.

An example is available in *(amitex_fftp)/cas_tests/comportements/polyxCC/comportement_umat*. The user can copy-paste this folder on his home directory (denoted here as *home*) and follow the procedure below.

The folder contains user material files and an additionnal *Makefile*. The name of the *umat* compatible procedure used for the integration of the behaviour law corresponds to the field ``Law`` to be filled in the *material.xml* file (see :ref:`Material-properties`).

The *Makefile* launch the compilation and create in the current folder a dynamic library *libUmatAmitex.so*. 
The full path towards this library must be used to fill the ``Lib`` field in the *material.xml* file (see :ref:`Material-properties`).

1. Adjust the following first lines of the Makefile ::

	# Compiler (ifort or gfortran)
	FC=ifort

	# Path to the librairie 2decomp_fft, (amitex)/lib_extern/2decomp_fft-1.5.847, for example:
	lDecomp=/home/gelebart/amitex_fftp/lib_extern/2decomp_fft-1.5.847

2. Launch the compilation to create the dynamic library ::
	
	make clean
	make

3. Adjust the ``Law`` and ``Lib`` fields of the desired *material.xml* files.

Using MFront
^^^^^^^^^^^^

**Interfacing AMITEX_FFTP with MFront libraries is straightforward.**

Of course, it assumes that you have installed MFront on your machine and that the LD_LIBRARY_PATH contains the MFront library path.

An example with a dynamic library created by MFront, specifying the *umat* interface , is given in *(amitex_fftp)/cas_tests/mazars*.

The dynamic library *(amitex_fftp)/cas_tests/mazars/src/libUmatBehaviour.so* is an output of MFront using the command: ::

     mfront --obuild --interface=umat mazars.mfront

where *mazars.mfront* is the Mazars behavior law [Mazars1990]_ delivered with MFront (see http://tfel.sourceforge.net/documentations.html for more information).

The new behaviour law can be used by filling correctly the ``Law`` and ``Lib`` fields of the desired *material.xml* files.

.. Note ::
	
	The ordering of the coefficients and internal variables in the *material.xml* file corresponds to their ordering in the *.dgibi* file (CAST3M finite element code input) located in the local *castem* folder generated together with the *src* folder when launching ``mfront --obuild --interface=umat file.mfront``.   

--------------------------------------------------------------------------------------------------------------------

.. _user output&algorithm:

User output & algorithm 
-------------------------

In the program, a standard output file *output.std* is printed after each loading increment. However, the user has also the possibility to modify this output through the *sortie_std* procedure.

As mentionned in section `algorithm-parameters`_, the user has the possibility to implement its own user algorithm.
 
The purpose of this section is to show how to build a new executable program *AMITEX_FFTP* incorporating a new user standard output and a new user algorithm.

A complete example, is available in *(amitex_fftp)/cas_tests/comportements/polyxCC*. The user can copy-paste this directory on his home directory (denoted here as *home*) and follow the procedure below.

The user output file and an additionnal *Makefile* are located in *home/polyxCC/sortie_std*.
The user algorithm file and an additionnal *Makefile* are located in *home/polyxCC/user_algo*.

The *Makefile* located in *home/polyxCC* has to be correctly filled with the correct path to libraries.

The executation of the *Makefile* creates a new executable *amitex_fftp* in the current folder.

--------------------------------------------------------------------------------------------------------------------

.. _Diffusion:

Diffusion 
---------

Since version 4.0.0, one variable stationnary diffusion problems, such as thermal diffusion, can be solved with **amitex_fftp**.

As the concepts used to describe the inputs of the code in mechanics are similar for diffusion, we only explain below some typical *xml* files.  


Material properties 
^^^^^^^^^^^^^^^^^^^

::

	<?xml version="1.0" encoding="UTF-8"?>
	<Materials>

	<!-- REFERENCE MATERIAL -->
	<Reference_MaterialD K0="14"/>

	<!-- MATERIAL 1 -->
	<Material numM="1" LibK="/home/gelebart/amitex_fftp/libAmitex/src/materiauxK/libUmatAmitexK.so" LawK="Fourier_iso">
                <CoeffK Index="1" Type="Constant_Zone" File="materiaux/coefficients/coeffK.txt" Format="ASCII"/>
	</Material>

	</Materials>

After defining the reference material property ``<Reference_MaterialD K0="14"/>``, the material behavior is given by the fields ``LibK`` and ``LawK``.

Here the simple isotropic Fourier behaviour is used and its coefficients are constant per zone, given in an ASCII file.


Algorithm parameters
^^^^^^^^^^^^^^^^^^^^

::

	<?xml version="1.0" encoding="UTF-8"?>
	<Algorithm_Parameters>

	<!-- GENERAL CASE  --> 
   	<Algorithm Type="Basic_Scheme">                         <!-- "Default" (Basic_Scheme) or "user" -->
       		<Convergence_Criterion Value="Default"/>        <!-- "Default" (1e-4) or positive value <1e-3 -->
       		<Convergence_Acceleration Value="True"/>        <!-- "True" ou "False" -->
   	</Algorithm>

	<!-- DIFFUSION  --> 
   	<Diffusion>    
       		<Filter Type="Default"/>                         <!-- "Default" (hexa) or "no_filter" or "hexa" or "octa" -->
       		<Stationary Value="true"/>                       <!-- "Tue" ("False" not implemented yet) -->
   	</Diffusion>
	</Algorithm_Parameters>


No comment, reading is enough.


Loading and output
^^^^^^^^^^^^^^^^^^

::

	<?xml version="1.0" encoding="UTF-8"?>
	<Loading_Output>

	<!-- OUPUT QUANTITIES -->
	<Output>
	    <vtk_FluxDGradD FluxD = "1" GradD="1"/>  
	    <Zone numM="1">
	    </Zone>
	</Output>
	
	<!-- SUCCESSIVE LOADINGS AND OUTPUT TIMES -->
	<Loading Tag="1">
	    <Time_Discretization Discretization="User" Nincr="1" />
	    <Time_List>1</Time_List>
	    <Output_vtkList>1</Output_vtkList>
	    <Output_zone Number="1"/>
	    
	    <x0 Driving="GradD" Evolution="Linear" Value="0.1"/>
	    <y0 Driving="FluxD" Evolution="Linear" Value="2" />
	    <z0 Driving="FluxD" Evolution="Linear" Value="3" />
	
	</Loading>

	</Loading_Output>
 

Here, the loading is made in one increment with a prescribed average temperature gradient in *x* direction and a prescribed average flux in *y* and *z* directions.

The **default output** are the same as in Mechanics, replacing the stress and strain tensors by the flux and temperature gradient vectors.

In addition, here, the **specific outputs** are : the flux and temperature gradient fields (vtk files), per zone statistical quantities for temperature gradient and flux.   


