GenChem
=======

Ideas for perl -> python
------------------------

The heart of GenChem.pl are the calls to::

   read_species()

   read_shorthand()

   read_reactions()


then various printouts to modules. I suggest simply making python equivalents for each, in turn. Thus, start with read_species, since this actually does a lot of useful work. I  have made some documentation of this below, but basically this routine reads the GenIn.species file, and assigns species to various characteristic types and groups.  The main types are as specified in the "adv" column, and this is used to assign species to the EMEP model's short-lived (adv=0 -> SHL) or advected species (adv=1 or 2 -> ADV) sets, with e.g. O3 ending up with an index IXADV_O3 within the advected set, and simply O3 within the total indices. A final main type is for organic aerosol (adv=2), which is both advected, and given some special treatments in the EMEP code.

In addition, the column for groups allows species to be allocated to groups, e.g. NO, NO2, HNO3 etc. all have the OXN group, and where deposited this results in a DRYDEP_OXN group. See the output CM_ChemGroups file mentioned below to see such how these groupins look.


The main first job though is just  to parse GenIn.species, and start collecting names ands groups in appropriate arrays. GenChem.pl also does nice tricks to calculate mol. wts if the chemical formula is given.


read_species - changes in python?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1) Changed integer to text flag for adv

GenChem currently assumes that first we have  short-lived species, then longer-lived, so it gets confused if the user starts to add new sets of short-lived later. The new code should allow for species in any order. Also, I am tempted to replace the current use of adv=0, adv=1, etc. by more explicit characters, e.g. "SHL" instead of adv=0, "ADV" instead of adv=1. The current use of adv=2 for SOA species can be replaced by using something in the group column
instead, I can work that out later.

 
2) Remove use of "SLOW"

Currently GenChem scans the GenIn.species file for the text marker SLOW, and
produces one CM_Reactions1.inc file for species above that, and another
CM_Reactions2.inc for the species below. I suggest we use the new "adv" type to indicate slowness, but still keep the CM_Reactions2.inc outputs.

This would again allow species files to be easily added (or 'catted' to one another, without regard to the position of short-lived, slow, or other types of species).


Background
----------

Usually, chemists write chemical equations in forms such as::

    NO  +  O   ->    NO2                   rate k1
    NO2 +  OH  ->    HNO3                  rate k2
    NO2        ->    LOSS                  rate Vg

where the first reaction has two reactants which react at rate k1 to produce
the product NO2, similarly for the 2nd equation. The third equation has
just one reactant, being lost at some rate Vg (deposition in this case).
For use in chemical transport models, these types of chemical equation need
to be turned intio different equations of type::

    d Ci/ dt = Pi - Li . Ci

where Ci is the concentration  of species i, Pi and Li are the production and
loss rates. Thus::

   d NO2/dt =  k1. NO.O - k2.OH.NO2 - Vg.NO2

Or, expressed as production and loss::

    P(NO2) = k1.NO.O
    L(NO2) = k2.OH + Vg


The aim of GenChem is to convert chemical equations into differential ones. GenChem also produces various modules to help the EMEP model. The three files read in by GenChem are:

	*  GenIn.species
	*  GenIn.shorthands
	*  GenIn.reactions

These files, and their processing are discussed below.


GenIn.species
-------------

This file lists the chemical species, and their properties. The file format is excel/gnumeric friendly, something like::

   Spec adv formula              MW DRY    WET   Extinc Cstar DeltaH Empty Groups  !Comments
   ---- --- -------------------- -- -----  ----- ------ ----- ------ ----- ------  ---------
   OD   0   O                    xx xx    xx    0      0     0       xx    xx      !
   HNO3 1   HNO3                 xx HNO3  HNO3  0      0     0       xx    OXN     !
   MPAN 1   CH2CH(CH3)C(=0)O2NO2 xx PAN   xx    0      0     0       xx    OXN     !( from macr degredation)
   ---- --- -------------------- -- -----  ----- ------ ----- ------ ----- ------  ---------

(File really has commas for separator.)

There are also some header lines which are worth reading. The columns are::

	Spec	species name in model

        adv	0 for short-lived, 

		1 for advected

		2 if SOA

		3 if slowly changing, see CODE-change comments below.

	formula	chemical formula (optional). 

		NOTE: if upper-case, formula assumed! 

	MW	mol. wt. (not needed if formula specified, or not emitted)

	DRY	surrogate species used for dry-deposition

	WET	surrogate species used for wet-deposition

	Extinc	extinction coefficient

	Cstar	saturation vapour pressure for organic aerosol (OA)

	DeltaH  Enthalpy (for OA)

	Empty 	Not used

	Groups	Very important! For example, DRYDEP, OXN


Thus, we might want to say that the species MPAN has the same dry deposition as PAN (thus PAN in DRY column), and belongs to the oxidised nitrogen group, OXN. 
Deposition is a little special in that GenChem also keep tracks of this, making
groups for e.g. DRYDEP_OXN which would


GenChem reads the file, and produces from this::

 CM_ChemDims.f90	- basic dimensions

 CM_ChemSpecs.f90	- assigns each species to arrays for short-lived, advected
			  or SOA (according to "adv" column), calculates
			  molecular weights, and assigns groups


 CM_ChemGroups.f90	- groups of species, e.g. collections of BVOC or OXN.

 CM_DryDep.inc		- mapping of real species to surrogate specified by DRY column 

 CM_WetDep.inc		- mapping of real species to surrogate specified by WET column 



        
CODE change. GenChem.pl has treatment for slow species by use of the #SLOW marker in the .csv file. Don't bother with this, but use the number 3 for adv instead.
