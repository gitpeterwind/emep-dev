GenChem
=======


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


This file lists the chemical species, and their properties. The file format is excel/gnumeric friendly:

   Spec,adv,forumula,MW,DRY,WET,Extinc,Cstar,DeltaH,Empty,Groups,,!Comments
   OD,0,O,xx,xx,xx,0,0,0,xx,xx,,!
   OP,0,O,xx,xx,xx,0,0,0,xx,xx,,!
   O3,1,O3,xx,O3,xx,0,0,0,xx,OX,,!
   NO,1,NO,xx,xx,xx,0,0,0,xx,NOX;OXN,,!
   NO2,1,NO2,xx,NO2,xx,0,0,0,xx,NOX;OX;OXN;daObs,,!
   PAN,1,CH3COO2NO2,xx,PAN,xx,0,0,0,xx,OXN,,!
   MPAN,1,CH2CH(CH3)C(=0)O2NO2,xx,PAN,xx,0,0,0,xx,OXN,,"!( from macr degredation)"
   NO3,1,NO3,xx,xx,xx,0,0,0,xx,OXN,,!
   N2O5,1,N2O5,xx,xx,xx,0,0,0,xx,OXN,,!
   ISONO3,1,isono3,xx,xx,xx,0,0,0,xx,OXN,,"!isoprene-NO3 adduct"
   HNO3,1,HNO3,xx,HNO3,HNO3,0,0,0,xx,OXN,,!

(Looks nicer in excel!)

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
