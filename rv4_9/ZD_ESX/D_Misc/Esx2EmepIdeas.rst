
EMEP + ESX ?
============

General idea:
-------------

Keep esx as a parallel system, with own config_esx settings. 

Towards end of DryDep  we have all the concentations (xn_2d) and meteo we should be in the Sub(il)=L types.

if ( uses%esx  )

   Take xn_2d  (conc., mol/cm3) from k=20=KMAX_MID (bottom) EMEP into esx xChem arrays::

       e.g.  xChem(SPEC, :) = xn_2d(SPEC, kESX: 20)


   Take meteo from L or Sub(iL) into Zmet::

       e.g. Zmet(:)%rh =  L%rh




Common routines:
----------------

   Easy? ::

     PhysicalConstants

     ModelConstants

   Moderate work? ::

     Solver

     ChemFunctions

   More work ::

     Kz  (needs hacking of EMEP BLPhysics)

