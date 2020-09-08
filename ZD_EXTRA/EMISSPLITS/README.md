June 2020 revised emissplits::

    - uses new .csv style (libreoffice friendly!)
    - Ref1 and Ref2 versions, to cope with condensables in PM splits
    - GNFR_CAMS (19 sector splits) provided here
    - SNAP-11 created from these GNFR_CAMS files.
    - revised NMVOC from CAMS71/Robert/Dave 
    - now part of emep-dev, in ZD_EXTRA/EMISSPLITS, but should be copied
      to DataDir/ZCMDIRS/EMISSPLITS (if not already there).
     

Directories::

  EMISPLITS/               - new central directory for emissplits::

     EmChem19/   - ok for EmChem19a, EmChem19p 

    /emissplits_gnfr   - updated system with new .csv style
                       - main directory for GNFR_CAMS 13/19-sector approach

         CAMS-REG-AP_v2.2.1_2015_REF1/ - based upon TNO's Ref1 emissions

         CAMS-REG-AP_v2.2.1_2015_REF2/ - based upon TNO's Ref2 emissions
  
         The only difference between REF1 and REF2 is in the PM splits, but
         complete sets of emissplits are given to make it easy in config
         setup. The NMVOC splits are identical for both, and taken from
         the CAMS71-II-scenarios_v1 data set which gave F1, F2, etc.

    /emissplits_snap   - updated system with new .csv style
                       - as gnfr above, but into SNAP-11
  
    
For sox, nox, co, nh3 the CAMS_GNFR files were simply reformated from the
earlier files.
