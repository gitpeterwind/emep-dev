Using Waf
=========

Waf_ is a Python-based standalone build system framework.  The build system itself is contained in 
the ``waf`` file, and its configuration lives in the ``wscript`` file, which is currently configured 
to compile all ``.f90`` files in this directory into one executable named ``esx_tester``.

To build and run the program using Waf::

    ./waf configure         # (or: python waf configure) Detect/configure environment
    ./waf build             # Build program
    ./build/esx_tester      # Run program

To use a specific compiler executable, the ``FC`` environment variable is obeyed::

    FC=ifort ./waf configure
    # or...
    FC=/usr/local/bin/gfortran-4.8 ./waf configure
    ./waf build
    ./build/esx_tester


.. _Waf: https://code.google.com/p/waf/
