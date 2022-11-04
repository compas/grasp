
=======================================================================

     Coupling:
     The program for searching optimal coupling scheme in atomic theory

     Gediminas Gaigalas
     Institute of Theoretical Physics and Astronomy,
     Vilnius University, Sauletekio ave. 3, LT-10222 Vilnius,
     Lithuania

     e-mail: gediminas.gaigalas@tfai.vu.lt

=======================================================================

The program is able to perform transformation from LS-coupling to other couplings
independently from other codes. It can also be included in the sequence of calculations
of Grasp2018 or Atsp2K packages.


===========================
INSTALLATION OF THE PROGRAM
===========================

Coupling program has been developed as (a new) component of the Grasp2018 package.
The steps below should be followed to ensure a proper installation in the bash shell.
The installation procedure assumes that the Grasp2018 package is already installed.

1. Go to the main directory GRASP2018 of the Grasp2018 package. Type

   >> source make-environment_xxx

   where xxx is the compiler name (see the README file of the Grasp2018 package.
   The Grasp2018 environment variables are now set.

2. Copy the file Coupling.tar.gz to the directory GRASP2018/src/appl of the Grasp2018
   package. Untar it by typing

   >> tar -zxvf Coupling.tar.gz

   A directory Coupling will appear.

3. In the grasp2018/src/appl/Coupling directory, execute the installation by typing

   >> make -f make-Coupling clean
   >> make -f make-Coupling

   This will generate an executable file Coupling in the directory GRASP2018/bin.


To ensure that the Coupling program is fully incorporated in future recompilations
of the entire grasp2K package do the following:

Go to the directory GRASP2018/src/appl and add COUPLING to the variable SUBDIR in
the Makefile.

It is also recommended that the user includes the path to the bin directory of the
Grasp2018 package in the system variable PATH to facilitate execution.


========================
EXECUTION OF THE PROGRAM
========================

The program is executed by typing

>> GRASP2018/bin/Coupling

where GRASP2018 denotes the full path to the Grasp2018 main directory, or simply by typing

>> Coupling

if the path to the directory GRASP2018/bin is included in the system variable PATH.

The subdirectory GRASP2018/src/appl/Coupling/Sample_Runs lists a number of examples
demonstrating the usage of the program.

The subdirectory GRASP2018/src/appl/Coupling/Sample_Output contains output files from this
three examples. To validate program operations these output files can be as references.


=====================
THE PROGRAM STRUCTURE
=====================

The program is written in Fortran 90 programming language. The main part of the program
consists of several subroutines collected into the module Coupling (file Coupling.f90).
The program's design is similar to the jj2lsj program from the Grasp2018 package.

The program itself has 32 separate modules:

   Coupling_main,
   Coupling_structures,
   Coupling_constants,
   Coupling_data,
   Coupling_inside_shell,
   Coupling_getdata_mchf,
   Coupling_evaluation,
   Coupling_transform_cfg_JJJK,
   Coupling_transform_cfg_LScLSJ3,
   Coupling_transform_cfg_LSjj1,
   Coupling_transform_cfg_LSjj2,
   Coupling_transform_cfg_LSjj3,
   Coupling_transform_cfg_LSJJ,
   Coupling_transform_cfg_LSJK3,
   Coupling_transform_cfg_LSLK3,
   Coupling_transform_cfg_LSLK,
   Coupling_transform_cfg_LSLS3,
   Coupling_transform_cfg_LSLScjj,
   Coupling_transform_cfg_LSLSJ3,
   Coupling_transform_JJJK,
   Coupling_transform_LScLSJ3,
   Coupling_transform_LSjj1,
   Coupling_transform_LSjj2,
   Coupling_transform_LSjj3,
   Coupling_transform_LSJJ,
   Coupling_transform_LSJK3,
   Coupling_transform_LSLK3,
   Coupling_transform_LSLK,
   Coupling_transform_LSLS3,
   Coupling_transform_LSLScjj,
   Coupling_transform_LSLSJ3,
   sqlsf.

Three modules from the jj2lsj program from GRASP2018 are also used:

   jj2lsj_data_1_C,
   jj2lsj_data_2_C,
   jj2lsj_data_3_C.


=======================================================================
THE SCHEME OF USE OF COUPLING IN THE SEQUENCE OF GRASP2018 CALCULATIONS
=======================================================================

>> rnucleus          # Generation of nuclear data

>> rcsfgenerate      # Generation of list of CSFs based on rules
                     # for excitations

>> rcsfinteract      # Reduction of a list to CSFs interacting with
                     # the multireference

>> rangular          # Angular integration

>> rwfnestimate      # Initial estimates of radial orbitals

>> rmcdhf            # Self-consistent field procedure

>> rci               # Relativistic RCI with optional transverse photon
                     # (Breit) interaction and vacuum polarization and
                     # self-energy (QED) corrections

>> jj2lsj            # Transform representation from jjto LSJ-coupling

>> Coupling          # Present atomic state functions in different
                     # coupling schemes

>> rhfs              # Hyperfine structure calculation

>> rhfsd             # Hyperfine structure calculation

>> ris               # Isotop shift calculation

>> rbiotransform     # Biorthonormal transformation

>> rtransition       # E1, M1, E2 transition calculations


================
ACKNOWLEDGEMENTS
================

The author wishes to thank Dr. Alexander Kramida and the Atomic Spectroscopy
Group of the National Institute of Standards and Technology, USA for their
support and encouragement. Part of the work was funded under a Guest Researcher
Agreement G-3-00334 at NIST.


===========
MIT LICENSE
===========
