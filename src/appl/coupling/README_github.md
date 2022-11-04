## Coupling - A program for searching optimal coupling scheme in atomic theory

Note: This version is an improved, further developed version of the code published in CPC.

### Introduction
The Coupling program, which is important not only for the Grasp2018 package but for the atom theory in general, is presented in this paper. This program is designed as a part of the Grasp2018 package. The Coupling performs the transformation of atomic state functions (ASFs) from a LSJ-coupled CSF basis into several other configuration state function (CSF) bases such as jj, LK, JK, JJ, and many others. It allows identification of the energy structure of practically any element in different coupling schemes and also allows selection of the most suitable one. In addition, examples of how to use the Coupling program are given in additional file Example_Calculation.pdf listed in source directory grasp2018/src/appl/Coupling/Sample_Runs.

### Program summary
   - Program Title: Coupling
   - Program Files doi: http://dx.doi.org/10.17632/7tv9n8g24w.1
   - Licensing provisions: MIT license
   - Programming language: Fortran 95
   - External routines/libraries used: Grasp2018 modules: jj2lsj_data_1_C, jj2lsj_data_2_C, jj2lsj_data_3_C and Atsp2K module: sqlsf.f.

### Nature of problem
The Coupling program is designed as a part of the Grasp2018 package for the computation of atomic state function transitions and for identification atomic properties in different coupling schemes.
Solution method: Spin-angular coupling transformation of multiconfiguration expansions obtained in LS and jj coupling schemes.
Additional comments including restrictions and unusual features: The transformations of the configuration state functions is supported for one, two, and three open shell structures including open s, p, d, and f-shells. For shells with l > 3 (i.e. beyond the f-shells), however, a proper transformation of the antisymmetrized shell states can be carried out only for the case of one or two equivalent electrons. The jj ←→ LS transformation matrices, which are applied internally by the program, are consistent with the previously published definitions of the reduced coefficients of fractional parentage and jj ←→ LS transformation matrices. The transformation from LSJ-coupling to all main coupling schemes is performed by the program Coupling. The number and types of coupling schemes depend on the number of coupled shells in LS-coupling.
