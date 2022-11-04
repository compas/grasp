!**********************************************************************************
!**********************************************************************************
!**                                                                              ** 
!**  **  ** *****  ****   ******* ***** ***** **     **       *       **     **  **
!**  **  ** **    **  **      **  **    **    ***   ***      ***      ***    **  **
!**  **  ** **    **         **   **    **    ** * * **     ** **     ** *   **  **
!**  ****** ****   ****     **    ****  ****  **  *  **    **   **    **  *  **  **
!**  **  ** **        **   **     **    **    **     **   *********   **   * **  **
!**  **  ** **    **  **  **      **    **    **     **  **       **  **    ***  **
!**  **  ** **     ****  *******  ***** ***** **     ** **         ** **     **  **
!**                                                                              **
!**                                                                              **
!**                  Relativistic Hyperfine and Zeeman Program                   **
!**                               GRASP2K Version                                **
!**                               Dynamic Storage                                **
!**                                                                              **
!**  =========================================================================   **
!**  Copyright (c) 2006 by M Andersson, P Jonsson                                **
!**  =========================================================================   **
!**   All rights reserved.  No part of this software or its accompanying docu-   **
!**   mentation may be reproduced,  stored in a retrieval system, or transmit-   **
!**   ted, in any form or by any  means, electronic, mechanical, photocopying,   **
!**   recording,  or otherwise,  without the  prior written  permission of the   **
!**   authors.                                                                   **
!**                                                                              **
!**                                Disclaimer                                    **
!**                                ==========                                    **
!**   The authors make no warranties,  express or implied,  that this software   **
!**   or its accompanying  documentation are  free of error or that  they will   **
!**   meet your  requirements  for any  particular  application.  Neither this   **
!**   software  nor its accompanying  documentation should  be relied upon for   **
!**   solving problems if an incorrect solution could result in injury or loss   **
!**   of, or damage to, property. If you use this software or its accompanying   **
!**   documentation, you do so entirely at your own risk; the authors disclaim   **
!**   all liability for direct or consequential damage.                          **
!**                                                                              **
!**********************************************************************************
!**********************************************************************************
!*                                                                                *
      PROGRAM HFSZEEMAN06
!*                                                                                *
!*   Entry routine for HFSZEEMAN. Controls the entire computation.                *
!*                                                                                *
!*   Call(s) to: [LIB92]: GETMIXA, GETMIXC, SETCSL, SETMC, SETCON.                *
!*               [HFSZEEMAN06]: CHKPLT, GETHFD, GETMIXBLOCK, HFSZEEMAN, SETDBG,   *
!*                               SETSUM, STRSUM.                                  *
!*               [NJGRAF]: FACTT.                                                 *
!*                                                                                *
!*      M. Andersson and P Jönsson                               2006             *
!*                                                                                *
!*   Translated by Wenxian Li F77 to F90 12/28/18                                 *
!*********************************************************************************
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE default_C
      USE iounit_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s 
!-----------------------------------------------
      USE getyn_I 
      USE setdbg_I 
      USE setmc_I 
      USE setcon_I 
      USE setsum_I 
      USE setcsla_I 
      USE gethfd_I 
      USE getmixblock_I 
      USE strsum_I 
      USE factt_I 
      USE hfszeeman_I 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s 
!-----------------------------------------------
      INTEGER   :: K, NCI, NCORE_NOT_USED, NOFFD 
      LOGICAL   :: YES 
      CHARACTER :: NAME*24 
!-----------------------------------------------

!
      WRITE (ISTDE, *)
      WRITE (ISTDE, *) 'HFSZEEMAN95'
      WRITE (ISTDE, *) 'This is the magnetic interaction program'
      WRITE (ISTDE, *) 'Input files:  isodata, name.c, name.(c)m, name.w'
      WRITE (ISTDE, *) 'Output files: name.(c)h, name.(c)gjhfs'
      WRITE (ISTDE,*)
      WRITE (ISTDE,*) 'HFSZEEMAN95: Execution begins ...'

      WRITE (ISTDE,*)
      WRITE (ISTDE,*) 'Default settings?'
      YES = GETYN ()
      WRITE (ISTDE,*)
      IF (YES) THEN
         NDEF = 0
      ELSE
         NDEF = 1
      ENDIF

   10 CONTINUE
      WRITE (ISTDE,*) 'Name of state'
      READ(*,'(A)') NAME
      K=INDEX(NAME,' ')
      IF (K == 1) THEN
         WRITE (ISTDE,*) 'Names may not start with a blank'
         GOTO 10
      ENDIF
      WRITE (ISTDE,*)
      WRITE (ISTDE,*) 'Mixing coefficients from a CI calc.?'
      YES = GETYN ()
      IF (YES) THEN
         NCI = 0
      ELSE
         NCI = 1
      ENDIF
      WRITE (ISTDE,*)
      WRITE (ISTDE,*) 'Calculate off-diagonal matrix elements? '
      YES = GETYN ()
      IF (YES) THEN
         NOFFD = 0
      ELSE
         NOFFD = 1
      ENDIF
!
!   Check compatibility of plant substitutions
!
!     CALL CHKPLT
!
!   Determine if there is to be any debug printout; this will be
!   made on the  .dbg  file
!
      CALL SETDBG
!
!   Perform machine- and installation-dependent setup
!
      CALL SETMC
!
!   Set up the physical constants
!
      CALL SETCON
!
!   Open the  .sum  file
!
      CALL SETSUM(NAME,NCI,NOFFD)
!
!   Open, check, load data from, and close, the  .csl  file
!
      CALL SETCSLA(NAME, NCORE_NOT_USED)
!
!   Get the remaining information
!
      CALL GETHFD(NAME)
!
!   Get the eigenvectors
!
!      WRITE(ISTDE,*) 'Block format?'
!      YES = GETYN ()
!      WRITE (ISTDE,*)
!      IF (YES) THEN
         CALL GETMIXBLOCK(NAME,NCI)
!      ELSE
!         IF (NCI == 0) THEN
!            CALL GETMIXC(NAME)
!         ELSE
!            CALL GETMIXA(NAME)
!         ENDIF
!      ENDIF
!
!   Append a summary of the inputs to the  .sum  file
!
      CALL STRSUM
!
!   Set up the table of logarithms of factorials
!
      CALL FACTT
!
!   Proceed with the HFSZEEMAN calculation
!
      CALL HFSZEEMAN(NOFFD)
!
!   Print completion message
!
      WRITE (ISTDE,*)
      WRITE (ISTDE,*) 'HFSZEEMAN95: Execution complete.'
!
      STOP
      END PROGRAM HFSZEEMAN06
