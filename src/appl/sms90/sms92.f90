!***********************************************************************
!***********************************************************************
!***********************************************************************
!**                                                                  ***
!**                                                                  ***
!**            Relativistic Isotope Shift Program                    ***
!**            GRASP92 Version Dynamic Storage                       ***
!**                                                                  ***
!**   ============================================================   ***
!**   Copyright (c) 1995 by P Jonsson, F A Parpia, and C F Fischer   ***
!**   ============================================================   ***
!**   All rights reserved.  No part of this software or its accom-   ***
!**   panying documentation may be reproduced, stored in a retrie-   ***
!**   val system,  or transmitted,  in any form or  by any  means,   ***
!**   electronic, mechanical,  photocopying,  recording, or other-   ***
!**   wise, without the prior written permission of the authors.     ***
!**                                                                  ***
!**                           Disclaimer                             ***
!**                           ==========                             ***
!**   The  authors make  no warranties,  express or implied,  that   ***
!**   this software or its  accompanying documentation are free of   ***
!**   error or that  they will meet your requirements for any par-   ***
!**   ticular application.  Neither this software nor its accompa-   ***
!**   nying documentation should be relied  upon for solving prob-   ***
!**   lems if an incorrect solution could result in injury or loss   ***
!**   of, or damage to, property. If you use this software or  its   ***
!**   accompanying documentation,  you do so entirely  at your own   ***
!**   risk;  the authors disclaim all liability for direct or con-   ***
!**   sequential damage.                                             ***
!**                                                                  ***
!***********************************************************************
!***********************************************************************
!***********************************************************************
!                                                                      *
      PROGRAM SMS92
!                                                                      *
!   Entry routine for SMS92. Controls the entire computation.          *
!                                                                      *
!   Call(s) to: [LIB92]: GETMIX, SETCSL, SETMC, SETCON.                *
!               [SMS92]: CHKPLT, GETSMD, SETDBG, SETSUM, SMS           *
!                        STRSUM.                                       *
!               [NJGRAF]: FACTT.                                       *
!                                                                      *
!   Written by Per Jonsson                Last revision: 17 Jan 1996   *
!   Modified by Gediminas Gaigalas for new spin-angular integration.   *
!                                         Last revision:    Nov 2017   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:07:11   1/ 3/07
!...Switches:
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I
      USE setdbg_I
      USE setmc_I
      USE setcon_I
      USE setsum_I
      USE setcsla_I
      USE getsmd_I
      USE getmixblock_I
      USE strsum_I
      USE factt_I
      USE sms_I
      IMPLICIT NONE
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
!...  /DEFAULT/
      COMMON /DEFAULT/ NDEF
      INTEGER   NDEF
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, NCI, NCORE_NOT_USED
      LOGICAL :: YES
      CHARACTER :: NAME*24
!-----------------------------------------------
!
      WRITE (6, *)
      WRITE (6, *) 'RSMS'
      WRITE (6, *) 'This is the isotope shift program'
      WRITE (6, *) 'Input files: isodata, name.c, name.(c)m, name.w'
      WRITE (6, *) 'Output files: name.(c)i'

      WRITE (6, *)
      WRITE (6, *) 'Default settings?'
      YES = GETYN()
      WRITE (6, *)
      IF (YES) THEN
         NDEF = 0
      ELSE
         NDEF = 1
      ENDIF

   10 CONTINUE
      WRITE (6, *) 'Name of state'
      READ (*, '(A)') NAME
      K = INDEX(NAME,' ')
      IF (K == 1) THEN
         WRITE (6, *) 'Names may not start with a blank'
         GO TO 10
      ENDIF
      WRITE (6, *)
      WRITE (6, *) 'Mixing coefficients from a CI calc.?'
      YES = GETYN()
      IF (YES) THEN
         NCI = 0
      ELSE
         NCI = 1
      ENDIF
!
!   Check compatibility of plant substitutions
!
!GG      CALL CHKPLT
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
      CALL SETSUM (NAME, NCI)
!
!   Open, check, load data from, and close, the  .csl  file
!
      CALL SETCSLA (NAME, NCORE_NOT_USED)
!
!   Get the remaining information
!
      CALL GETSMD (NAME)
!
!   Get the eigenvectors
!
!GG      WRITE (6, *) 'Block format?'
!GG      YES = GETYN()
!GG      WRITE (6, *)
!GG      IF (YES) THEN
         CALL GETMIXBLOCK (NAME, NCI)
!GG      ELSE
!GG         IF (NCI == 0) THEN
!GG            CALL GETMIXC (NAME)
!GG         ELSE
!GG            CALL GETMIXA (NAME)
!GG         ENDIF
!GG      ENDIF
!
!   Append a summary of the inputs to the  .sum  file
!
      CALL STRSUM
!
!   Set up the table of logarithms of factorials
!
      CALL FACTT
!
!   Proceed with the SMS calculation
!
      CALL SMS
!
!   Print completion message
!
      WRITE (6, *)
      WRITE (6, *) 'RSMS: Execution complete.'
!
      STOP
      END PROGRAM SMS92
