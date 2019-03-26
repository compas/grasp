!***********************************************************************
!                                                                      *
      INTEGER FUNCTION IQR (ISUBSH, ICSF)
!                                                                      *
!   IQR is the occupation of subshell ISUBSH in CSF  ICSF.             *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
!   Modified by G. Gaigalas                                 May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:35:54   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ORB_C,         ONLY: IQA
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: ISUBSH
      INTEGER  :: ICSF
!-----------------------------------------------
!
!
            IQR = IQA(isubsh,icsf)
!
      RETURN
      END FUNCTION IQR
