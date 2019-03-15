!***********************************************************************
!                                                                      *
      INTEGER FUNCTION IQ (ISUBSH, ICSF)
!                                                                      *
!   IQ is the occupation of subshell ISUBSH in CSF  ICSF.              *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 30 Oct 1992   *
!   Modified by G. Gaigalas                                 May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:48:38   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE parameter_def, ONLY: NNNW
      USE ORB_C,         ONLY: IQA, NCF
      USE IOUNIT_C,      ONLY: istde
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: ISUBSH
      INTEGER             :: ICSF
!-----------------------------------------------
!
      IQ = IQA(isubsh,icsf)
!
      RETURN
      END FUNCTION IQ
