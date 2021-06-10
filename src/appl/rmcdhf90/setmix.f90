!***********************************************************************
      SUBROUTINE SETMIX(NAME)
!
!   Opens the  .mix  file on stream 25; writes a header to this file.  *
!                                                                      *
!   Call(s) to: [LIB92]: LENGTH, OPENFL.                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 24 Dec 1992   *
!   Modified by Xinghong He               Last revision: 13 Jul 1998   *
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE def_C
      USE foparm_C
      USE mcpa_C
      USE mcpB_C
      USE orb_C
      USE hblock_C
      USE iounit_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE openfl_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER  :: NAME*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IERR
!-----------------------------------------------
!     POINTER (PCCMIN,ICCMIN(1))
!     POINTER (PNTRIQ,RIQDUM)
!
!
!     !...Wants nblock only
!     POINTER (pncfblk,ncfblk(0:*))
!
!-----------------------------------------------------------------------
      CALL OPENFL (25, NAME, 'UNFORMATTED', 'NEW', IERR)
      IF (IERR /= 0) THEN
         WRITE (ISTDE, *) 'Error when opening ', NAME(1:LEN_TRIM(NAME))
         ERROR STOP
      ENDIF
!
!   Write the file header
!
      WRITE (25) 'G92MIX'
      WRITE (25) NELEC, NCF, NW, 0, 0, NBLOCK
!     ...The above record will be overidden in matrix.f
!        with the final form of
!     WRITE (25)  NELEC, NCF, NW, nvectot, nvecsiz, nblock

      RETURN
      END SUBROUTINE SETMIX
