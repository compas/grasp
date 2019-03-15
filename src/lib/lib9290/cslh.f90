!***********************************************************************
!
      SUBROUTINE CSLH(NAME, NCORE, NBLKIN, IDBLK)
!  A container which calls setcsll to open, read <name> file to get
!     nblock, ncfblk(), idblk(), ncftot.
!  It then calls lib92/lodcsh to get
!     ncore, nelec, nw, np(), nak(), nkl(), nkj(), nh()
!  The file pointer points to the first CSL record after this routine.
!
!  Called by rcivu, mcpvu
!
!  Xinghong He 98-06-23
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:47:12   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE memory_man
      USE HBLOCK_C
      USE DEF_C
      USE ORB_C, NCFTOT=>NCF
!
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE setcsll_I
      USE lodcsh_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NCORE
      INTEGER  :: NBLKIN
      CHARACTER  :: NAME*(*)
      CHARACTER  :: IDBLK(*)*8
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IQADUM
      INTEGER :: IERR
!-----------------------------------------------
!

! node-0 does exactly the same as the serial code does

      CALL ALLOC (NCFBLK, NBLKIN + 1, 'NCFBLK', 'CSLH')

      CALL SETCSLL (21, NAME, NBLKIN, NBLOCK, NCFBLK, NCFTOT, IDBLK)
      CALL RALLOC (NCFBLK, NBLOCK + 1, 'NCFBLK', 'CSLH')
      REWIND (21)
      READ (21, *)
      !..Load header of <name> file
      CALL LODCSH (21, NCORE)

      RETURN
      END SUBROUTINE CSLH
