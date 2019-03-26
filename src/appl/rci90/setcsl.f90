!***********************************************************************
!
      SUBROUTINE SETCSL(NAME, NCORE, NBLKIN, IDBLK)
!
!  A container which calls setcsll to open, read <name>.c file to get
!     nblock, ncfblk(), idblk(), ncf (it is ncftot here).
!  It then calls lib92/lodcsh to get
!     ncore, nelec, nw, np(), nak(), nkl(), nkj(), nh()
!  The file pointer points to the first CSL record after this routine.
!
!  Xinghong He 98-06-23
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE memory_man
      USE hblock_C
      USE def_C
      USE orb_C, ONLY: ncftot=>ncf
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
!      USE setcsll_I
      USE lodcsh_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NCORE
      INTEGER  :: NBLKIN
      CHARACTER  :: NAME*(*)
      CHARACTER(LEN=8), DIMENSION(*) :: IDBLK
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IQADUM
!-----------------------------------------------
!
!     POINTER (ncfblk, ncfblk(0:*))
!-----------------------------------------------------------------------

      CALL ALLOC (NCFBLK, 0,  NBLKIN, 'NCFBLK', 'SETCSL' )
      CALL SETCSLL (21, NAME, NBLKIN, NBLOCK, NCFBLK(1), NCFTOT, IDBLK)
      CALL RALLOC (NCFBLK, 0, NBLOCK, 'NCFBLK', 'SETCSL' )

      REWIND (21)
      READ (21, *)

      !..Load header of <name>.c file
      CALL LODCSH (21, NCORE)

      RETURN
      END SUBROUTINE SETCSL
